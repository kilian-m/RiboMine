"""Ribo-seq QC and verdict for one BAM (RiboseQC-style).

From a genome BAM, the annotation index (`ribomine.qc.annotation`) and STAR's
Log.final.out, `qc()` measures:

  read length     share of reads in the footprint window (25-36 nt by default)
  periodicity     frame of the read 5' end within CDS, each read length shifted to
                  its own dominant frame and pooled: in-frame fraction (chance 1/3)
                  and total variation distance to uniform (0 = flat, 2/3 = one frame)
  regions         CDS / 5'UTR / 3'UTR / ncRNA / intron / intergenic / mito
  metagene        5'-end density around annotated start and stop codons
  top 5' locus    share of reads at the most covered 5' position
  mapping         STAR's unique / multimapping / unmapped fractions

and calls RIBO-SEQ, TI-SEQ or NOT RIBO-SEQ or LOW QUALITY, with reasons. The
thresholds are `qc.*` config keys; "genic" reads are all but intergenic and mito.

  periodic_strong  in-frame >= periodic_strong and TVD >= tvd_strong
  periodic_weak    in-frame >= periodic_min and TVD >= tvd_min
                   (both need >= min_cds_reads CDS reads)
  len_ok           modal length in footprint_len_lo..footprint_len_hi and
                   >= read_len_peak_frac_min of reads in that window
  cds_ok           CDS >= cds_enrich_min of genic reads
  strong_cds       CDS >= cds_strong of genic reads
  contaminant      top 5' locus > single_locus_max of reads
  locus_over_cds   top 5' locus >= locus_min_to_judge of reads and
                   >= locus_vs_cds x the CDS share of all reads
  low_unique       uniquely mapped fraction < min_unique_frac
  tiseq_like       start-codon peak >= tiseq_min_peak reads and
                   >= tiseq_ratio_min x the mean of the CDS body

  translating = cds_ok and not contaminant and not locus_over_cds and
                (periodic_strong or (len_ok and (periodic_weak or strong_cds)))
  usable      = translating and not low_unique and two hard gates on all reads
                that periodicity cannot excuse: CDS >= cds_region_min of all reads
                and >= read_len_min_frac of reads in the footprint window
  verdict     = TI-SEQ if usable and tiseq_like, RIBO-SEQ if usable, otherwise
                NOT RIBO-SEQ or LOW QUALITY

The verdict favours specificity: a real but messy library can be refused. It is
measured on a local alignment of untrimmed reads, because soft-clipping the 5'
construct gives a cleaner footprint 5' end than end-to-end alignment.
"""
from __future__ import annotations

import logging
import os
import pickle
from collections import Counter, defaultdict

import numpy as np
import pysam
from ncls import NCLS

from ..config import Config
from ..process.star import parse_log

LOG = logging.getLogger("ribomine.qc.verdict")

# Metagene window; a property of the measurement, not a calling threshold.
METAGENE_WIN = (-30, 60)          # nt around the start/stop codon, 5'-end based
# Mito share of reads above which the library is flagged as mitoribosome profiling.
MITO_DOMINANT = 0.30
# Antisense deposit: fewer than this fraction of the reads on a CDS are on its strand
# (ribo-seq libraries measure 98-100 % sense). Judged only when at least
# ANTISENSE_MIN_READS reads lie on a CDS.
ANTISENSE_SENSE_FRAC = 0.50
ANTISENSE_MIN_READS = 200


# --------------------------------------------------------------------------
def _load_index(path: str):
    """Load the annotation index. Returns (index, NCLS interval indexes per category
    and chromosome, sorted codon positions per (chrom, strand))."""
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"annotation index not found: {path} (build it with ribomine.qc.annotation)"
        )
    with open(path, "rb") as fh:
        idx = pickle.load(fh)
    ncls: dict[str, dict[str, NCLS]] = {}
    for cat in ("cds", "utr5", "utr3", "exon_nc", "gene"):
        ncls[cat] = {}
        for chrom, rec in idx[cat].items():
            s, e = rec["start"], rec["end"]
            if len(s):
                ncls[cat][chrom] = NCLS(s, e, np.arange(len(s), dtype=np.int64))
    # sorted start/stop codon positions per (chrom, strand) for the metagene
    codons: dict[str, dict[tuple, np.ndarray]] = {}
    for kind in ("start_codon", "stop_codon"):
        codons[kind] = {}
        for chrom, rec in idx[kind].items():
            pos, strand = rec["pos"], rec["strand"]
            for sv in (1, -1):
                m = strand == sv
                if m.any():
                    codons[kind][(chrom, sv)] = np.sort(pos[m])
    return idx, ncls, codons


def _overlap_ids(ncls_c, pos: int) -> list[int]:
    """ids of intervals covering the single base `pos` (NCLS is half-open)."""
    if ncls_c is None:
        return []
    return [i for _, _, i in ncls_c.find_overlap(pos, pos + 1)]


def _classify(idx, ncls, chrom: str, pos: int, strand: int):
    """(region, CDS frame or None, CDS at this base on either strand?) for a position.

    `strand` is the read strand (+1/-1); a CDS call requires a strand match. The
    third value lets `qc` count reads antisense to a CDS, which are otherwise
    classified as another region (mostly intron).
    """
    # The mitochondrion is its own region, so its reads never enter the nuclear
    # periodicity. The frame is still returned for mito CDS, which the caller pools
    # separately (mitoribosome profiling).
    mito = chrom in ("MT", "Mt", "chrM", "M")
    cds = ncls["cds"].get(chrom)
    ids = _overlap_ids(cds, pos)
    on_cds = bool(ids)                 # a CDS at this base, on either strand
    if ids:
        rec = idx["cds"][chrom]
        for i in ids:
            if rec["strand"][i] == strand:
                cs, ce, fr = rec["start"][i], rec["end"][i], int(rec["frame"][i])
                # GTF `frame` = bases from the feature start to the next codon start,
                # so the codon phase of x is (x - start - frame) % 3 (0 = first base)
                if strand == 1:
                    frame = (pos - cs - fr) % 3
                else:
                    frame = (ce - 1 - pos - fr) % 3
                return ("mito" if mito else "CDS"), frame, on_cds
    if mito:
        return "mito", None, on_cds     # on MT, outside its CDS
    for cat, tag in (("utr5", "5'UTR"), ("utr3", "3'UTR")):
        c = ncls[cat].get(chrom)
        if _overlap_ids(c, pos):
            return tag, None, on_cds
    if _overlap_ids(ncls["exon_nc"].get(chrom), pos):
        return "ncRNA", None, on_cds
    if _overlap_ids(ncls["gene"].get(chrom), pos):
        return "intron", None, on_cds
    return "intergenic", None, on_cds


def _read_5p(aln):
    """(chrom, 5'-genomic-position, strand, aligned-length) for a read."""
    if aln.is_reverse:
        return aln.reference_name, aln.reference_end - 1, -1, aln.query_alignment_length
    return aln.reference_name, aln.reference_start, 1, aln.query_alignment_length


def _tvd_uniform(counts3) -> float:
    """Total variation distance of a 3-frame distribution to uniform (1/3 each).
    0 = flat (no periodicity); 2/3 = all mass in one frame (perfect periodicity)."""
    t = sum(counts3)
    if not t:
        return 0.0
    p = [c / t for c in counts3]
    return 0.5 * sum(abs(pi - 1 / 3) for pi in p)


def periodicity(bam: str, index_path: str, *, max_reads: int = 200_000,
                label: str = "") -> dict:
    """3-nt periodicity of the final (deliverable) BAM. Reported only.

    `qc()` decides the verdict on a sample of untrimmed, locally aligned reads;
    this measures the trimmed, filtered reads of the final alignment. The two can
    differ, e.g. after a wrong trim.

    About `max_reads` reads are scored, taken at a fixed stride through the file,
    because a coordinate-sorted BAM starts with one chromosome. `max_reads <= 0`
    scores every read.
    """
    idx, ncls, _ = _load_index(index_path)
    bamf = pysam.AlignmentFile(bam, "rb")
    # The index counts alignments, not reads (a kept multimapper counts several
    # times); it is used only to set the stride.
    try:
        n_alignments = bamf.mapped
    except ValueError:                   # no index -- score everything
        n_alignments = 0
    stride = max(1, n_alignments // max_reads) if (n_alignments and max_reads > 0) else 1

    len_hist: Counter = Counter()
    region_hist: Counter = Counter()
    frame_by_len: dict[int, list[int]] = defaultdict(lambda: [0, 0, 0])
    seen = n = softclip = 0
    for aln in bamf:
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        seen += 1
        if (seen - 1) % stride:
            continue
        n += 1
        chrom, p5, strand, L = _read_5p(aln)
        len_hist[L] += 1
        ct = aln.cigartuples
        softclip += sum(ln for op, ln in ct if op == 4) if ct else 0
        region, frame, _ = _classify(idx, ncls, chrom, p5, strand)
        region_hist[region] += 1
        if region == "CDS" and frame is not None:
            frame_by_len[L][frame] += 1
    bamf.close()

    if n == 0:
        raise ValueError(f"no usable alignments in {bam}")

    # as in qc(): shift each read length to its own dominant frame before pooling,
    # so lengths whose 5' ends differ in phase do not cancel
    cds_reads = sum(sum(v) for v in frame_by_len.values())
    inframe = 0
    pooled = [0, 0, 0]
    for c3 in frame_by_len.values():
        if not sum(c3):
            continue
        dom = int(np.argmax(c3))
        inframe += c3[dom]
        for f in range(3):
            pooled[(f - dom) % 3] += c3[f]

    genic = n - region_hist.get("intergenic", 0) - region_hist.get("mito", 0)
    res = {
        "n_reads_in_bam": seen,          # primary alignments, counted
        "n_reads_scored": n,
        "stride": stride,
        # Aligned length (soft clips excluded). The mode (read_len_mode) is the
        # footprint peak and is robust to contamination in the tails.
        "mean_mapped_len": round(sum(L * c for L, c in len_hist.items()) / n, 1),
        # Mean soft-clipped bases per mapped read: near 0 after a clean trim, large
        # when residual construct was clipped by the aligner instead of trimmed.
        "mean_mapped_softclip": round(softclip / n, 2),
        "read_len_mode": max(len_hist, key=lambda k: len_hist[k]),
        "periodicity_inframe_frac": round(inframe / cds_reads, 3) if cds_reads else 0.0,
        "periodicity_tvd_uniform": round(_tvd_uniform(pooled), 3),
        "n_cds_reads": cds_reads,
        "cds_frac_of_genic": round(region_hist.get("CDS", 0) / genic, 3) if genic else 0.0,
        "region_frac": {k: round(v / n, 4)
                        for k, v in sorted(region_hist.items(), key=lambda x: -x[1])},
    }
    LOG.info("%s: BAM %s reads, in-frame %.0f%% (TVD %.2f) on %s CDS reads of %s scored",
             label or os.path.basename(bam), f"{res['n_reads_in_bam']:,}",
             100 * res["periodicity_inframe_frac"], res["periodicity_tvd_uniform"],
             f"{cds_reads:,}", f"{n:,}")
    return res


def qc(bam: str, index_path: str, cfg: Config, *, star_log: str = "", label: str = "",
       contam: dict | None = None, pileup: dict | None = None,
       total_reads: int | None = None) -> dict:
    """The ribo-seq verdict for one locally aligned BAM (see the module docstring).

    Scores the first `qc.max_reads_scored` primary alignments. `contam` / `pileup`
    are the stats dicts of the two upstream filters; `total_reads` is the read
    count of the whole run, used to project its number of usable reads.
    """
    # ---- thresholds (tuned on a 187-sample human cohort; see config.py)
    fp_len_lo = cfg["qc.footprint_len_lo"]
    fp_len_hi = cfg["qc.footprint_len_hi"]
    peak_frac_min = cfg["qc.read_len_peak_frac_min"]
    read_len_min_frac = cfg["qc.read_len_min_frac"]  # hard gate: reads in the footprint window
    cds_enrich_min = cfg["qc.cds_enrich_min"]
    cds_region_min = cfg["qc.cds_region_min"]    # hard gate: CDS share of all reads
    cds_strong_min = cfg["qc.cds_strong"]        # strong CDS enrichment excuses weak periodicity
    periodic_min = cfg["qc.periodic_min"]        # in-frame fraction (chance 0.33)
    periodic_strong_min = cfg["qc.periodic_strong"]
    tvd_min = cfg["qc.tvd_min"]
    tvd_strong_min = cfg["qc.tvd_strong"]
    min_cds_reads = cfg["qc.min_cds_reads"]      # need this many CDS reads to judge periodicity
    single_locus_max = cfg["qc.single_locus_max"]
    locus_vs_cds = cfg["qc.locus_vs_cds"]
    locus_min_to_judge = cfg["qc.locus_min_to_judge"]
    min_unique_frac = cfg["qc.min_unique_frac"]
    tiseq_ratio_min = cfg["qc.tiseq_ratio_min"]
    tiseq_min_peak = cfg["qc.tiseq_min_peak"]
    max_reads = cfg["qc.max_reads_scored"]
    # Strong periodicity may excuse an atypical footprint length but not an empty CDS:
    # periodicity is measured on CDS reads alone, so a few genuine footprints in an
    # otherwise non-coding library still score as strongly periodic.
    cds_floor = cds_enrich_min

    idx, ncls, codons = _load_index(index_path)
    bamf = pysam.AlignmentFile(bam, "rb")

    len_hist: Counter = Counter()
    region_hist: Counter = Counter()
    # per-length frame counts within CDS (5'-end frame)
    frame_by_len: dict[int, list[int]] = defaultdict(lambda: [0, 0, 0])
    # the same for the mitochondrial CDS, pooled separately; not used in the verdict
    mito_frame_by_len: dict[int, list[int]] = defaultdict(lambda: [0, 0, 0])
    # metagene: 5'-end offset (read - codon, translation dir) -> count
    meta = {"start_codon": Counter(), "stop_codon": Counter()}
    loc_hist: Counter = Counter()
    n = 0
    # reads that land on a coding base, split by whether they are on the gene's strand
    n_cds_sense = n_cds_anti = 0

    for aln in bamf:
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        if n >= max_reads:
            break
        n += 1
        chrom, p5, strand, L = _read_5p(aln)
        len_hist[L] += 1
        # 5' anchor including the soft clip (`pos - clip5`), as in fqdissect.profile:
        # unchanged when a chance match extends the alignment into the clip
        anchor = (aln.reference_end + (aln.cigartuples[-1][1] if aln.cigartuples[-1][0] == 4 else 0)) \
            if aln.is_reverse else \
            (aln.reference_start - (aln.cigartuples[0][1] if aln.cigartuples[0][0] == 4 else 0))
        loc_hist[(chrom, aln.is_reverse, anchor)] += 1

        region, frame, on_cds = _classify(idx, ncls, chrom, p5, strand)
        region_hist[region] += 1
        if region == "CDS" and frame is not None:
            frame_by_len[L][frame] += 1
        elif region == "mito" and frame is not None:
            mito_frame_by_len[L][frame] += 1
        # sense or antisense to the CDS at this base (see `_classify`)
        if on_cds and region != "mito":
            if region == "CDS":
                n_cds_sense += 1
            else:
                n_cds_anti += 1

        for kind in ("start_codon", "stop_codon"):
            arr = codons[kind].get((chrom, strand))
            if arr is None or not len(arr):
                continue
            j = np.searchsorted(arr, p5)
            for cand in (j - 1, j):
                if 0 <= cand < len(arr):
                    d = (p5 - arr[cand]) * strand   # translation-direction offset
                    if METAGENE_WIN[0] <= d <= METAGENE_WIN[1]:
                        meta[kind][int(d)] += 1
                        break

    if n == 0:
        raise ValueError(f"no usable alignments in {bam}")

    # ---- periodicity: align each length to its own dominant frame, then pool
    cds_reads = sum(sum(v) for v in frame_by_len.values())
    inframe = 0
    pooled = [0, 0, 0]
    per_len = {}
    for L, c3 in frame_by_len.items():
        tot = sum(c3)
        if tot == 0:
            continue
        dom = int(np.argmax(c3))
        inframe += c3[dom]
        for f in range(3):
            pooled[(f - dom) % 3] += c3[f]   # shift so dominant frame -> 0
        if tot >= 30:
            per_len[L] = {"n": tot, "dom_frame": dom, "inframe_frac": c3[dom] / tot}
    periodicity_frac = inframe / cds_reads if cds_reads else 0.0
    tvd = _tvd_uniform(pooled)

    # ---- the same on the mitochondrial CDS; reported only, to recognise
    # mitoribosome profiling
    mito_cds_reads = sum(sum(v) for v in mito_frame_by_len.values())
    mito_inframe = 0
    mito_pooled = [0, 0, 0]
    for L, c3 in mito_frame_by_len.items():
        if sum(c3) == 0:
            continue
        dom = int(np.argmax(c3))
        mito_inframe += c3[dom]
        for f in range(3):
            mito_pooled[(f - dom) % 3] += c3[f]
    mito_periodicity = mito_inframe / mito_cds_reads if mito_cds_reads else 0.0
    mito_tvd = _tvd_uniform(mito_pooled)

    # ---- region fractions (genic = all but intergenic and mito, for CDS enrichment)
    reg_frac = {k: v / n for k, v in region_hist.items()}
    genic = n - region_hist.get("intergenic", 0) - region_hist.get("mito", 0)
    cds_of_genic = region_hist.get("CDS", 0) / genic if genic else 0.0

    top5p_locus_frac = loc_hist.most_common(1)[0][1] / n if loc_hist else 0.0

    # ---- strand of the reads on coding genes
    n_on_cds = n_cds_sense + n_cds_anti
    cds_sense_frac = n_cds_sense / n_on_cds if n_on_cds else 0.0
    # Antisense-dominated (reverse-complemented or reverse-stranded) deposit: only the
    # sense strand is scored, so its CDS enrichment is near zero. Flagged so that the
    # verdict reason can name the strand.
    antisense_deposit = (n_on_cds >= ANTISENSE_MIN_READS
                         and cds_sense_frac < ANTISENSE_SENSE_FRAC)

    # ---- footprint length peak
    mode_len = max(len_hist, key=lambda k: len_hist[k])
    peak_frac = sum(len_hist[ln] for ln in range(fp_len_lo, fp_len_hi + 1)) / n

    mapping = {}
    if not star_log:
        cand = os.path.join(os.path.dirname(os.path.abspath(bam)), "Log.final.out")
        star_log = cand if os.path.exists(cand) else ""
    if star_log and os.path.exists(star_log):
        mapping = parse_log(star_log) or {}
    frac_unique = mapping.get("frac_unique")

    # ---- signals
    len_ok = fp_len_lo <= mode_len <= fp_len_hi and peak_frac >= peak_frac_min
    cds_ok = cds_of_genic >= cds_enrich_min
    enough_cds = cds_reads >= min_cds_reads
    # 3-nt periodicity, the hallmark of translating ribosomes
    periodic_strong = enough_cds and periodicity_frac >= periodic_strong_min and tvd >= tvd_strong_min
    periodic_weak = enough_cds and periodicity_frac >= periodic_min and tvd >= tvd_min
    strong_cds = cds_of_genic >= cds_strong_min
    # a majority at one locus = contaminant, not footprints
    contaminant = top5p_locus_frac > single_locus_max
    # poorly mapping library, unusable even if the mappable minority is periodic
    # (cohort median 42 % unique, 10th percentile 20 %)
    low_unique = frac_unique is not None and frac_unique < min_unique_frac
    cds_floor_ok = cds_of_genic >= cds_floor
    # One 5' position holding as many reads as the whole CDS: contaminant-dominated
    # even below single_locus_max. Relative to the library's own CDS share, so no
    # absolute cut-off (good libraries: top locus ~100x below the CDS share).
    cds_share = reg_frac.get("CDS", 0.0)
    locus_over_cds = (top5p_locus_frac >= locus_vs_cds * cds_share
                      and top5p_locus_frac >= locus_min_to_judge)

    # TI-seq: initiation inhibitors (harringtonine/LTM) hold ribosomes at start codons
    # while elongating ones run off. Ratio = start-codon peak (5' offset -15..-9) over
    # the mean of in-frame positions 5-15 codons downstream (floored at 1 read);
    # elongating ribo-seq reaches ~30.
    msc = meta["start_codon"]
    peak_off = max(range(-15, -8), key=lambda k: msc.get(k, 0)) if msc else -12
    peak_h = msc.get(peak_off, 0)
    body = [msc.get(peak_off + 3 * i, 0) for i in range(5, 16)
            if peak_off + 3 * i <= METAGENE_WIN[1]]
    body_mean = float(np.mean(body)) if body else 0.0
    start_ratio = peak_h / max(body_mean, 1.0)
    tiseq_like = peak_h >= tiseq_min_peak and start_ratio >= tiseq_ratio_min

    # ---- hard gates, measured on all reads. Unlike `cds_ok` (a fraction of genic
    # reads) and `len_ok`, strong periodicity does not excuse them.
    cds_region_ok = cds_share >= cds_region_min          # CDS share of all reads
    read_len_dominant = peak_frac >= read_len_min_frac   # share of reads in the footprint window

    # ---- verdict (see the module docstring)
    translating = (not contaminant) and cds_floor_ok and (not locus_over_cds) and (
        periodic_strong or (len_ok and cds_ok and (periodic_weak or strong_cds)))
    usable = translating and not low_unique and cds_region_ok and read_len_dominant
    if not usable:
        verdict = "NOT RIBO-SEQ or LOW QUALITY"
    elif tiseq_like:
        verdict = "TI-SEQ"
    else:
        verdict = "RIBO-SEQ"
    is_riboseq = usable                     # TI-seq counts as usable translation data

    reasons = [
        f"footprint length {'ok' if len_ok else 'ATYPICAL'} "
        f"(mode {mode_len} nt, {peak_frac:.0%} in {fp_len_lo}-{fp_len_hi})",
        f"read-length window {'ok' if read_len_dominant else 'LOW'} "
        f"({peak_frac:.0%} of reads in {fp_len_lo}-{fp_len_hi} nt, floor {read_len_min_frac:.0%})",
        f"CDS enrichment {'ok' if cds_ok else 'LOW'} ({cds_of_genic:.0%} of genic reads in CDS)",
        f"region composition {'ok' if cds_region_ok else 'LOW'} "
        f"({cds_share:.0%} of all reads in CDS, floor {cds_region_min:.0%})",
    ]
    if enough_cds:
        lvl = "STRONG" if periodic_strong else ("weak" if periodic_weak else "ABSENT")
        reasons.append(f"periodicity {lvl} (in-frame {periodicity_frac:.0%}, TVD {tvd:.2f})")
    else:
        reasons.append(f"periodicity untested (only {cds_reads} CDS reads)")
    if reg_frac.get("mito", 0.0) >= MITO_DOMINANT:
        # the signals above come from the nuclear reads only; say so
        reasons.append(
            f"MITOCHONDRIAL-DOMINATED ({reg_frac['mito']:.0%} of reads) -- this looks like "
            f"mitoribosome profiling. The verdict above was decided on the {cds_reads:,} "
            f"NUCLEAR CDS reads; the {mito_cds_reads:,} MT-CDS reads are "
            f"{mito_periodicity:.0%} in-frame (TVD {mito_tvd:.2f}) and are not scored")
    if antisense_deposit:
        reasons.append(
            f"ANTISENSE to the genes -- only {cds_sense_frac:.0%} of the {n_on_cds:,} reads "
            f"on a CDS are on its strand. This is a reverse-stranded library; the usual "
            f"cause is stranded RNA-seq (dUTP), not ribo-seq, whose footprints are pieces of "
            f"the mRNA and therefore sense. RiboMine scores the sense strand, so the CDS "
            f"enrichment and periodicity above are of the minority of reads that sit the "
            f"right way round -- read them as a floor, not a measurement")
    if frac_unique is not None:
        reasons.append(f"unique mapping {'LOW' if low_unique else 'ok'} ({frac_unique:.0%})")
    if contaminant:
        reasons.append(f"contaminant: {top5p_locus_frac:.0%} of reads at one locus")
    if locus_over_cds and not contaminant:
        reasons.append(f"contaminant-dominated: one locus holds {top5p_locus_frac:.0%} of reads, "
                       f"the whole CDS only {cds_share:.0%}")
    if not cds_floor_ok:
        reasons.append(f"CDS essentially empty ({cds_of_genic:.0%} of genic reads)")
    if tiseq_like:
        reasons.append(f"start-codon enrichment {start_ratio:.0f}x -> TI-seq")

    # headline reason for a NOT RIBO-SEQ or LOW QUALITY call
    verdict_reason = ""
    if not usable:
        why = []
        if contaminant:
            why.append(f"{top5p_locus_frac:.0%} of reads at one 5' locus (contaminant)")
        elif locus_over_cds:
            why.append(f"one 5' locus holds {top5p_locus_frac:.0%} of reads, as much as the "
                       f"entire CDS ({cds_share:.0%}) — contaminant-dominated")
        if not cds_floor_ok:
            why.append(f"only {cds_of_genic:.0%} of genic reads are in CDS (not translation)")
        # hard gates
        if not cds_region_ok:
            why.append(f"CDS is only {cds_share:.0%} of the region composition "
                       f"(need {cds_region_min:.0%})")
        if not read_len_dominant:
            why.append(f"only {peak_frac:.0%} of reads are footprint-length "
                       f"({fp_len_lo}-{fp_len_hi} nt; need {read_len_min_frac:.0%})")
        if antisense_deposit:
            msg = (f"the reads are ANTISENSE to the genes ({cds_sense_frac:.0%} of those on a "
                   f"CDS are on its strand), so the CDS is empty by construction")
            if len_ok:
                # footprint-length reads: the strand is the likely cause, so lead with it
                why.insert(0, msg + " -- a reverse-complemented deposit would look exactly "
                                     "like this, and its footprints ARE the right length")
            else:
                # not footprint-length either (typically reverse-stranded RNA-seq):
                # mention the strand, but not first
                why.append(msg + " -- but the reads are not footprint-length either, so the "
                                 "strand is not what disqualifies this run")
        if low_unique:
            why.append(f"only {frac_unique:.0%} of reads map uniquely (non-genomic / "
                       f"multimapping junk)")
        if not translating and not contaminant:
            if not (periodic_strong or periodic_weak):
                why.append("no 3-nt periodicity (not translating ribosomes)")
            elif not (len_ok and cds_ok):
                why.append("atypical footprint length and/or low CDS enrichment")
        verdict_reason = "; ".join(why) or "does not meet ribo-seq criteria"

    res = {
        "label": label or os.path.basename(bam),
        "n_reads_scored": n,
        "verdict": verdict,
        "verdict_reason": verdict_reason,
        "is_riboseq": is_riboseq,
        "reasons": reasons,
        "signals": {
            "footprint_len_ok": len_ok, "cds_enriched": cds_ok,
            "periodic_strong": periodic_strong, "periodic_weak": periodic_weak,
            "contaminant": contaminant, "low_unique_mapping": low_unique,
            "tiseq": tiseq_like,
            "locus_over_cds": locus_over_cds, "cds_floor_ok": cds_floor_ok,
            "antisense_deposit": bool(antisense_deposit),
            # hard gates (measured on all reads)
            "cds_region_ok": cds_region_ok, "read_len_dominant": read_len_dominant,
        },
        "start_codon_ratio": round(start_ratio, 1),
        "read_len_mode": mode_len,
        "read_len_peak_frac": round(peak_frac, 3),
        "periodicity_inframe_frac": round(periodicity_frac, 3),
        "periodicity_tvd_uniform": round(tvd, 3),
        "n_cds_reads": cds_reads,
        # mitochondrial CDS, pooled separately; meaningful for mitoribosome profiling
        "mito_periodicity_inframe_frac": round(mito_periodicity, 3),
        "mito_periodicity_tvd_uniform": round(mito_tvd, 3),
        "n_mito_cds_reads": mito_cds_reads,
        "mito_dominated": bool(reg_frac.get("mito", 0.0) >= MITO_DOMINANT),
        "cds_frac_of_genic": round(cds_of_genic, 3),
        # strandedness: a ribo-seq library is ~100% sense over CDS
        "cds_sense_frac": round(cds_sense_frac, 3),
        "n_reads_on_cds": n_on_cds,
        "antisense_deposit": bool(antisense_deposit),
        "top5p_locus_frac": round(top5p_locus_frac, 3),
        "region_frac": {k: round(v, 4) for k, v in sorted(reg_frac.items(), key=lambda x: -x[1])},
        "read_len_hist": dict(sorted(len_hist.items())),
        "frame_by_len": {int(L): v["inframe_frac"] for L, v in sorted(per_len.items())},
        "frame_by_len_dom": {int(L): v["dom_frame"] for L, v in sorted(per_len.items())},
        "frame_by_len_n": {int(L): v["n"] for L, v in sorted(per_len.items())},
        "pooled_frame_counts": pooled,
        "metagene_start": {int(k): v for k, v in sorted(meta["start_codon"].items())},
        "metagene_stop": {int(k): v for k, v in sorted(meta["stop_codon"].items())},
        "mapping": mapping,
        # The thresholds this verdict was decided with, named as in config.py. The QC
        # figure draws them, and a stored qc.json stays interpretable after the
        # config changes.
        "thresholds": {
            "footprint_len_lo": fp_len_lo,
            "footprint_len_hi": fp_len_hi,
            "read_len_peak_frac_min": peak_frac_min,
            "read_len_min_frac": read_len_min_frac,
            "cds_enrich_min": cds_enrich_min,
            "cds_region_min": cds_region_min,
            "cds_strong": cds_strong_min,
            "periodic_min": periodic_min,
            "periodic_strong": periodic_strong_min,
            "tvd_min": tvd_min,
            "tvd_strong": tvd_strong_min,
            "min_cds_reads": min_cds_reads,
            "single_locus_max": single_locus_max,
            "locus_vs_cds": locus_vs_cds,
            "locus_min_to_judge": locus_min_to_judge,
            "min_unique_frac": min_unique_frac,
            "tiseq_ratio_min": tiseq_ratio_min,
            "tiseq_min_peak": tiseq_min_peak,
        },
    }

    if contam:
        c = {
            "frac_rRNA_tRNA_etc": contam.get("frac_contaminant_structured_rna"),
            "frac_low_complexity": contam.get("frac_low_complexity"),
            "n_sampled": contam.get("n_input"),
        }
        # pile-ups are removed after mapping; express them as a fraction of the
        # sampled reads, like the other two
        if pileup:
            n_samp = contam.get("n_input") or 1
            c["frac_position_pileup"] = round(pileup.get("n_reads_removed", 0) / n_samp, 4)
            c["top_pileups"] = pileup.get("top_pileups", [])[:3]
            c["frac_kept"] = round(
                (contam.get("n_kept", 0) - pileup.get("n_reads_removed", 0)) / n_samp, 4
            )
        else:
            c["frac_kept"] = contam.get("frac_kept")
        res["contaminants"] = c

    # ---- projected usable reads: n_reads_scored (mapped, contaminant- and pile-up-
    #      filtered) as a fraction of the sample, scaled to the whole dataset
    n_sampled = res.get("contaminants", {}).get("n_sampled")
    if total_reads and n_sampled:
        frac = res["n_reads_scored"] / n_sampled
        res["usable"] = {
            "total_dataset_reads": total_reads,
            "n_sampled": n_sampled,
            "usable_frac_of_sample": round(frac, 4),
            "projected_usable_reads": int(round(frac * total_reads)),
        }

    LOG.info(
        "%s: %s [len mode %dnt, in-frame %.0f%%, TVD %.2f, CDS %.0f%% of genic, "
        "start-ratio %.0fx]",
        res["label"], verdict, mode_len, 100 * periodicity_frac, tvd,
        100 * cds_of_genic, start_ratio,
    )
    if antisense_deposit:
        LOG.warning(
            "%s: ANTISENSE DEPOSIT -- only %.0f%% of the %s reads that touch a CDS are on "
            "its strand. The deposit looks reverse-complemented, and every number above is "
            "measured on the sense strand, so its CDS enrichment and periodicity are of the "
            "minority. Re-check before believing this run's verdict either way.",
            res["label"], 100 * cds_sense_frac, f"{n_on_cds:,}")
    if verdict_reason:
        LOG.info("   reason: %s", verdict_reason)
    for r in reasons:
        LOG.debug("   - %s", r)
    return res
