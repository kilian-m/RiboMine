"""Is this a ribo-seq library, and how good is it? -- RiboseQC-style QC.

Given a genome BAM (+ the cached annotation index from `ribomine.qc.annotation`,
+ the STAR Log.final.out for mapping stats), compute the signatures that define a
ribosome-profiling dataset:

  * read-length distribution   -- footprints peak tightly at ~28-32 nt
  * 3-nt periodicity           -- the defining hallmark: P-site (here the read
                                  5' end) frames cluster in one frame within CDS.
                                  Reported as the in-frame fraction and the total
                                  variation distance of the frame distribution to
                                  uniform (0 = flat/no periodicity, 0.67 = perfect)
  * region composition         -- CDS / 5'UTR / 3'UTR / ncRNA / intron / intergenic
                                  / mito; ribo-seq is strongly CDS-enriched
  * start / stop metagene      -- 5'-end density around annotated start/stop codons
  * top mapping locus          -- a single over-represented species (contaminant)
  * mapping stats              -- input / uniquely mapped / multimapping / unmapped

then a three-way verdict -- RIBO-SEQ / TI-SEQ / NOT RIBO-SEQ or LOW QUALITY --
with the reasons stated. The verdict is deliberately INCLUSIVE so that ribo-seq
variants are not thrown away: strong periodicity is decisive on its own (it holds
for TI-seq / QTI-seq, whose initiating ribosomes are still in-frame, and for
atypically short or long footprints); failing that, footprint-like + CDS-enriched
reads with any periodicity (or strong CDS enrichment) pass. Only a
contaminant-dominated or badly-mapping library is refused outright.

The QC is measured on the LOCAL alignment because soft-clipping the 5' construct
/ RT base gives a cleaner footprint 5' end -- and hence sharper periodicity --
than an end-to-end alignment of the trimmed reads.

Every threshold comes from the config (`qc.*`); nothing is hard-coded here.
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

# The window the metagene is computed over: a property of the plot/measurement,
# not a calling threshold, so it is not a config key.
METAGENE_WIN = (-30, 60)          # nt around the start/stop codon, 5'-end based
# above this share of reads on the mitochondrion, the library is not a cytosolic
# ribo-seq with mito contamination -- it is (or is dominated by) mitoribosome profiling
MITO_DOMINANT = 0.30
# A ribo-seq footprint is a piece of the mRNA, so it is SENSE to the gene: a real library
# is ~100% sense over CDS (measured: 98-100% across a 100-run cohort). Below this, the
# deposit is reverse-complemented, and RiboMine -- which scores the sense strand -- would
# otherwise refuse it for "no CDS enrichment" and never say why. Judged only when enough
# reads actually touch a CDS to make the fraction mean anything.
ANTISENSE_SENSE_FRAC = 0.50
ANTISENSE_MIN_READS = 200


# --------------------------------------------------------------------------
def _load_index(path: str):
    """The pickled annotation index -> NCLS interval indexes + sorted codon arrays."""
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
    """(region, CDS frame, is a CDS here at all?) for a genomic position.

    strand is the read strand (+1/-1); CDS/UTR require a strand match (ribo-seq
    footprints map sense to the transcript).

    The third value is what makes an ANTISENSE deposit visible. A read that lands on a
    CDS but on the wrong strand fails that match, falls through every other test, and
    comes out as "intron" -- because a CDS sits inside a gene. So a reverse-complemented
    deposit does not look like a strand problem: it looks like a library with no CDS
    enrichment and a great deal of intronic signal, and it is REJECTED for exactly that,
    with nothing said about the strand. Reporting "a CDS was here, whatever the strand"
    costs nothing (the overlap has already been computed) and lets `qc` count the reads
    on each side of it.
    """
    # The mitochondrion is kept as its OWN region, never folded into CDS: in an ordinary
    # cytosolic library its reads are mostly degradation background, and pooling them
    # would dilute the nuclear periodicity that the verdict rests on.
    #
    # But it still has 13 protein-coding genes and a reading frame, and MITORIBOSOME
    # profiling exists -- SRR28710935 ("Monitoring mitochondrial translation") is 60%
    # mito, and its 4,490 MT-CDS reads are 50% in-frame. Returning frame=None there threw
    # away the entire experiment and left the verdict to be decided by 1,699 nuclear reads
    # of contamination. So: report the frame, and let the caller keep it in a separate
    # pool (it does) rather than pretend the mitochondrion does not translate.
    mito = chrom in ("MT", "Mt", "chrM", "M")
    cds = ncls["cds"].get(chrom)
    ids = _overlap_ids(cds, pos)
    on_cds = bool(ids)                 # a CDS is at this base, on ONE strand or the other
    if ids:
        rec = idx["cds"][chrom]
        for i in ids:
            if rec["strand"][i] == strand:
                cs, ce, fr = rec["start"][i], rec["end"][i], int(rec["frame"][i])
                # Ensembl `frame` = bases to remove from the feature start to reach the
                # first base of the next codon, so the codon-phase of position x is
                # (x - featurestart - frame) % 3 (0 = first base of a codon). The
                # annotation hands us the phase directly -- no per-transcript walk.
                if strand == 1:
                    frame = (pos - cs - fr) % 3
                else:
                    frame = (ce - 1 - pos - fr) % 3
                return ("mito" if mito else "CDS"), frame, on_cds
    if mito:
        return "mito", None, on_cds     # on MT, but not inside one of its 13 CDS
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
    """3-nt periodicity of a FINISHED BAM -- measured on the deliverable itself.

    The verdict in `qc()` is decided on the QC stage's own alignment: a 200k-read
    sample, untrimmed, aligned locally. That is the right measurement to CALL a
    library with, but it is not a measurement of what came out of the pipeline. This
    is: the trimmed reads, end-to-end aligned, contaminant- and pile-up-filtered,
    deduplicated if that was asked for. It decides nothing -- it is the number you
    look at to see whether the reads you are about to analyse are periodic.

    The two can legitimately disagree (a trim that cut the footprint boundary wrong
    shows up here and nowhere else), which is the whole reason for measuring twice.

    Reads are taken at a fixed STRIDE through the file, not from its front: the BAM
    is coordinate-sorted by now, so its first `max_reads` reads are the first
    chromosome, and periodicity there is not periodicity everywhere. `max_reads <= 0`
    scores every read.
    """
    idx, ncls, _ = _load_index(index_path)
    bamf = pysam.AlignmentFile(bam, "rb")
    # The .bai's count is of ALIGNMENTS, not of reads -- with multimappers kept
    # (mapping.multimap_nmax > 1) one read contributes several. It is exactly what the
    # stride needs, and NOT what "reads in the BAM" means, so the reads are counted in
    # the pass below and this is used for nothing else.
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

    # same construction as qc(): align each read length to its own dominant frame
    # before pooling, so a length whose 5' end sits one base off does not cancel out
    # a length that is in phase
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
        "n_reads_in_bam": seen,          # primary alignments = reads, counted, not inferred
        "n_reads_scored": n,
        "stride": stride,
        # The ALIGNED length (soft clips excluded), which for the deliverable BAM is the
        # whole read -- it is aligned end-to-end, so nothing is clipped. It stops being
        # the whole read the moment someone sets mapping.align_ends_type to Local, and
        # then the aligned length is the honest one: it is the part that is genomic.
        #
        # Reported as BOTH a mean and a mode, because they answer different questions. The
        # mode is where the footprint PEAK is, and a peak is a mode: contamination in the
        # tails cannot move it (SRR25706716 peaks at 28 nt with a median of 26, dragged
        # down by miRNA). The mean is what you want when the whole distribution matters.
        "mean_mapped_len": round(sum(L * c for L, c in len_hist.items()) / n, 1),
        # Mean soft-clip per mapped read -- the non-genomic bases local alignment shaved
        # off the trimmed read. This is what makes the mean_footprint_len/mean_mapped_len
        # gap READABLE: a large gap can mean either residual construct clipped off the
        # mapped reads (a trim miss, e.g. SRR18113808 pre-fix ~9 nt) OR merely that longer
        # unmappable read-through inflates mean_footprint_len (an average over ALL input
        # reads) while the mapped footprints are clean. This tells the two apart -- near 0
        # is clean (a kept non-templated RT base contributes ~its penetrance); large is
        # residual the trim missed and STAR clipped rather than aligned.
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
    """The ribo-seq verdict for one locally-aligned BAM.

    `contam` / `pileup` are the stats dicts of the two upstream filters (already
    parsed); `total_reads` is the read count of the whole run, used to project
    how many usable footprints the full dataset holds.
    """
    # ---- thresholds (tuned on a 187-sample human cohort; see config.py)
    fp_len_lo = cfg["qc.footprint_len_lo"]
    fp_len_hi = cfg["qc.footprint_len_hi"]
    peak_frac_min = cfg["qc.read_len_peak_frac_min"]
    cds_enrich_min = cfg["qc.cds_enrich_min"]
    cds_strong_min = cfg["qc.cds_strong"]        # strong CDS-specific enrichment: ribo-seq
                                                 # even if the (jitter-smeared) periodicity
                                                 # looks weak
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
    # Strong periodicity may excuse an atypical footprint LENGTH (SRR29327151: 22 nt,
    # but 78% of genic reads in CDS and 65% in-frame), never an empty CDS. Periodicity
    # is measured on CDS reads alone, so a handful of genuine footprints buried in a
    # library that is otherwise intron / intergenic / ncRNA still scores "strongly
    # periodic": it says the ribosomes are real, not that the library is (SRR10513632:
    # 2% of genic reads in CDS -- 51% intron, 43% intergenic -- yet 62% in-frame;
    # SRR27534376: 18%, with more ncRNA (35%) than CDS (14%), yet 70% in-frame).
    cds_floor = cds_enrich_min

    idx, ncls, codons = _load_index(index_path)
    bamf = pysam.AlignmentFile(bam, "rb")

    len_hist: Counter = Counter()
    region_hist: Counter = Counter()
    # per-length frame counts within CDS (5'-end frame)
    frame_by_len: dict[int, list[int]] = defaultdict(lambda: [0, 0, 0])
    # the same, for the 13 mitochondrial CDS -- kept SEPARATE so it can never dilute the
    # nuclear periodicity the verdict is decided on, but measured, because in a
    # mitoribosome-profiling library this is the experiment
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
        # leakage-invariant 5' anchor (`pos - clip5`), as in qc.profile: a chance
        # match that extends the alignment moves pos and clip5 together
        anchor = (aln.reference_end + (aln.cigartuples[-1][1] if aln.cigartuples[-1][0] == 4 else 0)) \
            if aln.is_reverse else \
            (aln.reference_start - (aln.cigartuples[0][1] if aln.cigartuples[0][0] == 4 else 0))
        loc_hist[(chrom, aln.is_reverse, anchor)] += 1

        region, frame, on_cds = _classify(idx, ncls, chrom, p5, strand)
        region_hist[region] += 1
        if region == "CDS" and frame is not None:
            frame_by_len[L][frame] += 1
        elif region == "mito" and frame is not None:
            mito_frame_by_len[L][frame] += 1        # a separate pool; decides nothing
        # Which side of a coding gene are the reads on? A footprint is a piece of the
        # mRNA, so it is SENSE, and every ribo-seq library is ~100% sense here. A read
        # sitting on a CDS on the wrong strand was counted above as "intron" (a CDS is
        # inside a gene), so without this it is invisible -- see `_classify`.
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

    # ---- the same measurement on the 13 mitochondrial CDS. Reported, never used in the
    # verdict: it is a diagnostic that tells a mitoribosome-profiling library apart from
    # an ordinary one whose mito reads are degradation background.
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

    # ---- region fractions (gene-body = everything except intergenic, for CDS enrichment)
    reg_frac = {k: v / n for k, v in region_hist.items()}
    genic = n - region_hist.get("intergenic", 0) - region_hist.get("mito", 0)
    cds_of_genic = region_hist.get("CDS", 0) / genic if genic else 0.0

    top5p_locus_frac = loc_hist.most_common(1)[0][1] / n if loc_hist else 0.0

    # ---- which strand of the coding genes is this library on?
    n_on_cds = n_cds_sense + n_cds_anti
    cds_sense_frac = n_cds_sense / n_on_cds if n_on_cds else 0.0
    # An antisense-dominated deposit is reverse-complemented -- some submitters deposit
    # the RC of what they sequenced. It is not a strand curiosity: RiboMine scores the
    # SENSE strand, so its CDS enrichment collapses to ~0 and its periodicity is measured
    # on whatever minority sits the right way round, and the library is then refused for
    # "no CDS enrichment" with nothing said about why. Measured on SRR5750390: 77% of the
    # reads that touch a CDS are on the wrong strand.
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
    # 3-nt periodicity is the DEFINING hallmark of translating ribosomes.
    periodic_strong = enough_cds and periodicity_frac >= periodic_strong_min and tvd >= tvd_strong_min
    periodic_weak = enough_cds and periodicity_frac >= periodic_min and tvd >= tvd_min
    strong_cds = cds_of_genic >= cds_strong_min
    # a majority at one locus = contaminant, not footprints
    contaminant = top5p_locus_frac > single_locus_max
    # below this uniquely-mapped fraction the library maps poorly (contaminant /
    # non-genomic junk) -- a quality red flag even when the mappable minority is
    # periodic (cohort median 42%, p10 20%; broken libraries < 10%)
    low_unique = frac_unique is not None and frac_unique < min_unique_frac
    cds_floor_ok = cds_of_genic >= cds_floor
    # One 5' coordinate carrying as many reads as the WHOLE coding transcriptome is
    # a contaminant-dominated library, even when that fraction is far below the
    # outright single-locus majority above: the reads went to a single molecule
    # instead of to translation (SRR27534376: 12% at one locus, 12% in CDS). The
    # test is scale-free -- it compares the pile to the library's own coding signal
    # instead of picking an absolute cut-off (good libraries: top locus 0.5% median
    # against 40-70% of reads in CDS, a ratio of ~100x).
    cds_share = reg_frac.get("CDS", 0.0)
    locus_over_cds = (top5p_locus_frac >= locus_vs_cds * cds_share
                      and top5p_locus_frac >= locus_min_to_judge)

    # TI-seq: initiation drugs (harringtonine/LTM) freeze ribosomes at start codons
    # and elongation runs off, so the start-codon metagene peak towers over a
    # near-empty CDS body (start-codon enrichment ratio; GTI-seq SRR618773 ~120 vs
    # elongating ribo-seq, which tops out ~30).
    msc = meta["start_codon"]
    peak_off = max(range(-15, -8), key=lambda k: msc.get(k, 0)) if msc else -12
    peak_h = msc.get(peak_off, 0)
    body = [msc.get(peak_off + 3 * i, 0) for i in range(5, 16)
            if peak_off + 3 * i <= METAGENE_WIN[1]]
    body_mean = float(np.mean(body)) if body else 0.0
    start_ratio = peak_h / max(body_mean, 1.0)
    tiseq_like = peak_h >= tiseq_min_peak and start_ratio >= tiseq_ratio_min

    # ---- three-way verdict: RIBO-SEQ / TI-SEQ / NOT RIBO-SEQ or LOW QUALITY.
    # Translating-ribosome data (kept INCLUSIVE for variants) requires periodicity
    # or footprint-like+CDS-enriched reads, and no single-locus contamination; a
    # library that maps poorly (low unique mapping = non-genomic/multimapping junk)
    # is not usable even if the mappable minority is periodic.
    translating = (not contaminant) and cds_floor_ok and (not locus_over_cds) and (
        periodic_strong or (len_ok and cds_ok and (periodic_weak or strong_cds)))
    usable = translating and not low_unique
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
        f"CDS enrichment {'ok' if cds_ok else 'LOW'} ({cds_of_genic:.0%} of genic reads in CDS)",
    ]
    if enough_cds:
        lvl = "STRONG" if periodic_strong else ("weak" if periodic_weak else "ABSENT")
        reasons.append(f"periodicity {lvl} (in-frame {periodicity_frac:.0%}, TVD {tvd:.2f})")
    else:
        reasons.append(f"periodicity untested (only {cds_reads} CDS reads)")
    if reg_frac.get("mito", 0.0) >= MITO_DOMINANT:
        # Say this out loud. Everything above was measured on the NUCLEAR reads, which in
        # a mitoribosome library are the minority and arguably the contamination.
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

    # a single headline reason for a NOT-RIBO-SEQ-or-LOW-QUALITY call
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
        if antisense_deposit:
            msg = (f"the reads are ANTISENSE to the genes ({cds_sense_frac:.0%} of those on a "
                   f"CDS are on its strand), so the CDS is empty by construction")
            if len_ok:
                # The reads ARE footprint-length, so the strand is the whole story and
                # nothing else about this run can be read until it is settled: an empty CDS
                # is exactly what a reverse-complemented deposit looks like, and "no CDS
                # enrichment" would send the reader hunting for a biological answer to a
                # bookkeeping problem. Lead with it.
                why.insert(0, msg + " -- a reverse-complemented deposit would look exactly "
                                     "like this, and its footprints ARE the right length")
            else:
                # The reads are the wrong length for a footprint, so this is not ribo-seq
                # whichever way round it is -- in a random 100-run cohort all 9 antisense
                # deposits read 66-150 nt, i.e. reverse-stranded RNA-seq that the query's
                # recall net swept in. Say it, but do NOT lead with it: flipping the strand
                # would not turn a 150 nt read into a ribosome footprint, and a reader who
                # starts there has been sent to fix the wrong thing.
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
        },
        "start_codon_ratio": round(start_ratio, 1),
        "read_len_mode": mode_len,
        "read_len_peak_frac": round(peak_frac, 3),
        "periodicity_inframe_frac": round(periodicity_frac, 3),
        "periodicity_tvd_uniform": round(tvd, 3),
        "n_cds_reads": cds_reads,
        # mitochondrial translation, measured the same way but pooled separately. In an
        # ordinary cytosolic library these are a handful of degradation reads and mean
        # nothing; in a MITORIBOSOME-profiling library they are the entire experiment,
        # and the nuclear numbers above are then measuring the contamination.
        "mito_periodicity_inframe_frac": round(mito_periodicity, 3),
        "mito_periodicity_tvd_uniform": round(mito_tvd, 3),
        "n_mito_cds_reads": mito_cds_reads,
        "mito_dominated": bool(reg_frac.get("mito", 0.0) >= MITO_DOMINANT),
        "cds_frac_of_genic": round(cds_of_genic, 3),
        # strandedness: a real ribo-seq library is ~100% sense over CDS
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
        # Record what this verdict was decided against. The QC figure draws these as
        # its reference lines (the footprint-length window, the periodicity cut-offs),
        # so it audits the decision that was actually made rather than falling back on
        # its own module defaults -- a figure captioned with a verdict must not
        # contradict it. And a qc.json kept on disk stays interpretable after someone
        # re-tunes the config. Keys are named exactly as config.py names them.
        "thresholds": {
            "footprint_len_lo": fp_len_lo,
            "footprint_len_hi": fp_len_hi,
            "read_len_peak_frac_min": peak_frac_min,
            "cds_enrich_min": cds_enrich_min,
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
        # pile-up removal happens after mapping, so its fraction is of the mapped
        # reads; express it as a fraction of the sampled reads for one clean total
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

    # ---- projected usable reads: the reads that map (contaminant- and pile-up-
    #      filtered, = n_reads_scored) scaled from the sample to the whole dataset.
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
