"""Turn locally-aligned reads into the positional statistics that reveal the
read architecture.

The one idea this module rests on
--------------------------------
For a read aligned in local mode, STAR's soft-clip boundary is *not* a reliable
marker of where the footprint starts. If the last base of a 5' UMI happens to
match the genome base next to the footprint (probability 1/4), STAR extends the
alignment by one and the clip shrinks. So clip lengths leak.

What does *not* leak is the genomic coordinate that the alignment implies for
read position 0. Extending the alignment leftwards by one decrements both
`pos` and `clip5`, so `pos - clip5` is invariant. Anchoring on that, every read
position p maps to a fixed genomic base, and we can ask the only question that
matters:

    does read base p match the genome base the alignment implies for it?

Footprint bases match ~99% of the time. UMI, adapter and barcode bases are not
genomic, so they match ~25% of the time by chance. The transition between those
two regimes is the architecture.

Anchors
-------
5' side: read position 0. The 5' construct has the same length in every read,
so this anchor is exact.

3' side: read position 0 is useless because the footprint length varies. The
alignment end is no good either -- STAR extends until it hits a mismatch, so
the base just past it mismatches *by construction* (match rate 0.00). The two
honest anchors are the **adapter start** and the **read's own 3' end**, and
this module profiles around both.

Reads carrying indels or splice junctions break the read-pos -> genome-index
map and are skipped (a couple of percent).

Everything here is descriptive. Thresholds and calling live in
`ribomine.arch.infer` so they can be tuned without re-reading BAMs -- which is
why `profile_bam` takes no config.
"""
from __future__ import annotations

import logging
import os
from collections import Counter
from dataclasses import dataclass

import numpy as np
import pysam

LOG = logging.getLogger("ribomine.qc.profile")

# --- internals of the profiler: deliberately NOT config keys ----------------
# These describe the shape of the measurement (how many positions to profile,
# how long a de-novo seed is), not a tunable calling threshold.
MAXP5 = 24            # read positions profiled from the 5' end
MAXW = 30             # window profiled around the 3' anchors
SEED_K = 12           # k-mer length used to discover an adapter de novo
MIN_CORE_ENT = 1.2    # bits; aligned cores below this are homopolymer junk
MAX_CORE_BASE = 0.70  # ... as is any core dominated by a single base
BASES = "ACGT"
B2I = {b: i for i, b in enumerate(BASES)}

# 3' adapters used in ribo-seq / small-RNA libraries.
ADAPTER_PANEL = [
    ("illumina_truseq", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"),
    ("illumina_smallrna_ra3", "TGGAATTCTCGGGTGCCAAGG"),
    ("ingolia_linker", "CTGTAGGCACCATCAAT"),
    ("nextera", "CTGTCTCTTATACACATCT"),
    ("polyA", "AAAAAAAAAAAAAAAAAAAA"),
    ("polyG", "GGGGGGGGGGGGGGGGGGGG"),
]
_PANEL_SEQ = dict(ADAPTER_PANEL)

# An adapter must be visible in this share of the reads before we will anchor the 3'
# profile on it. The architecture caller applies the SAME gate to the same quantity
# (`frac_anchored`, via architecture.min_anchor_frac), so the two must move together:
# if this one stayed frozen, lowering the config key would be inert -- the anchor would
# never be selected here in the first place, and the caller would silently fall back to
# the weaker read-3'-end anchor. Hence the caller passes its value in, and this is only
# the default.
MIN_PANEL_FRAC = 0.15

# ... unless the adapter is hiding behind a long insert. When the molecule is about
# as long as the read, only the short-footprint minority sequences far enough to
# reach the adapter -- but those reads read the construct out perfectly well. An
# adapter seen in this many reads, sitting at a FIXED distance from the footprint
# end (`MIN_GAP_CONC`), is a real adapter no matter how small its share.
MIN_ANCHOR_READS = 300
MIN_GAP_CONC = 0.30

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# base -> 0..3 for the composition counts; anything else (N) -> 4, counted in the
# denominator but in no base's numerator, exactly as the reference's B2I.get(b)
_CODE = np.full(256, 4, dtype=np.uint8)
_CODE[np.frombuffer(BASES.encode("ascii"), dtype=np.uint8)] = np.arange(4, dtype=np.uint8)


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def max_mm(overlap: int) -> int:
    """cutadapt-style error budget: 0 errors below 9 nt, then ~12%."""
    return int(0.12 * overlap)


def find_adapter(read: str, adapter: str, min_start: int, min_overlap: int = 7) -> int:
    """Leftmost i >= min_start where read[i:] is a prefix of `adapter`.

    Handles the common case where the adapter runs off the end of the read, so
    only its first few bases are visible.
    """
    n, m = len(read), len(adapter)
    for i in range(max(0, min_start), n - min_overlap + 1):
        ov = min(m, n - i)
        if ov < min_overlap:
            break
        budget = max_mm(ov)
        mm = 0
        ok = True
        for a, b in zip(read[i : i + ov], adapter[:ov]):
            if a != b:
                mm += 1
                if mm > budget:
                    ok = False
                    break
        if ok:
            return i
    return -1


def cigar_is_simple(ct) -> bool:
    """True when the read has no I/D/N, so read pos -> genome pos is a shift."""
    return all(op in (0, 4, 7, 8) for op, _ in ct)


def low_complexity(s: str) -> bool:
    """Adapter dimers and poly-A tails map to genomic homopolymers. Their
    'footprint' is meaningless, and they drag every profile toward noise."""
    if not s:
        return True
    counts = [s.count(b) for b in BASES]
    n = len(s)
    if max(counts) / n > MAX_CORE_BASE:
        return True
    p = np.array([c / n for c in counts if c])
    return float(-(p * np.log2(p)).sum()) < MIN_CORE_ENT


@dataclass(slots=True)
class _Read:
    """One usable alignment, reduced to what the profiles need."""

    seq: str            # read as sequenced (5'->3' of the original read)
    n: int              # len(seq)
    codes: np.ndarray   # uint8, 0..3 = ACGT, 4 = other
    matches: np.ndarray # bool, read base p == implied genome base
    valid: np.ndarray   # bool, position p has an implied genome base at all
    aln_s: int          # first aligned read position (5' clip length)
    aln_e: int          # first read position past STAR's aligned block
    fp_end: int         # ... extended right through chance matches: footprint end


def load_reads(bam_path: str, fasta_path: str, max_reads: int):
    """Read the BAM once; keep the read, the match mask and the two boundaries.

    Also count, per used read, the leakage-invariant genomic coordinate of the
    5' end (`pos - clip5`, which does not move when a chance match extends the
    alignment). Genuine footprints spread across thousands of loci, so no single
    5' coordinate holds more than ~1-2% of reads; a large share at one coordinate
    means the reads are a single over-represented species (an adapter dimer that
    aligns to an adapter-like locus, a spike-in, a super-abundant contaminant) --
    not footprints.
    """
    fa = pysam.FastaFile(fasta_path)
    bam = pysam.AlignmentFile(bam_path, "rb")
    recs: list[_Read] = []
    loc_hist: Counter = Counter()
    n_seen = n_skip_cigar = n_skip_ctx = n_skip_lowcomp = 0

    for aln in bam:
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        n_seen += 1
        if len(recs) >= max_reads:
            break
        ct = aln.cigartuples
        if not ct or not cigar_is_simple(ct):
            n_skip_cigar += 1
            continue

        seq = aln.query_sequence.upper()
        L = len(seq)
        clipL = ct[0][1] if ct[0][0] == 4 else 0
        clipR = ct[-1][1] if ct[-1][0] == 4 else 0

        # every read position must have an implied genome base, so the context
        # has to reach `clip` bases beyond the aligned block on either side
        pad = min(250, L + 10)
        lo, hi = aln.reference_start - pad, aln.reference_end + pad
        if lo < 0 or hi > fa.get_reference_length(aln.reference_name):
            n_skip_ctx += 1
            continue
        ctx = fa.fetch(aln.reference_name, lo, hi).upper()

        if aln.is_reverse:
            read, ctxr = revcomp(seq), revcomp(ctx)
            clip5, clip3 = clipR, clipL
        else:
            read, ctxr = seq, ctx
            clip5, clip3 = clipL, clipR

        if low_complexity(read[clip5 : L - clip3]):
            n_skip_lowcomp += 1
            continue

        # off = pad - clip5 is leakage-invariant: read pos p -> ctxr[off + p].
        # Compared byte-wise (both strings are upper-case), so an N in the read
        # matches an N in the genome, exactly as the character comparison did.
        off = pad - clip5
        rb = np.frombuffer(read.encode("ascii"), dtype=np.uint8)
        cb = np.frombuffer(ctxr.encode("ascii"), dtype=np.uint8)
        i = off + np.arange(L, dtype=np.int64)
        valid = (i >= 0) & (i < cb.size)
        matches = np.zeros(L, dtype=bool)
        matches[valid] = rb[valid] == cb[i[valid]]

        aln_s, aln_e = clip5, L - clip3
        # extend right through chance matches: a definition of the footprint end
        # that depends on the genome, not on STAR's scoring
        e = aln_e
        while e < L and valid[e] and matches[e]:
            e += 1

        recs.append(_Read(seq=read, n=L, codes=_CODE[rb], matches=matches, valid=valid,
                          aln_s=aln_s, aln_e=aln_e, fp_end=e))
        # leakage-invariant 5' genomic anchor (implied coordinate of read pos 0)
        anchor5 = (aln.reference_end + clip5) if aln.is_reverse \
            else (aln.reference_start - clip5)
        loc_hist[(aln.reference_name, aln.is_reverse, anchor5)] += 1

    return recs, loc_hist, n_seen, n_skip_cigar, n_skip_ctx, n_skip_lowcomp


def profile_window(recs: list[_Read], anchors: list[int | None], lo: int, hi: int):
    """Composition + genome-match rate at offsets [lo, hi) from a per-read anchor.

    Composition and match rate have different denominators: a read position
    always has a base, but it only has an implied genome base when the fetched
    context reached that far. Counting context-less positions as mismatches
    would push the match rate below the 1/4 chance floor.

    (Vectorised per read over the window; the counts are integer-valued, so the
    result is identical to the reference's position-by-position accumulation.)
    """
    w = hi - lo
    comp = np.zeros((w, 4))
    match = np.zeros(w)
    n_comp = np.zeros(w)
    n_match = np.zeros(w)
    for rec, a in zip(recs, anchors):
        if a is None:
            continue
        base = a + lo                      # read position of window index 0
        j0 = max(0, -base)                 # ... clipped to the read
        j1 = min(w, rec.n - base)
        if j1 <= j0:
            continue
        sl = slice(base + j0, base + j1)
        n_comp[j0:j1] += 1
        codes = rec.codes[sl]
        for b in range(4):
            comp[j0:j1, b] += codes == b   # an N contributes to no base
        n_match[j0:j1] += rec.valid[sl]
        match[j0:j1] += rec.matches[sl]    # matches is False wherever valid is False
    return comp, match, n_comp, n_match


def discover_seed(recs: list[_Read], n_min: int):
    """Most common k-mer among the non-genomic 3' tails."""
    cnt: Counter = Counter()
    n_room = 0
    for rec in recs:
        if rec.n - rec.fp_end < SEED_K:
            continue
        n_room += 1
        read = rec.seq
        seen = set()
        for i in range(rec.fp_end, rec.n - SEED_K + 1):
            km = read[i : i + SEED_K]
            if km not in seen and "N" not in km:
                seen.add(km)
                cnt[km] += 1
    if not cnt:
        return None, 0, n_room, []
    top = cnt.most_common(8)
    kmer, c = top[0]
    if c < n_min:
        return None, c, n_room, top
    return kmer, c, n_room, top


def rate(a: np.ndarray, b: np.ndarray) -> list:
    return np.where(b > 0, a / np.maximum(b, 1), np.nan).tolist()


def comp_frac(c: np.ndarray) -> list:
    s = c.sum(axis=1, keepdims=True)
    return np.where(s > 0, c / np.maximum(s, 1), np.nan).tolist()


def profile_bam(bam: str, fasta: str, *, label: str, max_reads: int = 150_000,
                min_panel_frac: float = MIN_PANEL_FRAC) -> dict:
    """Positional match/composition profile of a locally-aligned BAM.

    Purely descriptive: every number here is a measurement, and the architecture
    call (`ribomine.arch.infer`) is made from this dict alone, never from the BAM.

    `min_panel_frac` is the share of reads an adapter must appear in before it is
    used as the 3' anchor; the pipeline passes `architecture.min_anchor_frac` so
    that this gate and the caller's cannot drift apart (see MIN_PANEL_FRAC).
    """
    recs, loc_hist, n_seen, n_skip_cigar, n_skip_ctx, n_skip_lc = load_reads(
        bam, fasta, max_reads
    )
    n = len(recs)
    if n == 0:
        raise ValueError(f"no usable alignments in {bam}")

    top5p_count = loc_hist.most_common(1)[0][1] if loc_hist else 0
    top5p_locus_frac = top5p_count / n
    n_distinct_5p_loci = len(loc_hist)

    len_hist = Counter(r.n for r in recs)
    clip5_hist = Counter(r.aln_s for r in recs)
    clip3_hist = Counter(r.n - r.aln_e for r in recs)
    fp_len_hist = Counter(r.fp_end - r.aln_s for r in recs)
    tail_hist = Counter(r.n - r.fp_end for r in recs)   # non-genomic 3' length

    # ---- 5' side: anchor on read position 0 (exact, no leakage)
    c5, m5, n5c, n5 = profile_window(recs, [0] * n, 0, MAXP5)

    # ---- 3' side, anchor A: the read's own 3' end (offset 0 = last base)
    ct, mt, ntc, nt_ = profile_window(recs, [r.n - 1 for r in recs], -(MAXW - 1), 1)
    ct, mt, nt_ = ct[::-1], mt[::-1], nt_[::-1]   # index 0 = last base of read

    # ---- 3' side, anchor B: the adapter start
    panel_hits: dict[str, int] = {}
    panel_starts: dict[str, list[int]] = {}
    for name, seq in ADAPTER_PANEL:
        starts = [find_adapter(r.seq, seq, min_start=r.fp_end - 2) for r in recs]
        panel_hits[name] = sum(1 for s in starts if s >= 0)
        panel_starts[name] = starts

    seed_kmer, seed_count, n_room, seed_top = discover_seed(
        recs, n_min=max(50, int(0.10 * n))
    )

    # The 3' adapter is the one ligated *nearest the insert*, not the most
    # frequent. Ingolia-style libraries sequence straight through the ligated
    # linker into the TruSeq read primer, so the reads carry both; the linker at
    # the smaller read position is the real adapter, the outer match is
    # sequencing read-through. Among adapters seen in >= 15% of reads (excluding
    # the homopolymer artefacts), pick the one with the smallest median start.
    def median_start(name: str) -> int:
        s = sorted(x for x in panel_starts[name] if x >= 0)
        return s[len(s) // 2] if s else 10 ** 9

    def gap_concentration(name: str) -> float:
        """How sharply the adapter sits at a fixed distance from the footprint end.

        A real ligated adapter is separated from the footprint by a construct of
        FIXED length (a UMI, a barcode, or nothing), so the gap has a sharp mode. A
        chance 7-mer match to a panel sequence lands anywhere, so its gaps scatter.
        This is what lets a minority-of-reads adapter be trusted (below).
        """
        gaps = Counter(s - r.fp_end for r, s in zip(recs, panel_starts[name])
                       if s >= 0 and s >= r.fp_end - 2)
        tot = sum(gaps.values())
        return max(gaps.values()) / tot if tot else 0.0

    def qualifies(name: str) -> bool:
        if name in ("polyA", "polyG"):
            return False
        hits = panel_hits[name]
        if hits >= min_panel_frac * n:
            return True
        # A MINORITY of reads showing the adapter does not mean there is no adapter.
        # When the insert is about as long as the read -- 2 nt UMI + ~32 nt footprint
        # + 10 nt construct in a 46 nt read -- the adapter simply falls off the end of
        # most reads, and only the short-footprint minority reaches it. Those reads
        # still read the construct out exactly (the architecture is a per-molecule
        # property, not a per-read one), and they are what an adapter-anchored profile
        # needs. Rejecting them forces the caller onto the read's own 3' end, which
        # SMEARS the construct across positions and fabricates a long "UMI".
        # So accept a minority adapter on two conditions: enough reads to measure a
        # construct from, and a gap to the footprint end that is actually fixed.
        return hits >= MIN_ANCHOR_READS and gap_concentration(name) >= MIN_GAP_CONC

    qual = [name for name in panel_hits if qualifies(name)]
    if qual:
        best_panel = min(qual, key=median_start)
        use_panel = True
    else:
        # nothing but a homopolymer: accept it so infer() can refuse cleanly
        best_panel = max(panel_hits, key=lambda k: panel_hits[k]) if panel_hits else None
        use_panel = best_panel is not None and panel_hits[best_panel] >= min_panel_frac * n

    anchor_kind, anchor_seq = "none", ""
    anchors: list[int | None] = [None] * n
    if use_panel:
        anchor_kind, anchor_seq = f"panel:{best_panel}", _PANEL_SEQ[best_panel]
        anchors = [s if s >= 0 else None for s in panel_starts[best_panel]]
    elif seed_kmer:
        anchor_kind, anchor_seq = "denovo", seed_kmer
        anchors = [find_adapter(r.seq, seed_kmer, min_start=r.fp_end - 2) for r in recs]
        anchors = [s if s >= 0 else None for s in anchors]

    n_anchored = sum(a is not None for a in anchors)
    if n_anchored:
        ca, ma, nac, na = profile_window(recs, anchors, -MAXW, MAXW)
    else:
        ca = np.full((2 * MAXW, 4), np.nan)
        ma = np.full(2 * MAXW, np.nan)
        nac = na = np.zeros(2 * MAXW)

    # Distance from the footprint end to the adapter start. A fixed-length
    # construct (UMI, barcode) gives a sharp mode; a poly(A) tail or any other
    # variable-length insert gives a broad distribution.
    gap_hist = Counter(
        a - r.fp_end for r, a in zip(recs, anchors)
        if a is not None and a >= r.fp_end - 2
    )

    out = {
        "label": label or os.path.basename(bam),
        "n_alignments_seen": n_seen,
        "n_used": n,
        "n_skipped_indel_or_junction": n_skip_cigar,
        "n_skipped_contig_edge": n_skip_ctx,
        "n_skipped_low_complexity": n_skip_lc,
        "frac_low_complexity": n_skip_lc / max(n + n_skip_lc, 1),
        # concentration of reads on a single 5' genomic coordinate: high => a
        # single over-represented species (adapter dimer / contaminant / spike-in)
        "top5p_locus_frac": top5p_locus_frac,
        "n_distinct_5p_loci": n_distinct_5p_loci,
        "read_len_hist": dict(sorted(len_hist.items())),
        "clip5_hist": dict(sorted(clip5_hist.items())),
        "clip3_hist": dict(sorted(clip3_hist.items())),
        "footprint_len_hist": dict(sorted(fp_len_hist.items())),
        "tail_len_hist": dict(sorted(tail_hist.items())),
        "fpend_to_adapter_gap_hist": dict(sorted(gap_hist.items())),
        # 5'-anchored, index 0 = first base of the read
        "p5_comp": comp_frac(c5),
        "p5_match": rate(m5, n5),
        "p5_n": n5.tolist(),
        # read-3'-end-anchored, index 0 = last base of the read, growing inwards
        "t3_comp": comp_frac(ct),
        "t3_match": rate(mt, nt_),
        "t3_n": nt_.tolist(),
        # adapter-anchored, index MAXW = first base of the adapter
        "anchor_kind": anchor_kind,
        "anchor_seq": anchor_seq,
        "anchor_offset": MAXW,
        "n_anchored": n_anchored,
        "frac_anchored": n_anchored / n,
        "adap_comp": comp_frac(ca),
        "adap_match": rate(ma, na),
        "adap_n": na.tolist(),
        "panel_hits": panel_hits,
        "panel_frac": {k: v / n for k, v in panel_hits.items()},
        "seed_kmer": seed_kmer or "",
        "seed_count": seed_count,
        "seed_room": n_room,
        "seed_top": [[k, c] for k, c in seed_top],
        "bases": BASES,
    }
    LOG.info(
        "%s: %d reads | anchor=%s (%.0f%% of reads) | skipped %d indel/junction, "
        "%d low-complexity",
        out["label"], n, anchor_kind, 100 * out["frac_anchored"], n_skip_cigar, n_skip_lc,
    )
    return out
