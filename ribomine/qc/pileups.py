"""Remove contaminant read pile-ups by position; no adapter sequences needed.

Adapter and primer dimers, small RNAs (miRNA) and other fixed contaminants stack
on a single genomic 5' position as one molecule of (nearly) one length. A highly
translated codon can be covered as deeply, but its footprints spread over lengths
(~26-34 nt). A 5' position is removed by either of two rules (`contaminants.*`):

  fixed-length  >= pileup_min_count reads and >= pileup_min_frac of all reads at
                the position, and the modal aligned length holds >= pileup_len_conc
                of them. The default of 0.75 also catches mature miRNAs, whose
                length concentration is 0.75-0.85; no CDS or start-codon position
                was found in that band.
  dominant      >= pileup_max_frac of the surviving reads at the position, whatever
                its length spread (a pile whose insert length varies). Applied to a
                fixpoint: fractions are recomputed after each removal, so a second
                pile cannot hide behind a larger one. The default of 10 % is well
                above real signal: in libraries passing QC the top 5' position
                holds a median 0.5 % of reads (3.8 % at the 90th percentile).
"""
from __future__ import annotations

import os
from collections import Counter, defaultdict

import pysam

from ..config import Config
from ..utils import LOG


def filter_bam(in_bam: str, out_bam: str, cfg: Config, *, label: str = "") -> dict:
    """Drop every primary alignment whose 5' end sits on a pile-up position.

    Unmapped, secondary and supplementary records are written unchanged.
    Returns the stats dict for `Sample.pileup_json`.
    """
    label = label or os.path.basename(in_bam)
    if not os.path.isfile(in_bam):
        raise ValueError(f"input BAM not found: {in_bam}")

    min_count = int(cfg.get("contaminants.pileup_min_count", 30))
    min_frac = float(cfg.get("contaminants.pileup_min_frac", 0.005))
    len_conc = float(cfg.get("contaminants.pileup_len_conc", 0.75))
    max_frac = float(cfg.get("contaminants.pileup_max_frac", 0.10))

    bad, n = _find_pileups(in_bam, min_count=min_count, min_frac=min_frac,
                           len_conc=len_conc, max_frac=max_frac)

    os.makedirs(os.path.dirname(os.path.abspath(out_bam)) or ".", exist_ok=True)
    bam = pysam.AlignmentFile(in_bam, "rb")
    out = pysam.AlignmentFile(out_bam, "wb", template=bam)
    n_out = n_removed = 0
    try:
        for a in bam:
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                out.write(a)
                continue
            if _five_prime(a) in bad:
                n_removed += 1
            else:
                out.write(a)
                n_out += 1
    finally:
        out.close()
        bam.close()

    top = sorted(bad.items(), key=lambda kv: -kv[1][0])[:8]
    stats = {
        "label": label,
        "n_reads_in": n,
        "n_reads_removed": n_removed,
        "n_reads_kept": n_out,
        "frac_removed": round(n_removed / max(n, 1), 4),
        "n_pileup_positions": len(bad),
        "n_dominant_positions": sum(1 for v in bad.values() if v[3] == "dominant"),
        "top_pileups": [{"chrom": k[0], "strand": "-" if k[1] else "+", "pos": k[2],
                         "reads": v[0], "frac": round(v[0] / max(n, 1), 4),
                         "footprint_len": v[1], "len_conc": round(v[2], 3),
                         "rule": v[3]}
                        for k, v in top],
    }
    msg = (f"{label}: removed {n_removed:,}/{n:,} reads "
           f"({100 * stats['frac_removed']:.0f}%) at {len(bad)} pile-up position(s)")
    if top:
        msg += (f"; top {top[0][0][0]}:{top[0][0][2]} "
                f"{100 * top[0][1][0] / max(n, 1):.0f}% ({top[0][1][1]}nt)")
    LOG.info("%s", msg)
    return stats


def _five_prime(a) -> tuple[str, bool, int]:
    """(chrom, is_reverse, position) of the read's 5' aligned base."""
    return (a.reference_name, a.is_reverse,
            a.reference_end - 1 if a.is_reverse else a.reference_start)


def _find_pileups(bam_path: str, *, min_count: int, min_frac: float,
                  len_conc: float, max_frac: float) -> tuple[dict, int]:
    """Find the pile-up positions. Returns ({position: (reads, modal length,
    length concentration, rule)}, number of primary alignments).

    Two passes over the BAM to bound memory on full-run BAMs: pass 1 counts reads
    per 5' position, pass 2 builds length histograms only for the candidates,
    i.e. positions with >= max(min_count, min_frac * n) reads.
    """
    # Insertion order (first occurrence in the BAM) matters: the dominance sweep
    # takes the first qualifying position, not the largest.
    count: dict[tuple[str, bool, int], int] = {}
    names: dict[str, str] = {}
    n = 0
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        for a in bam:
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                continue
            n += 1
            chrom, rev, p5 = _five_prime(a)
            # share one str per chromosome across keys (pysam returns a new one each time)
            k = (names.setdefault(chrom, chrom), rev, p5)
            count[k] = count.get(k, 0) + 1
    finally:
        bam.close()

    thr = max(min_count, min_frac * n)
    # keep only the candidates (first-seen order is preserved)
    count = {k: c for k, c in count.items() if c >= thr}
    lens = _length_hists(bam_path, set(count))
    if any(sum(lens[k].values()) != c for k, c in count.items()):
        raise RuntimeError(f"{bam_path} changed between the two pile-up passes "
                           f"(position counts disagree)")

    bad: dict = {}
    for k, c in count.items():
        ld = lens[k]
        lmode = max(ld.values()) / c            # length concentration at this position
        if lmode >= len_conc:
            bad[k] = (c, max(ld, key=ld.get), lmode, "fixed-length")

    # Dominance, to a fixpoint: remove a position holding >= max_frac of the
    # surviving reads, recompute the shares, and repeat.
    while True:
        alive = n - sum(v[0] for v in bad.values())
        if alive <= 0:
            break
        hit = None
        for k, c in count.items():
            if k in bad:
                continue
            if c / alive >= max_frac:
                hit = k
                break
        if hit is None:
            break
        ld = lens[hit]
        c = count[hit]
        bad[hit] = (c, max(ld, key=ld.get), max(ld.values()) / c, "dominant")
    return bad, n


def _length_hists(bam_path: str, positions: set) -> dict:
    """Pass 2: aligned-length histogram for each candidate position.

    Lengths are added in BAM order; modal-length ties go to the length seen first.
    """
    lens: dict = defaultdict(Counter)
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        for a in bam:
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                continue
            k = _five_prime(a)
            if k in positions:
                lens[k][a.query_alignment_length] += 1
    finally:
        bam.close()
    return lens
