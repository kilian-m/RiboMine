"""Remove contaminant read pile-ups by position -- no adapter sequences needed.

Adapter/primer dimers and fixed contaminants collapse onto a *single* genomic
5' position as one identical, fixed-length molecule (e.g. the TruSeq adapter maps
to chr15:56,885,929 in many libraries, 40-70 % of reads, all 22-27 nt, all
identical). A genuinely translated codon can also be highly covered, but its
ribosome footprints always **spread across lengths** (~26-34 nt) -- so length
concentration at a position is what separates a fixed contaminant from real
translation.

A 5' position is dropped when it is BOTH:
  * over-represented -- >= pileup_min_count reads and >= pileup_min_frac of all
    reads, and
  * a single fixed length -- the modal footprint length is >= pileup_len_conc of
    the reads there (a translated position never is; its footprints vary in
    length).

The fixed-length test misses a pile whose *insert* varies in length: the reads
then align at 21/22/23 nt from the same 5' base and look length-diverse, even
though they are one molecule (SRR2096968: 76 % of the library at chr15:56,885,929
at length concentration 0.67; the Y:3,367,775 tRNA fragment aligns at both 25 and
29 nt, concentration 0.60). So a position is also dropped on sheer dominance:

  * >= pileup_max_frac of ALL surviving reads at one 5' coordinate, whatever its
    length spread.

10 % is far above anything real: across the 187-sample cohort the top 5' position
of a library that passes ribo-seq QC holds a median 0.5 % of reads and 3.8 % at
the 90th percentile, while every pile above ~10 % is demonstrably junk -- adapter
dimers, poly(T)+P7 primer dimers, adapter-homologous force-alignments (100 % of
their reads carry a mismatch) and tRNA/rRNA fragments that slipped past the
bowtie contaminant filter. No translated codon holds a tenth of a library.

The dominance test is applied to a fixpoint: fractions are recomputed over the
surviving reads and the sweep repeats, so a second pile cannot hide behind a
larger one (SRR10491340 has three).

Widening pileup_len_conc to a +-1 nt window would catch these piles too, but it
keys on the wrong property (length spread, not abundance) and clips genuine
footprint peaks.

This is fully data-driven and adapter-agnostic.
"""
from __future__ import annotations

import os
from collections import Counter, defaultdict

import pysam

from ..config import Config
from ..utils import LOG


def filter_bam(in_bam: str, out_bam: str, cfg: Config, *, label: str = "") -> dict:
    """Drop every read whose 5' end sits on a contaminant pile-up position.

    Unmapped / secondary / supplementary records are passed through untouched:
    the pile-up rule is about primary alignment positions, and the QC stage still
    wants the mapping statistics intact. Returns the dict for `Sample.pileup_json`.
    """
    label = label or os.path.basename(in_bam)
    if not os.path.isfile(in_bam):
        raise ValueError(f"input BAM not found: {in_bam}")

    min_count = int(cfg.get("contaminants.pileup_min_count", 30))
    min_frac = float(cfg.get("contaminants.pileup_min_frac", 0.005))
    len_conc = float(cfg.get("contaminants.pileup_len_conc", 0.85))
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
    """The read's 5' genomic base: the pile-up collapses onto exactly this."""
    return (a.reference_name, a.is_reverse,
            a.reference_end - 1 if a.is_reverse else a.reference_start)


def _find_pileups(bam_path: str, *, min_count: int, min_frac: float,
                  len_conc: float, max_frac: float) -> tuple[dict, int]:
    """The pile-up positions, by the fixed-length rule and then by dominance.

    Two passes over the BAM, deliberately. Both rules can only ever fire on a
    position holding >= max(min_count, min_frac * n) reads -- they test that
    first -- so the *length histogram*, which is the expensive part, is only ever
    needed for the handful of positions that clear the count. Keeping a length
    Counter for every distinct 5' position instead costs a dict object per
    position: harmless on the 200k-read QC sample, but this also runs on the
    FULL-run BAM (stage 4, tens of millions of distinct 5' positions) times
    `project.jobs` concurrent samples -- tens of GB of live dicts, i.e. an OOM.

    So: pass 1 counts reads per position (one int each, and it hands us `n`),
    the count dict is then narrowed to the candidates -- at most 1/min_frac of
    them, 200 by default -- and pass 2 re-reads the BAM to build the length
    histograms for just those. Re-reading a file that is sitting on disk is
    cheap next to running out of memory.

    The algorithm is untouched: same counts, same candidate set, same length
    histograms where they are consulted, and the same first-seen iteration order
    that the dominance sweep breaks its ties on -- so the calls are identical.
    """
    # A plain dict of int, not a Counter of Counters. Insertion order (= first
    # occurrence in the BAM) is load-bearing: the dominance sweep below takes the
    # FIRST position that qualifies, not the largest.
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
            # one shared str per chromosome: pysam builds `reference_name` fresh on
            # every access, and these keys are held for the whole pass
            k = (names.setdefault(chrom, chrom), rev, p5)
            count[k] = count.get(k, 0) + 1
    finally:
        bam.close()

    thr = max(min_count, min_frac * n)
    # Only these can ever be a pile-up; everything else is dropped here, which is
    # what makes pass 2 (and the rest of this function) cost nothing. The dict
    # comprehension preserves the first-seen order of the survivors.
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

    # Dominance, to a fixpoint: a position holding >= max_frac of the *surviving*
    # reads is one molecule, not translation. Removing it raises every other
    # position's share, so re-measure and sweep again -- otherwise a second pile
    # stays hidden behind a larger one.
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
    """Pass 2: footprint-length histogram for the candidate positions ONLY.

    Reading the BAM a second time to fill a few hundred Counters, rather than
    carrying tens of millions of them through pass 1 (see `_find_pileups`).
    Lengths are accumulated in BAM order, so each histogram is identical -- down
    to its own insertion order, which `max(ld, key=ld.get)` breaks modal-length
    ties on -- to the one the single-pass version built.
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
