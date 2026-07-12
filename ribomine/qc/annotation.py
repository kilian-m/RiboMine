"""Preprocess an Ensembl GTF into a compact cached index for ribo-seq QC.

The full GTF is ~1 GB; parsing it on every run is wasteful. It is parsed once
into per-chromosome integer arrays (pickled) that the QC stage loads in a
fraction of a second and turns into NCLS interval indexes.

What is stored, per chromosome:
  cds        CDS intervals with strand and reading frame (GTF col 8)
  utr5/utr3  5' / 3' UTR intervals
  exon_nc    exons of non-coding genes (gene_biotype != protein_coding)
  gene       gene spans (a position inside a gene but outside any exon = intron)
  start/stop codon A-positions, strand-aware (translation-direction first base)

Region priority at a position: MT chrom > CDS > 5'UTR > 3'UTR > ncRNA exon >
intron (inside a gene) > intergenic.
"""
from __future__ import annotations

import os
import pickle
import re
import tempfile
import time
from collections import defaultdict
from typing import Any

import numpy as np

from ..config import Config
from ..utils import LOG

FEATURES = {"CDS", "five_prime_utr", "three_prime_utr", "exon", "gene",
            "start_codon", "stop_codon"}
_BIOTYPE = re.compile(r'gene_biotype "([^"]+)"')


def _parse_gtf(gtf_path: str) -> dict[str, Any]:
    """Parse the GTF into the cache dict. Python lists first, arrays at the end."""
    cds = defaultdict(lambda: ([], [], [], []))       # start,end,strand,frame
    utr5 = defaultdict(lambda: ([], []))
    utr3 = defaultdict(lambda: ([], []))
    exon_nc = defaultdict(lambda: ([], []))
    gene = defaultdict(lambda: ([], []))
    startc = defaultdict(lambda: ([], []))            # pos,strand
    stopc = defaultdict(lambda: ([], []))

    t0 = time.time()
    n = 0
    with open(gtf_path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            feat = f[2]
            if feat not in FEATURES:
                continue
            n += 1
            if n % 2_000_000 == 0:
                LOG.debug("  ...%s feature lines (%.0fs)", f"{n:,}", time.time() - t0)
            chrom = f[0]
            s = int(f[3]) - 1        # GTF is 1-based inclusive -> 0-based half-open
            e = int(f[4])
            strand = 1 if f[6] == "+" else -1
            if feat == "CDS":
                a = cds[chrom]
                a[0].append(s); a[1].append(e); a[2].append(strand)
                a[3].append(int(f[7]) if f[7] in ("0", "1", "2") else 0)
            elif feat == "five_prime_utr":
                a = utr5[chrom]; a[0].append(s); a[1].append(e)
            elif feat == "three_prime_utr":
                a = utr3[chrom]; a[0].append(s); a[1].append(e)
            elif feat == "gene":
                a = gene[chrom]; a[0].append(s); a[1].append(e)
            elif feat == "exon":
                m = _BIOTYPE.search(f[8])
                if m and m.group(1) != "protein_coding":
                    a = exon_nc[chrom]; a[0].append(s); a[1].append(e)
            elif feat == "start_codon":
                pos = s if strand == 1 else e - 1     # first base in translation dir
                a = startc[chrom]; a[0].append(pos); a[1].append(strand)
            elif feat == "stop_codon":
                pos = s if strand == 1 else e - 1
                a = stopc[chrom]; a[0].append(pos); a[1].append(strand)

    def pack_iv(d, with_meta=False):
        out = {}
        for c, cols in d.items():
            rec = {"start": np.asarray(cols[0], np.int64),
                   "end": np.asarray(cols[1], np.int64)}
            if with_meta:
                rec["strand"] = np.asarray(cols[2], np.int8)
                rec["frame"] = np.asarray(cols[3], np.int8)
            out[c] = rec
        return out

    def pack_pts(d):
        return {c: {"pos": np.asarray(cols[0], np.int64),
                    "strand": np.asarray(cols[1], np.int8)}
                for c, cols in d.items()}

    cache = {
        "gtf": os.path.abspath(gtf_path),
        "gtf_mtime": os.path.getmtime(gtf_path),
        "cds": pack_iv(cds, with_meta=True),
        "utr5": pack_iv(utr5),
        "utr3": pack_iv(utr3),
        "exon_nc": pack_iv(exon_nc),
        "gene": pack_iv(gene),
        "start_codon": pack_pts(startc),
        "stop_codon": pack_pts(stopc),
    }
    LOG.info("parsed %s feature lines in %.0fs; %s CDS across %d chroms",
             f"{n:,}", time.time() - t0,
             f"{sum(len(v['start']) for v in cache['cds'].values()):,}",
             len(cache["cds"]))
    return cache


def build_index(gtf: str, out_pkl: str) -> str:
    """Parse `gtf` and write the pickled index to `out_pkl`. Returns `out_pkl`.

    Written to a temp file in the destination directory and `os.replace`d into
    position, so a half-written pickle is never visible to a concurrent worker
    (and a killed build leaves no file that `ensure_index` would trust).
    """
    if not os.path.isfile(gtf):
        raise ValueError(f"GTF not found: {gtf}")
    outdir = os.path.dirname(os.path.abspath(out_pkl)) or "."
    os.makedirs(outdir, exist_ok=True)

    cache = _parse_gtf(gtf)

    fd, tmp = tempfile.mkstemp(dir=outdir, suffix=".idx.tmp")
    try:
        with os.fdopen(fd, "wb") as fh:
            pickle.dump(cache, fh, protocol=4)
        os.replace(tmp, out_pkl)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    LOG.info("wrote %s (%.1f MB)", out_pkl, os.path.getsize(out_pkl) / 1e6)
    return out_pkl


def ensure_index(cfg: Config) -> str:
    """Return the annotation index path, building it if it is absent or stale.

    Safe to call from every worker of a process pool: the build goes to a
    private temp file and is `os.replace`d into place, so concurrent builders
    cannot see (or write) a partial pickle -- the worst case is that the same
    index is parsed twice and the second one atomically overwrites an identical
    first.
    """
    gtf = cfg.ref("gtf")
    if not gtf or not os.path.isfile(gtf):
        raise ValueError(f"reference.gtf does not exist: {gtf!r}")
    out_pkl = cfg.annotation_index

    # a pickle older than the GTF was built from a different annotation
    if os.path.isfile(out_pkl) and os.path.getmtime(out_pkl) >= os.path.getmtime(gtf):
        LOG.debug("annotation index up to date: %s", out_pkl)
        return out_pkl

    LOG.info("building annotation index from %s", gtf)
    return build_index(gtf, out_pkl)
