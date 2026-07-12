"""UMI deduplication of the mapped reads (`umi_tools dedup`).

`arch.trim` moved the random-templated content of the read -- the 5' and 3' UMI
blocks the architecture found -- into the read name as `_<UMI>`, umi_tools style.
This stage collapses reads that share a mapping position *and* a UMI: those are
PCR copies of one molecule.

It is OFF by default (`process.umi_dedup`), and that default is deliberate. The
CALLER decides whether to deduplicate at all; this module's job is to refuse when
deduplication would be wrong. Two things have to be true before it is meaningful:

1. **The reads must actually carry a UMI.** If the architecture found no
   random-templated content (`dedup_umi_len == 0`), `arch.trim` wrote no `_<UMI>`
   tag, and there is nothing to distinguish two independent ribosome footprints
   that start at the same codon from a PCR duplicate pair. Running umi_tools
   anyway would silently collapse genuine footprints into one read -- at ribo-seq
   depth, where a well-expressed CDS has hundreds of reads on the same start,
   that quietly destroys the signal it is supposed to clean. So: raise.
2. **The INPUT BAM must be coordinate-sorted and indexed** -- umi_tools walks it
   by position. We check and say so if it is not.

What this module does NOT do is index its output: `out_bam` is a temporary the
caller renames over the final BAM, so a `.bai` built here would be orphaned by
that rename (and would then be a *stale* index sitting next to the final BAM).
Indexing the final artefact belongs to whoever owns it -- `pipeline.process_sample`
sorts and indexes after the rename.
"""
from __future__ import annotations

import os

import pysam

from ..config import Config
from ..utils import LOG, nonempty, require_tools, run

# How many reads to look at when asking "do these reads carry a UMI?" -- the tag
# is written by our own trimmer, so it is either there for (essentially) every
# read or for none; a head sample settles it.
UMI_HEAD_READS = 10_000
UMI_TAG_MIN_FRAC = 0.5
UMI_BASES = set("ACGTN")
MAX_UMI_TAG = 24        # longer than any real UMI: a read name that merely happens
                        # to end in "_<something long>" is not a UMI tag
UMI_SEPARATOR = "_"


def dedup(bam: str, out_bam: str, cfg: Config, *, log: str = "") -> dict:
    """`umi_tools dedup` on a coordinate-sorted+indexed BAM whose read names end
    in `_<UMI>`. Returns {'n_in', 'n_out', 'frac_kept', 'method'}.

    `out_bam` keeps the input's coordinate sort but is left UNINDEXED on purpose:
    it is the caller's temporary, and the caller indexes the final BAM.

    Raises ValueError if the reads carry no UMI (see the module docstring) or if
    the input BAM is not coordinate-sorted and indexed.
    """
    require_tools("umi_tools", "samtools")
    if not nonempty(bam):
        raise ValueError(f"no BAM to deduplicate: {bam}")
    method = str(cfg.get("process.umi_dedup_method", "directional"))

    _require_sorted_indexed(bam)
    _require_umis(bam)

    os.makedirs(os.path.dirname(os.path.abspath(out_bam)) or ".", exist_ok=True)
    cmd = [
        "umi_tools", "dedup",
        "-I", bam,
        "-S", out_bam,
        "--method", method,
        "--extract-umi-method", "read_id",     # the UMI is in the name, not a tag
        "--umi-separator", UMI_SEPARATOR,
    ]
    if log:
        cmd += ["-L", log]
    LOG.info("umi_tools dedup (--method %s): %s", method, os.path.basename(bam))
    run(cmd, log_to=log or None)
    if not nonempty(out_bam):
        raise RuntimeError(f"umi_tools dedup produced no output: {out_bam}")

    # Count both BAMs with samtools rather than scrape the umi_tools log: the log's
    # wording is version-dependent, a scrape breaks silently (and yields a
    # nonsensical frac_kept), whereas `samtools view -c` is exact.
    # No `samtools index out_bam` here -- see the module docstring: the caller
    # renames out_bam over the final BAM and indexes that.
    n_in = _count(bam)
    n_out = _count(out_bam)
    stats = {
        "n_in": n_in,
        "n_out": n_out,
        "frac_kept": (n_out / n_in) if n_in else 0.0,
        "method": method,
    }
    LOG.info("dedup: %d -> %d reads (%.1f%% kept)",
             n_in, n_out, 100 * stats["frac_kept"])
    return stats


def _require_sorted_indexed(bam: str) -> None:
    with pysam.AlignmentFile(bam, "rb") as fh:
        so = (fh.header.to_dict().get("HD") or {}).get("SO", "")
        if so != "coordinate":
            raise ValueError(
                f"{bam} is not coordinate-sorted (@HD SO:{so or 'unsorted'}); "
                f"umi_tools dedup needs a sorted+indexed BAM -- run star.sort_index first"
            )
        if not fh.has_index():
            raise ValueError(f"{bam} has no index -- run samtools index (or star.sort_index)")


def _require_umis(bam: str) -> None:
    """Refuse to deduplicate reads that carry no UMI.

    A read name ends in `_<UMI>` only if the architecture found random-templated
    content; with `dedup_umi_len == 0` there is nothing to deduplicate *on*, and
    umi_tools would collapse every read sharing a start position into one.
    """
    n_seen = n_umi = 0
    with pysam.AlignmentFile(bam, "rb") as fh:
        for aln in fh.head(UMI_HEAD_READS):
            n_seen += 1
            if _umi_of(aln.query_name):
                n_umi += 1
    if not n_seen:
        raise ValueError(f"{bam} contains no reads")
    frac = n_umi / n_seen
    if frac < UMI_TAG_MIN_FRAC:
        raise ValueError(
            f"{os.path.basename(bam)}: only {frac:.0%} of the first {n_seen} reads carry a "
            f"`_<UMI>` read-name tag -- this library has no random-templated content "
            f"(the architecture found dedup_umi_len == 0), so there is nothing to "
            f"deduplicate on. umi_tools would collapse genuine duplicate footprints "
            f"into one read. Set process.umi_dedup = false for this sample."
        )
    LOG.debug("UMI tag present on %.0f%% of the first %d reads", 100 * frac, n_seen)


def _umi_of(name: str) -> str:
    """The `_<UMI>` suffix of a read name, or "" if the read carries none."""
    head, sep, tail = name.rpartition(UMI_SEPARATOR)
    if not sep or not head or not tail or len(tail) > MAX_UMI_TAG:
        return ""
    return tail if set(tail.upper()) <= UMI_BASES else ""


def _count(bam: str) -> int:
    p = run(["samtools", "view", "-c", bam])
    return int((p.stdout or "0").strip() or 0)
