"""UMI deduplication of the mapped reads (UMICollapse, or umi_tools).

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
   that start at the same codon from a PCR duplicate pair. Running a deduplicator
   anyway would silently collapse genuine footprints into one read -- at ribo-seq
   depth, where a well-expressed CDS has hundreds of reads on the same start,
   that quietly destroys the signal it is supposed to clean. So: raise.
2. **The INPUT BAM must be coordinate-sorted and indexed** -- both tools walk it
   by position, and UMICollapse's `--two-pass` requires it.

Which tool
----------
`process.umi_dedup_tool` picks the backend; both implement the same algorithms
(directional / adjacency / cluster) over the same `_<UMI>` read-name tag, and on
a 101k-read test BAM they returned the *same* 81,182 reads.

* **umicollapse** (the default) indexes the UMIs at each position in an n-gram
  BK-tree rather than comparing every pair, so its cost at a position grows far
  more gently than umi_tools' all-pairs scan. That is exactly the shape of
  ribo-seq: a well-translated start codon piles hundreds of distinct UMIs onto
  one coordinate, which is umi_tools' worst case and UMICollapse's design case.
* **umi_tools** is kept because it is the reference implementation everyone knows,
  and because it has methods (`unique`, `percentile`) UMICollapse does not.

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

TOOLS = ("umicollapse", "umi_tools")

# umi_tools' method names are the vocabulary `process.umi_dedup_method` speaks, in
# both backends -- the algorithms are the same, only UMICollapse's names for them
# are shorter. It has no equivalent of `unique` or `percentile`.
UMICOLLAPSE_ALGO = {"directional": "dir", "adjacency": "adj", "cluster": "cc"}


def dedup(bam: str, out_bam: str, cfg: Config, *, log: str = "") -> dict:
    """Deduplicate a coordinate-sorted+indexed BAM whose read names end in `_<UMI>`.

    Returns {'n_in', 'n_out', 'frac_kept', 'method', 'tool'}.

    `out_bam` keeps the input's coordinate sort but is left UNINDEXED on purpose:
    it is the caller's temporary, and the caller indexes the final BAM.

    Raises ValueError if the reads carry no UMI (see the module docstring), if the
    input BAM is not coordinate-sorted and indexed, or if the configured tool and
    method are not a pair that exists.
    """
    if not nonempty(bam):
        raise ValueError(f"no BAM to deduplicate: {bam}")
    tool = str(cfg.get("process.umi_dedup_tool", "umicollapse"))
    method = str(cfg.get("process.umi_dedup_method", "directional"))
    if tool not in TOOLS:
        raise ValueError(f"process.umi_dedup_tool must be one of {TOOLS}, got {tool!r}")

    _require_sorted_indexed(bam)
    _require_umis(bam)
    os.makedirs(os.path.dirname(os.path.abspath(out_bam)) or ".", exist_ok=True)

    if tool == "umicollapse":
        _umicollapse(bam, out_bam, method, cfg, log=log)
    else:
        _umi_tools(bam, out_bam, method, log=log)
    if not nonempty(out_bam):
        raise RuntimeError(f"{tool} produced no output: {out_bam}")

    # Count both BAMs with samtools rather than scrape the tool's log: the wording is
    # version-dependent, a scrape breaks silently (and yields a nonsensical
    # frac_kept), whereas `samtools view -c` is exact.
    # No `samtools index out_bam` here -- see the module docstring: the caller
    # renames out_bam over the final BAM and indexes that.
    n_in = _count(bam)
    n_out = _count(out_bam)
    stats = {
        "n_in": n_in,
        "n_out": n_out,
        "frac_kept": (n_out / n_in) if n_in else 0.0,
        "method": method,
        "tool": tool,
    }
    LOG.info("dedup (%s, %s): %d -> %d reads (%.1f%% kept)",
             tool, method, n_in, n_out, 100 * stats["frac_kept"])
    return stats


def _umicollapse(bam: str, out_bam: str, method: str, cfg: Config, *, log: str) -> None:
    """UMICollapse, on a BAM whose read names carry the UMI."""
    require_tools("umicollapse", "samtools")
    algo = UMICOLLAPSE_ALGO.get(method)
    if algo is None:
        raise ValueError(
            f"process.umi_dedup_method={method!r} has no UMICollapse equivalent "
            f"(it implements {', '.join(sorted(UMICOLLAPSE_ALGO))}). Either pick one of "
            f"those, or set process.umi_dedup_tool='umi_tools', which has {method!r}."
        )
    # --two-pass streams the BAM in coordinate order instead of holding it in memory,
    # which is what keeps a deep run inside the heap -- and it is why the input has to
    # be sorted (it is; _require_sorted_indexed just checked). The output comes back in
    # coordinate order, so the caller does not have to re-sort it.
    cmd = [
        "umicollapse", "bam",
        "-i", bam,
        "-o", out_bam,
        "--umi-sep", UMI_SEPARATOR,     # the UMI is in the read name, not a tag
        "--algo", algo,
        "--two-pass",
    ]
    LOG.info("umicollapse (--algo %s): %s", algo, os.path.basename(bam))
    run(cmd, env=_jvm_env(cfg), log_to=log or None)


def _jvm_env(cfg: Config) -> dict:
    """The environment UMICollapse has to be run in, because its bioconda launcher
    will not otherwise do the right thing.

    That launcher forwards only `-Xm*` arguments to the JVM, so `-Xss` -- which
    UMICollapse's BK-tree recursion needs raised, and which its own README passes --
    cannot be given on the command line at all. `_JAVA_OPTIONS` is the only route in,
    and setting it also replaces the launcher's fixed 4 GB heap, which a deep BAM
    outgrows. And if `TEMP` happens to be set, the launcher silently appends `-log`
    and `-temp_folder` arguments that the jar does not accept and dies on -- so it
    must not be set for the child.
    """
    mem = int(cfg.get("process.umi_dedup_mem_gb", 8))
    env = dict(os.environ)
    env["_JAVA_OPTIONS"] = f"-Xms1g -Xmx{mem}g -Xss64m"
    env.pop("TEMP", None)
    return env


def _umi_tools(bam: str, out_bam: str, method: str, *, log: str) -> None:
    """umi_tools dedup -- the reference implementation, kept as the alternative."""
    require_tools("umi_tools", "samtools")
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


def _require_sorted_indexed(bam: str) -> None:
    with pysam.AlignmentFile(bam, "rb") as fh:
        so = (fh.header.to_dict().get("HD") or {}).get("SO", "")
        if so != "coordinate":
            raise ValueError(
                f"{bam} is not coordinate-sorted (@HD SO:{so or 'unsorted'}); "
                f"UMI deduplication needs a sorted+indexed BAM -- run star.sort_index first"
            )
        if not fh.has_index():
            raise ValueError(f"{bam} has no index -- run samtools index (or star.sort_index)")


def _require_umis(bam: str) -> None:
    """Refuse to deduplicate reads that carry no UMI.

    A read name ends in `_<UMI>` only if the architecture found random-templated
    content; with `dedup_umi_len == 0` there is nothing to deduplicate *on*, and
    the deduplicator would collapse every read sharing a start position into one.
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
            f"deduplicate on. Deduplicating anyway would collapse genuine duplicate "
            f"footprints into one read. Set process.umi_dedup = false for this sample."
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
