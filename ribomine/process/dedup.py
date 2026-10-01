"""UMI deduplication of the mapped reads (UMICollapse or umi_tools).

Trimming (fqdissect, with cutadapt) writes the UMI into the read name as
`_<UMI>`. This stage collapses reads that share a mapping position and a UMI,
i.e. PCR copies of one molecule. It is off by default (`process.umi_dedup`) and
checks two preconditions:

1. The reads carry a UMI. Without one (`dedup_umi_len == 0`), independent
   footprints at the same position cannot be told from PCR duplicates and would
   be collapsed, so `dedup` raises.
2. The input BAM is coordinate-sorted and indexed, as both tools require.

`process.umi_dedup_tool` picks the backend; both implement the directional,
adjacency and cluster methods. umicollapse (default) indexes the UMIs at a
position in a BK-tree instead of comparing all pairs, which suits ribo-seq,
where a start codon can stack hundreds of distinct UMIs on one coordinate.
umi_tools additionally offers `unique` and `percentile`.

The output is left unindexed: the caller (`pipeline.process_sample`) renames it
over the final BAM and indexes that.
"""
from __future__ import annotations

import os

import pysam

from ..config import Config
from ..utils import LOG, nonempty, require_tools, run

# Reads inspected to decide whether a BAM carries UMIs. The trimmer writes the
# suffix for every read or for none, so the head of the file is enough.
UMI_HEAD_READS = 10_000
UMI_TAG_MIN_FRAC = 0.5  # minimum fraction of those reads with a UMI suffix
UMI_BASES = set("ACGTN")
MAX_UMI_TAG = 24        # nt; a longer "_<suffix>" is not taken for a UMI
UMI_SEPARATOR = "_"

TOOLS = ("umicollapse", "umi_tools")

# `process.umi_dedup_method` uses umi_tools' method names for both backends; this
# maps them to UMICollapse's. UMICollapse has no `unique` or `percentile`.
UMICOLLAPSE_ALGO = {"directional": "dir", "adjacency": "adj", "cluster": "cc"}


def dedup(bam: str, out_bam: str, cfg: Config, *, log: str = "") -> dict:
    """Deduplicate a coordinate-sorted, indexed BAM whose read names end in `_<UMI>`.

    Returns {'n_in', 'n_out', 'frac_kept', 'method', 'tool'}. `out_bam` keeps the
    coordinate order and is left unindexed; the caller indexes the final BAM.

    Raises ValueError if the reads carry no UMI, if the input BAM is not
    coordinate-sorted and indexed, or if the tool does not implement the method.
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

    # Count with samtools rather than parse the tools' logs, whose wording depends
    # on the version.
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
    # --two-pass streams the BAM in coordinate order instead of holding it in
    # memory; it needs sorted input and writes sorted output.
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
    """Environment for UMICollapse's bioconda launcher.

    The launcher forwards only `-Xm*` arguments to the JVM, so the larger thread
    stack (`-Xss`) that the BK-tree recursion needs has to go through
    `_JAVA_OPTIONS`, which also replaces the launcher's fixed 4 GB heap
    (`process.umi_dedup_mem_gb`). `TEMP` is unset because the launcher otherwise
    appends `-log` and `-temp_folder` arguments that the jar rejects.
    """
    mem = int(cfg.get("process.umi_dedup_mem_gb", 8))
    env = dict(os.environ)
    env["_JAVA_OPTIONS"] = f"-Xms1g -Xmx{mem}g -Xss64m"
    env.pop("TEMP", None)
    return env


def _umi_tools(bam: str, out_bam: str, method: str, *, log: str) -> None:
    """umi_tools dedup, on a BAM whose read names carry the UMI."""
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
    """Raise unless the first reads carry a `_<UMI>` read-name suffix.

    Without a UMI (`dedup_umi_len == 0`) the deduplicator would collapse every
    read sharing a start position into one.
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
