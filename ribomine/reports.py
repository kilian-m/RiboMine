"""The cross-sample tables: one row per run, one file per end point.

For two of the three end points (`qc`, `architecture`) the TSV *is* the product,
so the columns are chosen the way a biologist reads them -- identity first, then
the verdict/call, then the numbers that back it -- and every column is a scalar
you can sort on. Nothing here re-derives biology: each stage already wrote its
JSON into the per-sample directory (`utils.Sample`), and a report only joins
them.

A run that never reached a stage still gets a row. A table with 187 rows and 12
blanks is honest; a table with 175 rows silently hides which runs fell over.

Where the numbers come from (this matters, and the pipeline is deliberate about
it -- see docs / PIPELINE.md §13):

* periodicity, region composition and the mapping fractions in `qc_summary.tsv`
  come from the QC stage's **local** alignment of the untrimmed sampled reads.
  Soft-clipping the 5' construct / RT base gives a cleaner footprint 5' end, and
  hence sharper periodicity, than an end-to-end alignment of trimmed reads. The
  same permissiveness inflates multimapping, which is why
* the mapping numbers in `mapping_summary.tsv` come from the **end-to-end**
  alignment of the trimmed reads -- the deliverable BAM.
"""
from __future__ import annotations

import logging
import os
from typing import Any

from .arch.infer import architecture_string
from .config import Config
from .utils import Sample, nonempty, read_json, read_tsv, write_tsv

LOG = logging.getLogger("ribomine.reports")

# --- column order: identity -> call -> evidence -> provenance ----------------
QC_COLUMNS = [
    "run_accession", "study_accession", "verdict", "is_riboseq", "verdict_reason",
    "n_reads_sampled", "n_reads_scored", "total_run_reads", "projected_usable_reads",
    "read_len_mode", "read_len_peak_frac", "periodicity_inframe", "periodicity_tvd",
    "cds_frac_of_genic", "start_codon_ratio", "top5p_locus_frac",
    "mito_dominated", "mito_periodicity_inframe", "n_mito_cds_reads",
    "frac_rRNA_tRNA_etc", "frac_low_complexity", "frac_position_pileup",
    "frac_uniquely_mapped", "frac_multimapping", "frac_unmapped",
    "region_CDS", "region_5UTR", "region_3UTR", "region_ncRNA", "region_intron",
    "region_intergenic", "region_mito",
    "qc_plot",
]

ARCH_COLUMNS = [
    "run_accession", "status", "reason", "architecture", "deposit_state",
    "trim_5p", "umi5_len", "rt_len", "rt_penetrance", "barcode5_seq",
    "template_switch5_seq", "umi3_len", "barcode3_seq", "adapter_name",
    "adapter_seq", "polyA_tail", "trim_3p_construct", "dedup_umi_len",
    "dedup_umi_mask", "footprint_len_mode", "read_len_mode", "genomic_plateau",
    "frac_anchored", "n_reads_used", "flags", "arch_plot",
]

PROCESS_COLUMNS = [
    "run_accession", "verdict", "architecture", "download_route", "download_mb_per_s",
    "fastq_bytes", "n_reads_raw", "n_reads_after_trim", "frac_trimmed_out", "frac_no_adapter",
    "mean_len_before_trim", "mean_len_after_trim", "n_contaminant_removed",
    "frac_contaminant", "n_reads_into_mapping", "n_uniquely_mapped",
    "frac_uniquely_mapped", "frac_multimapping", "frac_unmapped", "umi_dedup",
    "n_reads_after_dedup", "frac_duplicates", "bam", "bam_bytes",
]

# region_frac keys as the annotation index emits them -> flat column names
REGIONS = [
    ("region_CDS", "CDS"),
    ("region_5UTR", "5'UTR"),
    ("region_3UTR", "3'UTR"),
    ("region_ncRNA", "ncRNA"),
    ("region_intron", "intron"),
    ("region_intergenic", "intergenic"),
    ("region_mito", "mito"),
]

# a stage that never produced its JSON: said out loud, not left blank
NOT_RUN = "NOT RUN"


# --- the three tables --------------------------------------------------------
def qc_tsv(cfg: Config, accs: list[str]) -> str:
    """Write `<workdir>/qc/qc_summary.tsv` -- the QC end point's deliverable."""
    meta = _candidates(cfg)
    fails = _failures(cfg)
    rows = [_qc_row(cfg, acc, meta.get(acc, {}), fails.get(acc, "")) for acc in accs]
    path = os.path.join(cfg.dir("qc"), "qc_summary.tsv")
    write_tsv(path, rows, columns=QC_COLUMNS)
    n_verdict = sum(1 for r in rows if r["verdict"] != NOT_RUN)
    LOG.info("qc_summary.tsv: %d runs (%d with a verdict) -> %s", len(rows), n_verdict, path)
    return path


def arch_tsv(cfg: Config, accs: list[str]) -> str:
    """Write `<workdir>/architecture/architecture.tsv` -- the architecture end
    point's deliverable: the read layout of every run, one line each."""
    fails = _failures(cfg)
    rows = [_arch_row(cfg, acc, fails.get(acc, "")) for acc in accs]
    path = os.path.join(cfg.dir("architecture"), "architecture.tsv")
    write_tsv(path, rows, columns=ARCH_COLUMNS)
    n_ok = sum(1 for r in rows if r["status"] == "ok")
    LOG.info("architecture.tsv: %d runs (%d called) -> %s", len(rows), n_ok, path)
    return path


def process_tsv(cfg: Config, accs: list[str]) -> str:
    """Write `<workdir>/mapping_summary.tsv` -- the BAM end point's deliverable:
    what came out of every run, from bytes downloaded to reads in the BAM."""
    rows = [_process_row(cfg, acc) for acc in accs]
    path = os.path.join(cfg.workdir, "mapping_summary.tsv")
    write_tsv(path, rows, columns=PROCESS_COLUMNS)
    n_bam = sum(1 for r in rows if r.get("bam"))
    LOG.info("mapping_summary.tsv: %d runs (%d BAMs) -> %s", len(rows), n_bam, path)
    return path


def summary_line(cfg: Config, accs: list[str]) -> str:
    """The paragraph the CLI prints when the run finishes."""
    keep = set(cfg.get("pipeline.keep_verdicts") or [])
    verdicts: dict[str, int] = {}
    n_pass = n_ok = n_undet = n_arch_missing = n_bam = 0
    for acc in accs:
        s = Sample(acc, cfg.workdir)
        qc = read_json(s.qc_json) or {}
        v = qc.get("verdict", NOT_RUN)
        verdicts[v] = verdicts.get(v, 0) + 1
        if v in keep:
            n_pass += 1
        call = read_json(s.arch_json) or {}
        status = call.get("status")
        if status == "ok":
            n_ok += 1
        elif status:
            n_undet += 1
        else:
            n_arch_missing += 1
        if nonempty(s.bam):
            n_bam += 1

    order = ["RIBO-SEQ", "TI-SEQ", "NOT RIBO-SEQ or LOW QUALITY", NOT_RUN]
    breakdown = ", ".join(
        f"{verdicts[v]} {v}" for v in order + sorted(set(verdicts) - set(order))
        if verdicts.get(v)
    )
    end = cfg.get("pipeline.end", "bam")
    parts = [
        f"{len(accs)} run(s), {n_pass} passed QC.",
        f"Verdicts: {breakdown or 'none'}.",
    ]
    if end in ("architecture", "bam"):
        parts.append(
            f"Architecture called for {n_ok}, undetermined for {n_undet}"
            + (f", not attempted for {n_arch_missing}" if n_arch_missing else "")
            + "."
        )
    if end == "bam":
        parts.append(f"{n_bam} BAM(s) written.")
    parts.append("Tables: " + ", ".join(_tables(cfg, end)) + ".")
    return " ".join(parts)


def _tables(cfg: Config, end: str) -> list[str]:
    """The deliverables that exist, in end-point order."""
    out = [os.path.join(cfg.workdir, "qc", "qc_summary.tsv")]
    if end in ("architecture", "bam"):
        out.append(os.path.join(cfg.workdir, "architecture", "architecture.tsv"))
    if end == "bam":
        out.append(os.path.join(cfg.workdir, "mapping_summary.tsv"))
    return [p for p in out if os.path.exists(p)]


# --- one row per table -------------------------------------------------------
def _qc_row(cfg: Config, acc: str, meta: dict, failure: str) -> dict[str, Any]:
    s = Sample(acc, cfg.workdir)
    qc = read_json(s.qc_json) or {}
    contam = read_json(s.contam_json) or {}
    pileup = read_json(s.pileup_json) or {}

    # the run's total read count: from the QC projection if the stage had it,
    # else straight from the ENA metadata we queried with
    total = _some(_get(qc, "usable.total_dataset_reads"), _int(meta.get("read_count")))
    row: dict[str, Any] = {
        "run_accession": acc,
        "study_accession": meta.get("study_accession", ""),
        "total_run_reads": total,
    }
    if not qc:
        # no verdict: say so in the verdict column rather than leaving a blank
        # row that reads like a pass
        row["verdict"] = NOT_RUN
        row["verdict_reason"] = failure or "QC produced no output for this run"
        row["n_reads_sampled"] = contam.get("n_input")
        return row

    cm = qc.get("contaminants") or {}
    mapping = qc.get("mapping") or {}
    reg = qc.get("region_frac") or {}

    n_sampled = _some(cm.get("n_sampled"), _get(qc, "usable.n_sampled"),
                      contam.get("n_input"))
    n_scored = qc.get("n_reads_scored")
    projected = _get(qc, "usable.projected_usable_reads")
    if projected is None and total and n_sampled and n_scored is not None:
        # the QC stage runs before the ENA read count is necessarily known; if we
        # have both now, project here rather than leave the column empty
        projected = int(round(n_scored / n_sampled * total))

    frac_pileup = cm.get("frac_position_pileup")
    if frac_pileup is None and n_sampled and pileup.get("n_reads_removed") is not None:
        frac_pileup = round(pileup["n_reads_removed"] / n_sampled, 4)

    row.update({
        "verdict": qc.get("verdict", NOT_RUN),
        "is_riboseq": qc.get("is_riboseq"),
        "verdict_reason": qc.get("verdict_reason", ""),
        "n_reads_sampled": n_sampled,
        "n_reads_scored": n_scored,
        "projected_usable_reads": projected,
        "read_len_mode": qc.get("read_len_mode"),
        "read_len_peak_frac": qc.get("read_len_peak_frac"),
        "periodicity_inframe": qc.get("periodicity_inframe_frac"),
        "periodicity_tvd": qc.get("periodicity_tvd_uniform"),
        "cds_frac_of_genic": qc.get("cds_frac_of_genic"),
        "start_codon_ratio": qc.get("start_codon_ratio"),
        "top5p_locus_frac": qc.get("top5p_locus_frac"),
        # a mitoribosome-profiling library: the numbers above were measured on its
        # NUCLEAR reads, which are the minority. These say what its mito reads do.
        "mito_dominated": qc.get("mito_dominated"),
        "mito_periodicity_inframe": qc.get("mito_periodicity_inframe_frac"),
        "n_mito_cds_reads": qc.get("n_mito_cds_reads"),
        "frac_rRNA_tRNA_etc": _some(cm.get("frac_rRNA_tRNA_etc"),
                                    contam.get("frac_contaminant_structured_rna")),
        "frac_low_complexity": _some(cm.get("frac_low_complexity"),
                                     contam.get("frac_low_complexity")),
        "frac_position_pileup": frac_pileup,
        # local, permissive alignment -- multimapping is inflated by design here
        "frac_uniquely_mapped": mapping.get("frac_unique"),
        "frac_multimapping": mapping.get("frac_multimapping"),
        "frac_unmapped": mapping.get("frac_unmapped"),
        "qc_plot": _rel(cfg, _plot(cfg, s.qc_plot)),
    })
    for col, key in REGIONS:
        row[col] = reg.get(key)
    return row


def _arch_row(cfg: Config, acc: str, failure: str) -> dict[str, Any]:
    s = Sample(acc, cfg.workdir)
    call = read_json(s.arch_json) or {}
    row: dict[str, Any] = {"run_accession": acc}
    if not call:
        qc = read_json(s.qc_json) or {}
        keep = set(cfg.get("pipeline.keep_verdicts") or [])
        verdict = qc.get("verdict")
        if verdict and verdict not in keep:
            # dropped on purpose between the stages -- not a failure
            row.update(status="not called", reason=f"QC verdict: {verdict}")
        elif failure:
            row.update(status="error", reason=failure)
        else:
            row.update(status=NOT_RUN,
                       reason="architecture produced no output for this run")
        return row

    fn = call.get("functional") or {}
    row.update({
        "status": call.get("status", ""),
        "reason": call.get("reason", ""),
        "architecture": _arch_string(call),
        "deposit_state": call.get("deposit_state"),
        # the trim plan: what actually comes off the read
        "trim_5p": fn.get("trim_5p"),
        "umi5_len": call.get("umi5_len"),
        # the enzymatic RT addition is KEPT (it is a footprint base), so it is
        # reported apart from the 5' overhead that is trimmed
        "rt_len": fn.get("footprint_retains_rt_nt"),
        "rt_penetrance": call.get("rt_penetrance"),
        "barcode5_seq": call.get("barcode5_seq"),
        "template_switch5_seq": call.get("template_switch5_seq"),
        "umi3_len": call.get("umi3_len"),
        "barcode3_seq": call.get("barcode3_seq"),
        "adapter_name": call.get("adapter3_name"),
        "adapter_seq": call.get("adapter3_seq"),
        "polyA_tail": call.get("polyA_tail"),
        "trim_3p_construct": fn.get("trim_3p_construct"),
        "dedup_umi_len": fn.get("dedup_umi_len"),
        "dedup_umi_mask": fn.get("dedup_umi_mask"),
        "footprint_len_mode": call.get("footprint_len_mode"),
        "read_len_mode": call.get("read_len_mode"),
        "genomic_plateau": call.get("genomic_plateau"),
        "frac_anchored": call.get("frac_anchored"),
        "n_reads_used": call.get("n_reads_used"),
        "flags": call.get("flags") or [],
        "arch_plot": _rel(cfg, _plot(cfg, s.arch_plot)),
    })
    return row


def _process_row(cfg: Config, acc: str) -> dict[str, Any]:
    s = Sample(acc, cfg.workdir)
    proc = read_json(s.process_json) or {}
    trim = read_json(s.trim_json) or {}
    qc = read_json(s.qc_json) or {}
    call = read_json(s.arch_json) or {}

    dl = _block(proc, "download", key="route")
    cm = _block(proc, "contaminants", "contam", key="n_contaminant_rRNA_tRNA_etc")
    mp = _block(proc, "mapping", "star", key="n_unique")
    dd = _block(proc, "dedup", "umi_dedup", key="n_out")

    n_in = trim.get("n_reads_in")
    n_out = trim.get("n_reads_out")
    # what trimming threw away: reads left shorter than process.min_len once the
    # construct came off (an adapter-dimer-heavy library loses a lot here)
    frac_trimmed_out = round(1 - n_out / n_in, 4) if n_in and n_out is not None else None

    n_contam = None
    frac_contam = None
    if cm:
        n_contam = cm.get("n_input", 0) - cm.get("n_kept", 0)
        if cm.get("n_input"):
            frac_contam = round(n_contam / cm["n_input"], 4)

    n_map_in = _some(mp.get("n_input"), cm.get("n_kept"), n_out)

    frac_dup = None
    if dd.get("frac_kept") is not None:
        frac_dup = round(1 - dd["frac_kept"], 4)
    elif dd.get("n_in") and dd.get("n_out") is not None:
        frac_dup = round(1 - dd["n_out"] / dd["n_in"], 4)

    fastq_bytes = dl.get("bytes")
    if fastq_bytes is None and nonempty(s.full_fastq):
        fastq_bytes = os.path.getsize(s.full_fastq)

    bam = s.bam if nonempty(s.bam) else ""
    return {
        "run_accession": acc,
        "verdict": qc.get("verdict", ""),
        "architecture": _arch_string(call),
        "download_route": dl.get("route"),
        "download_mb_per_s": dl.get("mb_per_s"),
        "fastq_bytes": fastq_bytes,
        "n_reads_raw": n_in,
        "n_reads_after_trim": n_out,
        "frac_trimmed_out": frac_trimmed_out,
        # reads whose fixed 3' scaffold ran off the end of the read: they cannot be cut
        # at the footprint boundary and are discarded. High here = the library's
        # molecules are as long as its reads, and most of its depth is unusable.
        "frac_no_adapter": trim.get("frac_no_adapter"),
        "mean_len_before_trim": trim.get("mean_len_in"),
        "mean_len_after_trim": trim.get("mean_len_out"),
        "n_contaminant_removed": n_contam,
        "frac_contaminant": frac_contam,
        "n_reads_into_mapping": n_map_in,
        "n_uniquely_mapped": mp.get("n_unique"),
        # end-to-end alignment of the trimmed reads: these three sum to 1
        "frac_uniquely_mapped": mp.get("frac_unique"),
        "frac_multimapping": mp.get("frac_multimapping"),
        "frac_unmapped": mp.get("frac_unmapped"),
        # the pipeline records whether dedup RAN; `dd` only says whether it produced
        # stats, which is not the same thing (dedup off => no block, but also no 'no')
        "umi_dedup": proc.get("umi_dedup") if proc else None,
        "n_reads_after_dedup": dd.get("n_out"),
        "frac_duplicates": frac_dup,
        "bam": _rel(cfg, bam),
        "bam_bytes": os.path.getsize(bam) if bam else None,
    }


# --- joins and lookups -------------------------------------------------------
def _candidates(cfg: Config) -> dict[str, dict]:
    """run_accession -> its row in `<workdir>/meta/candidates.tsv`.

    The study accession and the run's total read count live in the ENA metadata
    the query stage wrote. A run that came from a bare accession list or a local
    FASTQ has no such row -- its cells stay empty, which is the truth.
    """
    path = os.path.join(cfg.dir("meta"), "candidates.tsv")
    if not nonempty(path):
        return {}
    try:
        return {r["run_accession"]: r for r in read_tsv(path) if r.get("run_accession")}
    except (OSError, KeyError) as exc:
        LOG.warning("could not read %s: %s", path, exc)
        return {}


def _failures(cfg: Config) -> dict[str, str]:
    """run_accession -> why the pipeline recorded it as failed, from
    `<workdir>/failed.tsv`. Best-effort: the table is for humans, and a missing
    or oddly-shaped failure log must never take a report down."""
    path = os.path.join(cfg.workdir, "failed.tsv")
    if not nonempty(path):
        return {}
    out: dict[str, str] = {}
    try:
        for row in read_tsv(path):
            acc = _pick(row, "run_accession", "accession", "acc", "sample", "run")
            if not acc:
                continue
            err = _pick(row, "error", "reason", "message", "exception", "detail") or "failed"
            stage = _pick(row, "stage", "step")
            out[acc] = f"{stage}: {err}" if stage else err
    except OSError as exc:
        LOG.warning("could not read %s: %s", path, exc)
    return out


def _arch_string(call: dict) -> str:
    """The human-readable one-liner, e.g.
    `5'-[UMI,2nt]-[footprint,~30nt]-[UMI,5nt]-[barcode,AGCTA]-[TruSeq]-3'`."""
    if not call:
        return ""
    try:
        return architecture_string(call)
    except (KeyError, TypeError, ValueError) as exc:
        LOG.debug("no architecture string for %s: %s", call.get("run_accession"), exc)
        return f"{call.get('status', '?')}: {call.get('reason', '')}"


def _block(d: dict, *names: str, key: str) -> dict:
    """The sub-dict one stage wrote, wherever the process JSON nested it.

    The processing stage merges several tools' stats (download, contaminant
    filter, mapping, dedup) into one JSON. Locating a block by a key only that
    block owns keeps the report robust to how they are nested.
    """
    for name in names:
        blk = d.get(name)
        if isinstance(blk, dict) and key in blk:
            return blk
    if key in d:                       # flattened into the top level
        return d
    for v in d.values():               # nested under a name we did not guess
        if isinstance(v, dict) and key in v:
            return v
    return {}


def _get(d: dict, dotted: str, default=None):
    cur: Any = d
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


def _some(*vals):
    """The first value that was actually measured. Not `a or b`: a measured
    fraction of 0.0 (no low-complexity reads at all) is an answer, not a miss."""
    for v in vals:
        if v is not None:
            return v
    return None


def _pick(row: dict, *keys: str) -> str:
    for k in keys:
        v = row.get(k)
        if v:
            return str(v).strip()
    return ""


def _int(v) -> int | None:
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _plot(cfg: Config, plot: Any) -> str:
    """The sample's plot for the configured format, if it was drawn."""
    if not cfg.get("plots.enabled", True):
        return ""
    p = plot(cfg.get("plots.format", "png"))
    return p if nonempty(p) else ""


def _rel(cfg: Config, path: str) -> str:
    """Paths in the tables are relative to the workdir (`bams/SRR1.bam`): short
    enough to read, and the table survives the workdir being moved."""
    return os.path.relpath(path, cfg.workdir) if path else ""
