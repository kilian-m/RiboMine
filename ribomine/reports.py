"""The cross-sample tables: one row per run, one file per end point.

Each stage already wrote its JSON into the per-sample directory (`utils.Sample`);
a report only joins them. A run that never reached a stage still gets a row, so
the table shows which runs fell over. `counts/gene_counts.tsv` is the genes x runs
matrix of STAR's own gene counts.

Where the numbers come from:

* `qc_summary.tsv` -- the QC stage's local alignment of the untrimmed read sample
  (soft-clipping the 5' construct gives the sharpest periodicity, but inflates
  multimapping);
* `mapping_summary.tsv` -- mapping numbers from the alignment of the trimmed
  reads, and periodicity measured again on the finished BAM. When the two
  periodicities disagree, the trim is the suspect.
"""
from __future__ import annotations

import logging
import os
from typing import Any

from .architecture import structure_string
from .config import Config
from .process import counts
from .utils import Sample, nonempty, read_json, read_tsv, write_tsv

LOG = logging.getLogger("ribomine.reports")

# --- column order: identity -> call -> evidence -> provenance ----------------
QC_COLUMNS = [
    "run_accession", "study_accession", "verdict", "is_riboseq", "verdict_reason",
    "n_reads_sampled", "n_reads_scored", "total_run_reads", "projected_usable_reads",
    "read_len_mode", "read_len_peak_frac", "periodicity_inframe", "periodicity_tvd",
    "cds_frac_of_genic", "start_codon_ratio", "top5p_locus_frac",
    # footprints are sense to the gene; well under 1 = a reverse-complemented deposit
    "cds_sense_frac", "antisense_deposit",
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
    # the headline: is the run usable?
    "run_accession", "verdict", "architecture",
    "n_mapped",              # reads in the final BAM
    "mean_footprint_len",    # trimmed, contaminant-free reads, as fed to STAR
    "mean_mapped_len",       # aligned length of what reached the BAM
    "mean_mapped_softclip",  # large = residue the trim missed
    "periodicity_tvd",
    # --- the rest of the BAM's own measurement
    "periodicity_inframe", "read_len_mode", "n_cds_reads", "cds_frac_of_genic",
    "n_reads_scored_periodicity",
    # --- where the reads came from and what was thrown away on the way
    "download_route", "download_mb_per_s", "fastq_bytes",
    "n_reads_raw", "n_reads_after_trim", "frac_trimmed_out", "frac_no_adapter",
    # mean_len_after_trim still includes the contaminants; mean_footprint_len does not
    "mean_len_before_trim", "mean_len_after_trim",
    "n_contaminant_removed", "frac_contaminant",
    "n_reads_into_mapping", "n_uniquely_mapped",
    "frac_uniquely_mapped", "frac_multimapping", "frac_unmapped",
    "umi_dedup", "n_reads_after_dedup", "frac_duplicates",
    # --- STAR's gene counts (taken during alignment: before pile-up filter and dedup)
    "n_reads_in_genes", "frac_reads_in_genes", "n_genes_detected",
    "n_ambiguous", "n_no_feature", "sense_over_antisense",
    "bam", "bam_bytes",
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

# a stage that never produced its JSON
NOT_RUN = "NOT RUN"


# --- the three tables --------------------------------------------------------
def qc_tsv(cfg: Config, accs: list[str]) -> str:
    """Write `<workdir>/qc/qc_summary.tsv`."""
    meta = _candidates(cfg)
    fails = _failures(cfg)
    rows = [_qc_row(cfg, acc, meta.get(acc, {}), fails.get(acc, "")) for acc in accs]
    path = os.path.join(cfg.dir("qc"), "qc_summary.tsv")
    write_tsv(path, rows, columns=QC_COLUMNS)
    n_verdict = sum(1 for r in rows if r["verdict"] != NOT_RUN)
    LOG.info("qc_summary.tsv: %d runs (%d with a verdict) -> %s", len(rows), n_verdict, path)
    return path


def arch_tsv(cfg: Config, accs: list[str]) -> str:
    """Write `<workdir>/architecture/architecture.tsv`: the read layout of every run."""
    fails = _failures(cfg)
    rows = [_arch_row(cfg, acc, fails.get(acc, "")) for acc in accs]
    path = os.path.join(cfg.dir("architecture"), "architecture.tsv")
    write_tsv(path, rows, columns=ARCH_COLUMNS)
    n_ok = sum(1 for r in rows if r["status"] == "ok")
    LOG.info("architecture.tsv: %d runs (%d called) -> %s", len(rows), n_ok, path)
    return path


def process_tsv(cfg: Config, accs: list[str]) -> str:
    """Write `<workdir>/mapping_summary.tsv`: from bytes downloaded to reads in the BAM."""
    rows = [_process_row(cfg, acc) for acc in accs]
    path = os.path.join(cfg.workdir, "mapping_summary.tsv")
    write_tsv(path, rows, columns=PROCESS_COLUMNS)
    n_bam = sum(1 for r in rows if r.get("bam"))
    LOG.info("mapping_summary.tsv: %d runs (%d BAMs) -> %s", len(rows), n_bam, path)
    return path


def counts_tsv(cfg: Config, accs: list[str]) -> str:
    """Write `<workdir>/counts/gene_counts.tsv`, the genes x runs read-count matrix,
    from STAR's per-run `ReadsPerGene.out.tab`. Returns "" when no run has counts
    (an index built without a GTF cannot count)."""
    if not accs:
        return ""
    cols = [(a, counts.path(Sample(a, cfg.workdir).star_final)) for a in accs]
    have = [(a, p) for a, p in cols if nonempty(p)]
    if not have:
        LOG.warning("no run produced gene counts -- no read-count matrix. (STAR only "
                    "counts reads into genes when its index was built with the GTF: "
                    "--sjdbGTFfile at genome-generate.)")
        return ""
    if len(have) < len(cols):
        missing = [a for a, p in cols if not nonempty(p)]
        LOG.warning("read-count matrix: no gene counts for %d run(s) -- they are not "
                    "columns in it: %s", len(missing), ", ".join(missing[:5]))
    path = os.path.join(cfg.dir("counts"), "gene_counts.tsv")
    return counts.matrix(path, cfg.ref("star_index"), have)


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
        out.append(os.path.join(cfg.workdir, "counts", "gene_counts.tsv"))
    return [p for p in out if os.path.exists(p)]


# --- one row per table -------------------------------------------------------
def _qc_row(cfg: Config, acc: str, meta: dict, failure: str) -> dict[str, Any]:
    s = Sample(acc, cfg.workdir)
    qc = read_json(s.qc_json) or {}
    contam = read_json(s.contam_json) or {}
    pileup = read_json(s.pileup_json) or {}

    # the run's total read count: from the QC projection, else the ENA metadata
    total = _some(_get(qc, "usable.total_dataset_reads"), _int(meta.get("read_count")))
    row: dict[str, Any] = {
        "run_accession": acc,
        "study_accession": meta.get("study_accession", ""),
        "total_run_reads": total,
    }
    if not qc:
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
        # the QC stage may have run before the ENA read count was known
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
        "cds_sense_frac": qc.get("cds_sense_frac"),
        "antisense_deposit": qc.get("antisense_deposit"),
        "start_codon_ratio": qc.get("start_codon_ratio"),
        "top5p_locus_frac": qc.get("top5p_locus_frac"),
        # mitoribosome profiling: the numbers above are from the nuclear reads
        "mito_dominated": qc.get("mito_dominated"),
        "mito_periodicity_inframe": qc.get("mito_periodicity_inframe_frac"),
        "n_mito_cds_reads": qc.get("n_mito_cds_reads"),
        "frac_rRNA_tRNA_etc": _some(cm.get("frac_rRNA_tRNA_etc"),
                                    contam.get("frac_contaminant_structured_rna")),
        "frac_low_complexity": _some(cm.get("frac_low_complexity"),
                                     contam.get("frac_low_complexity")),
        "frac_position_pileup": frac_pileup,
        # local, permissive alignment: multimapping is inflated here
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
        "trim_5p": fn.get("trim_5p"),
        "umi5_len": call.get("umi5_len"),
        # the RT addition is kept with the footprint, so it is not part of trim_5p
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

    # measured on the deliverable BAM, and on STAR's gene counts
    per = proc.get("periodicity") or {}
    ct = proc.get("counts") or {}

    # a BAM deleted on purpose (`keep.bam: false`) is not a BAM that was never made
    bam = s.bam if nonempty(s.bam) else ""
    bam_bytes = os.path.getsize(bam) if bam else proc.get("bam_bytes")
    bam_cell = _rel(cfg, bam)
    if not bam and proc.get("keep", {}).get("bam") is False:
        bam_cell = "deleted (keep.bam=false)"

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
        # reads dropped because the 3' adapter ran off the end of the read
        "frac_no_adapter": trim.get("frac_no_adapter"),
        "mean_len_before_trim": trim.get("mean_len_in"),
        "mean_len_after_trim": trim.get("mean_len_out"),
        "n_contaminant_removed": n_contam,
        "frac_contaminant": frac_contam,
        "n_reads_into_mapping": n_map_in,
        "n_uniquely_mapped": mp.get("n_unique"),
        # alignment of the trimmed reads: these three sum to 1
        "frac_uniquely_mapped": mp.get("frac_unique"),
        "frac_multimapping": mp.get("frac_multimapping"),
        "frac_unmapped": mp.get("frac_unmapped"),
        "umi_dedup": proc.get("umi_dedup") if proc else None,
        "n_reads_after_dedup": dd.get("n_out"),
        "frac_duplicates": frac_dup,
        # --- measured on the finished BAM, after every filter
        "n_mapped": per.get("n_reads_in_bam"),
        "mean_footprint_len": mp.get("avg_input_len"),
        "mean_mapped_len": per.get("mean_mapped_len"),
        "mean_mapped_softclip": per.get("mean_mapped_softclip"),
        "read_len_mode": per.get("read_len_mode"),
        "periodicity_inframe": per.get("periodicity_inframe_frac"),
        "periodicity_tvd": per.get("periodicity_tvd_uniform"),
        "n_cds_reads": per.get("n_cds_reads"),
        "cds_frac_of_genic": per.get("cds_frac_of_genic"),
        "n_reads_scored_periodicity": per.get("n_reads_scored"),
        # --- STAR's gene counts for this run
        "n_reads_in_genes": ct.get("n_in_genes"),
        "frac_reads_in_genes": ct.get("frac_in_genes"),
        "n_genes_detected": ct.get("n_genes_detected"),
        "n_ambiguous": ct.get("n_ambiguous"),
        "n_no_feature": ct.get("n_no_feature"),
        "sense_over_antisense": ct.get("sense_over_antisense"),
        "bam": bam_cell,
        "bam_bytes": bam_bytes,
    }


# --- joins and lookups -------------------------------------------------------
def _candidates(cfg: Config) -> dict[str, dict]:
    """run_accession -> its row in `<workdir>/meta/candidates.tsv` (empty for runs
    that came from an accession list or a local FASTQ)."""
    path = os.path.join(cfg.dir("meta"), "candidates.tsv")
    if not nonempty(path):
        return {}
    try:
        return {r["run_accession"]: r for r in read_tsv(path) if r.get("run_accession")}
    except (OSError, KeyError) as exc:
        LOG.warning("could not read %s: %s", path, exc)
        return {}


def _failures(cfg: Config) -> dict[str, str]:
    """run_accession -> why it failed, from `<workdir>/failed.tsv`. Best-effort."""
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
    return structure_string(call) if call else ""


def _block(d: dict, *names: str, key: str) -> dict:
    """The sub-dict one stage wrote into the process JSON, located by a key only
    that block owns."""
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
    """The first value that is not None (a measured 0.0 counts)."""
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
    """Paths in the tables are relative to the workdir."""
    return os.path.relpath(path, cfg.workdir) if path else ""
