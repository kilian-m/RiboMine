"""Read architecture, delegated to fqdissect.

fqdissect (https://github.com/kilian-m/fqdissect) reads the structure of a read --
UMIs, barcodes, RT additions, poly(A) tail, adapter -- out of a local alignment of
the untrimmed reads, and applies it with cutadapt. RiboMine hands it the
contaminant- and pile-up-filtered QC alignment and keeps the results in its own
per-sample layout. This module is the only place fqdissect is called from.
"""
from __future__ import annotations

import logging
from dataclasses import fields

from fqdissect import infer, plot, profile, structure, trim

from .config import Config
from .utils import read_json

# fqdissect's progress lines carry no sample name; in a pool they cannot be told apart
logging.getLogger("fqdissect").setLevel(logging.WARNING)


def profile_bam(bam: str, cfg: Config, *, label: str) -> dict:
    """Positional match/composition statistics of a locally-aligned BAM."""
    return profile.profile_bam(bam, cfg.ref("genome_fasta"), label=label,
                               max_reads=cfg["qc.sample_reads"],
                               min_anchor_frac=cfg["architecture.min_anchor_frac"])


def call(prof: dict, cfg: Config) -> dict:
    """The architecture call for one profile: `status` is 'ok' or 'undetermined'."""
    thr = infer.Thresholds(**{f.name: cfg[f"architecture.{f.name}"]
                              for f in fields(infer.Thresholds)})
    out = infer.infer(prof, thr)
    out["structure"] = structure.structure_string(out, prof)
    out["segments"] = structure.segments(out, prof)
    return out


def structure_string(call: dict) -> str:
    """`5'-[UMI,2nt]-[footprint,~30nt]-[UMI,5nt]-[barcode,AGCTA]-[TruSeq,...]-3'`"""
    return call.get("structure") or structure.structure_string(call)


def plot_call(prof: dict, call: dict, out_path: str, *, dpi: int = 110) -> str:
    return plot.plot_structure(prof, call, out_path, dpi=dpi)


def trim_fastq(call: dict, in_fq: str, out_fq: str, cfg: Config) -> dict:
    """Trim `in_fq` to the footprint with cutadapt; UMIs go to the read name
    (`@name_UMI`). Returns fqdissect's statistics plus the read lengths and the
    share of reads dropped for not showing the 3' adapter."""
    kw = dict(min_len=cfg["process.min_len"], threads=cfg["project.threads"],
              min_overlap=cfg["process.adapter_min_overlap"],
              discard_untrimmed=cfg["process.discard_untrimmed"])
    stats = trim.trim_fastq(call, in_fq, out_fq, **kw)

    # cutadapt's own per-pass reports hold the numbers fqdissect does not summarise
    stages, _ = trim.build_pipeline(call, in_fq, out_fq, **kw)
    passes = [read_json(s["json"], {}) for s in stages]
    n_in, n_out = stats["n_reads_in"], stats["n_reads_out"]
    if n_in:
        if kw["discard_untrimmed"]:
            no_adapter = sum(p.get("read_counts", {}).get("filtered", {})
                             .get("discard_untrimmed") or 0 for p in passes)
            stats["frac_no_adapter"] = round(no_adapter / n_in, 4)
        bp_in = passes[0].get("basepair_counts", {}).get("input")
        stats["mean_len_in"] = round(bp_in / n_in, 1) if bp_in else None
    if n_out:
        bp_out = passes[-1].get("basepair_counts", {}).get("output")
        stats["mean_len_out"] = round(bp_out / n_out, 1) if bp_out else None
    return stats
