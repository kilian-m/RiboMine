"""Stage orchestration: what runs, in what order, and what a resume skips.

    query ──► qc ──► architecture ──► bam
      ▲        ▲
      │        └── start="accessions" (a list) or start="fastq" (a local directory)
      └── start="query"

A run stops after the stage named by `pipeline.end`. The QC stage samples,
filters, aligns (locally) and profiles the reads; the architecture stage reads
only the profile JSON, so it costs milliseconds and can be re-run over a cohort
without touching a BAM.

Per-sample work runs in a process pool. A sample that fails is recorded in
failed.tsv and does not take the run down.
"""
from __future__ import annotations

import os
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool

from . import architecture, reports
from .config import Config, stages_to_run
from .process import counts
from .process import dedup as proc_dedup
from .process import star
from .qc import annotation, contaminants, pileups, plot as qc_plot, verdict
from .sra import download, query
from .utils import (LOG, Sample, human, nonempty, read_json, read_lines, read_tsv,
                    rm, write_json, write_tsv)


# --- the sample list: the three start points -------------------------------
def resolve_inputs(cfg: Config) -> tuple[list[str], dict[str, str]]:
    """(accessions, source_map). `source_map[acc]` is what the reads come from:
    an accession (streamed from ENA) or a local FASTQ path."""
    start = cfg["pipeline.start"]

    if start == "fastq":
        d = cfg._abs(cfg["pipeline.fastq_dir"])
        exts = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
        src = {}
        for f in sorted(os.listdir(d)):
            ext = next((e for e in exts if f.endswith(e)), None)
            if ext:
                src[f[: -len(ext)]] = os.path.join(d, f)
        if not src:
            raise RuntimeError(f"no FASTQ files in {d}")
        LOG.info("start=fastq: %d local FASTQ file(s) in %s", len(src), d)
        return list(src), src

    if start == "accessions":
        accs = read_lines(cfg._abs(cfg["pipeline.accession_list"]))
        LOG.info("start=accessions: %d run accession(s)", len(accs))
        return accs, {a: a for a in accs}

    # start == "query": the candidates TSV drives everything downstream
    tsv = os.path.join(cfg.dir("meta"), "candidates.tsv")
    if cfg["pipeline.resume"] and nonempty(tsv):
        LOG.info("resume: reusing %s (delete it to re-query)", tsv)
    else:
        tsv = query.run_query(cfg)
    accs = [r["run_accession"] for r in read_tsv(tsv)]
    LOG.info("start=query: %d candidate run(s) from %s", len(accs), tsv)
    return accs, {a: a for a in accs}


def _ensure_fasta_index(cfg: Config) -> None:
    """Build the genome `.fai` before the pool starts. pysam builds a missing
    index in place, so concurrent workers would corrupt it."""
    fasta = cfg.ref("genome_fasta")
    if nonempty(fasta + ".fai"):
        return
    import pysam

    LOG.info("indexing the genome FASTA (once): %s", fasta)
    pysam.faidx(fasta)


def _run_meta(cfg: Config) -> dict[str, dict]:
    """run_accession -> its candidates.tsv row (read_count, study, ...), if any."""
    tsv = os.path.join(cfg.workdir, "meta", "candidates.tsv")
    if not nonempty(tsv):
        return {}
    return {r["run_accession"]: r for r in read_tsv(tsv)}


# --- stage qc: sample -> contaminant filter -> local alignment -> profile -> verdict
def qc_sample(acc: str, cfg: Config, src: str, meta: dict) -> dict:
    s = Sample(acc, cfg.workdir)
    resume = cfg["pipeline.resume"]
    threads = cfg["project.threads"]

    if resume and nonempty(s.qc_json) and nonempty(s.profile_json):
        LOG.info("[%s] qc: cached", acc)
        return read_json(s.qc_json)

    t0 = time.time()

    # 1. a uniform read sample, streamed -- nothing is downloaded in full here
    if not (resume and nonempty(s.sample_fastq)):
        st = download.sample_reads(src, s.sample_fastq, n=cfg["qc.sample_reads"],
                                   scan=cfg["qc.scan_reads"], seed=cfg["project.seed"])
        LOG.info("[%s] sampled %s reads of %s scanned", acc,
                 human(st["n_sampled"]), human(st["n_scanned"]))
        if st.get("sorted_warning"):
            LOG.warning("[%s] %s", acc, st["sorted_warning"])

    # 2. contaminants (rRNA/tRNA/snRNA/Mt, low complexity). The reads stay untrimmed:
    #    the architecture is read out of the soft clips.
    if not (resume and nonempty(s.filtered_fastq) and nonempty(s.contam_json)):
        cst = contaminants.filter_fastq(s.sample_fastq, s.filtered_fastq, cfg,
                                        label=acc, threads=threads, log=s.log)
        write_json(s.contam_json, cst)
    cst = read_json(s.contam_json, {})

    # 3. local alignment: everything non-genomic lands in the soft clips
    if not (resume and nonempty(s.local_bam)):
        star.align_local(s.filtered_fastq, s.star_local, cfg, threads=threads, log=s.log)

    # 4. position pile-ups (adapter dimers, fixed contaminants)
    if not (resume and nonempty(s.pileup_bam) and nonempty(s.pileup_json)):
        pst = pileups.filter_bam(s.local_bam, s.pileup_bam, cfg, label=acc)
        write_json(s.pileup_json, pst)
    pst = read_json(s.pileup_json, {})

    # 5. the positional profile the architecture is called from
    if not (resume and nonempty(s.profile_json)):
        write_json(s.profile_json, architecture.profile_bam(s.pileup_bam, cfg, label=acc))

    # 6. the verdict
    q = verdict.qc(s.pileup_bam, cfg.annotation_index, cfg, star_log=s.local_log,
                   label=acc, contam=cst, pileup=pst, total_reads=_total_reads(acc, meta))
    write_json(s.qc_json, q)
    if cfg["plots.enabled"]:
        qc_plot.plot_qc(q, s.qc_plot(cfg["plots.format"]), dpi=cfg["plots.dpi"])

    # 7. the working data. Every number is already in the JSONs, so the sample and
    #    its alignments go unless `keep` says otherwise; kept BAMs are sorted and
    #    indexed (only now -- the steps above read them in file order).
    if cfg["keep.qc_bam"]:
        for b in (s.local_bam, s.pileup_bam):
            if nonempty(b):
                star.ensure_sorted_indexed(b, threads=min(threads, 8))
    else:
        rm(s.local_bam, s.local_bam + ".bai", s.pileup_bam, s.pileup_bam + ".bai")
    if not cfg["keep.qc_fastq"]:
        rm(s.sample_fastq, s.filtered_fastq)

    LOG.info("[%s] %s  (in-frame %.0f%%, CDS %.0f%% of genic, %.0fs)", acc, q["verdict"],
             100 * q["periodicity_inframe_frac"], 100 * q["cds_frac_of_genic"],
             time.time() - t0)
    return q


def _total_reads(acc: str, meta: dict) -> int | None:
    """Reads in the whole run -- used to project how much usable data it holds."""
    try:
        return int((meta.get(acc) or {})["read_count"])
    except (KeyError, TypeError, ValueError):
        pass
    from .sra import metadata
    try:
        return metadata.read_count(acc)
    except Exception:  # noqa: BLE001  -- a missing projection must not fail the sample
        return None


# --- stage architecture: the call, from the profile JSON alone --------------
def arch_sample(acc: str, cfg: Config) -> dict:
    s = Sample(acc, cfg.workdir)
    prof = read_json(s.profile_json)
    if prof is None:
        raise RuntimeError(f"{acc}: no profile ({s.profile_json}); run the qc stage first")

    call = architecture.call(prof, cfg)
    write_json(s.arch_json, call)
    if cfg["plots.enabled"]:
        architecture.plot_call(prof, call, s.arch_plot(cfg["plots.format"]),
                               dpi=cfg["plots.dpi"])
    LOG.info("[%s] %s", acc, call["structure"][:120])
    return call


# --- stage bam: download -> trim -> filter -> map (-> dedup) ----------------
def is_processed(cfg: Config, s: Sample) -> bool:
    """Has the bam stage finished for this sample? With `keep.bam` off the BAM is
    deleted on purpose, so the process JSON is the evidence, not the BAM."""
    return nonempty(s.process_json) and (nonempty(s.bam) or not cfg["keep.bam"])


def process_sample(acc: str, cfg: Config, src: str) -> dict:
    s = Sample(acc, cfg.workdir)
    resume = cfg["pipeline.resume"]
    threads = cfg["project.threads"]

    if resume and is_processed(cfg, s):
        LOG.info("[%s] bam: cached", acc)
        return read_json(s.process_json)

    call = read_json(s.arch_json)
    if call is None:
        raise RuntimeError(f"{acc}: no architecture call; run the architecture stage first")
    if call["status"] != "ok" and not cfg["architecture.process_undetermined"]:
        raise RuntimeError(
            f"{acc}: architecture is '{call['status']}' ({call.get('reason', '')[:60]}); "
            f"set architecture.process_undetermined=true to trim it with a best-effort plan")

    # each step's own stats, verbatim, in its own block (`reports.process_tsv` reads them)
    info: dict = {"run_accession": acc}

    # 1. the full dataset; a local FASTQ needs no download
    if os.path.exists(src):
        full = src
        info["download"] = {"route": "local", "bytes": os.path.getsize(src)}
    else:
        dl = download.download_full(acc, s.full_fastq, cfg, log=s.log)
        full = s.full_fastq
        info["download"] = dl
        LOG.info("[%s] downloaded %s via %s at %.1f MB/s", acc,
                 human(dl.get("bytes")), dl["route"], dl.get("mb_per_s") or 0)

    # 2. trim to the footprint (the RT base stays); UMIs go to the read name
    if not (resume and nonempty(s.trimmed_fastq) and nonempty(s.trim_json)):
        tst = architecture.trim_fastq(call, full, s.trimmed_fastq, cfg)
        write_json(s.trim_json, tst)
        LOG.info("[%s] trimmed: %s of %s reads kept", acc,
                 human(tst["n_reads_out"]), human(tst["n_reads_in"]))
    info["trim"] = read_json(s.trim_json, {})

    # 3. contaminants again, on the whole run
    to_map = s.trimmed_fastq
    contam_full = os.path.join(s.dir, f"{acc}.contam_full.json")
    if cfg["process.filter_contaminants"]:
        if not (resume and nonempty(s.clean_fastq) and nonempty(contam_full)):
            cst = contaminants.filter_fastq(s.trimmed_fastq, s.clean_fastq, cfg,
                                            label=acc, threads=threads, log=s.log)
            write_json(contam_full, cst)
        info["contaminants"] = read_json(contam_full, {})
        to_map = s.clean_fastq

    # 4. the deliverable alignment; STAR counts reads into genes as it aligns
    bam = star.align_final(to_map, s.star_final, cfg, threads=threads, log=s.log)
    info["mapping"] = star.parse_log(os.path.join(s.star_final, "Log.final.out"))

    gene_counts = counts.path(s.star_final)
    if nonempty(gene_counts):
        info["counts"] = counts.read_counts(gene_counts, label=acc)[1]

    # 5. pile-up removal, on the unsorted BAM STAR just wrote
    if cfg["process.filter_pileups"]:
        pb = os.path.join(s.star_final, "aligned.pileup.bam")
        pst = pileups.filter_bam(bam, pb, cfg, label=acc)
        write_json(os.path.join(s.dir, f"{acc}.pileup_full.json"), pst)
        info["pileup"] = pst
        rm(bam)
        bam = pb

    # 6. sort + index
    star.sort_index(bam, s.bam, threads=min(threads, 8))
    rm(bam)

    # 7. optional UMI deduplication (off by default: without a UMI it would
    #    collapse genuine footprints that share a start position)
    info["umi_dedup"] = bool(cfg["process.umi_dedup"])
    if cfg["process.umi_dedup"]:
        tmp = s.bam + ".dedup.bam"
        info["dedup"] = proc_dedup.dedup(s.bam, tmp, cfg, log=s.log)
        os.replace(tmp, s.bam)
        rm(s.bam + ".bai", tmp + ".bai")
        star.ensure_sorted_indexed(s.bam, threads=min(threads, 8))
        LOG.info("[%s] UMI dedup: %s -> %s reads (%.0f%% duplicates)", acc,
                 human(info["dedup"]["n_in"]), human(info["dedup"]["n_out"]),
                 100 * (1 - info["dedup"]["frac_kept"]))
        if "counts" in info:
            LOG.warning("[%s] the gene counts are NOT deduplicated: STAR counts while it "
                        "aligns, and UMI dedup happens after.", acc)

    # 8. periodicity of the finished BAM. The QC verdict was measured on a sample of
    #    the untrimmed reads; a mis-trimmed footprint boundary shows up only here.
    info["periodicity"] = verdict.periodicity(
        s.bam, cfg.annotation_index, max_reads=cfg["qc.max_reads_scored"], label=acc)

    info["bam"] = s.bam
    info["bam_bytes"] = os.path.getsize(s.bam)
    info["keep"] = {k: bool(cfg[f"keep.{k}"])
                    for k in ("bam", "fastq", "trimmed_fastq", "clean_fastq")}
    write_json(s.process_json, info)

    # 9. housekeeping, last: everything above was measured before its input goes
    if not cfg["keep.fastq"] and not os.path.exists(src):
        rm(s.full_fastq)                       # never the user's own input FASTQ
    if not cfg["keep.trimmed_fastq"]:
        rm(s.trimmed_fastq)
    if not cfg["keep.clean_fastq"]:
        rm(s.clean_fastq)
    if not cfg["keep.bam"]:
        rm(s.bam, s.bam + ".bai")              # STAR's gene count table stays

    LOG.info("[%s] BAM %s (%.0f%% uniquely mapped)", acc, human(info["bam_bytes"]),
             100 * (info.get("mapping", {}).get("frac_unique") or 0))
    return info


# --- the driver -------------------------------------------------------------
def _pool(cfg: Config, fn, items, stage: str) -> tuple[list[str], list[dict]]:
    """Run `fn(acc, ...)` over `items` in a process pool. One sample's failure is
    recorded and the rest carry on."""
    jobs = max(1, int(cfg["project.jobs"]))
    ok: list[str] = []
    failed: list[dict] = []
    if not items:
        return ok, failed

    LOG.info("── stage %s: %d dataset(s), %d job(s) x %d thread(s)",
             stage, len(items), jobs, cfg["project.threads"])
    t0 = time.time()
    if jobs == 1:
        for args in items:
            acc = args[0]
            try:
                fn(*args)
                ok.append(acc)
            except Exception as exc:  # noqa: BLE001
                LOG.error("[%s] %s failed: %s", acc, stage, exc)
                failed.append({"run_accession": acc, "stage": stage, "error": str(exc),
                               "traceback": traceback.format_exc(limit=3)})
    else:
        # A worker that is KILLED (OOM, a segfault in STAR) breaks the whole pool and
        # every pending future fails with BrokenProcessPool. Retry what did not run:
        #   * at half the width -- a memory wave fits in a narrower pool;
        #   * and if a single worker still dies, the head of the queue is the dataset
        #     that kills it: fail that one alone and return to full width.
        remaining = list(items)
        workers = jobs
        while remaining:
            done: set[str] = set()
            broke = False
            with ProcessPoolExecutor(max_workers=workers) as ex:
                futs = {ex.submit(fn, *args): args[0] for args in remaining}
                for fut in as_completed(futs):
                    acc = futs[fut]
                    try:
                        fut.result()
                        ok.append(acc)
                        done.add(acc)
                    except BrokenProcessPool:
                        broke = True          # never ran, or was the one killed: retry
                    except Exception as exc:  # noqa: BLE001
                        LOG.error("[%s] %s failed: %s", acc, stage, exc)
                        failed.append({"run_accession": acc, "stage": stage, "error": str(exc)})
                        done.add(acc)
            remaining = [args for args in remaining if args[0] not in done]
            if not broke or not remaining:
                break
            if workers > 1:
                workers = max(1, workers // 2)
                LOG.error("stage %s: a worker was killed (pool broken); retrying %d "
                          "remaining dataset(s) at %d job(s)", stage, len(remaining), workers)
                continue
            poison = remaining[0]
            LOG.error("stage %s: [%s] repeatedly killed its worker process (pool broken); "
                      "recording it failed and continuing with %d other dataset(s)",
                      stage, poison[0], len(remaining) - 1)
            failed.append({"run_accession": poison[0], "stage": stage,
                           "error": "killed its worker process (pool broken); a crash or "
                                    "OOM on this dataset's reads"})
            remaining = remaining[1:]
            workers = jobs

    LOG.info("── stage %s done: %d ok, %d failed (%.0fs)",
             stage, len(ok), len(failed), time.time() - t0)
    return sorted(ok), failed


def run(cfg: Config) -> dict:
    cfg.validate()
    stages = stages_to_run(cfg["pipeline.start"], cfg["pipeline.end"])
    LOG.info("RiboMine: %s -> %s   (stages: %s)",
             cfg["pipeline.start"], cfg["pipeline.end"], " -> ".join(stages))
    LOG.info("workdir: %s", cfg.workdir)

    accs, src = resolve_inputs(cfg)
    if not accs:
        LOG.warning("no datasets to process")
        return {"n_input": 0}

    all_failed: list[dict] = []
    summary: dict = {"n_input": len(accs), "workdir": cfg.workdir}

    # Everything the workers share is built here, before the pool forks: two
    # processes building the same file at once corrupt it.
    annotation.ensure_index(cfg)
    contaminants.ensure_index(cfg)
    _ensure_fasta_index(cfg)
    meta = _run_meta(cfg)

    # one copy of the STAR index in shared memory for the whole batch
    star.load_genome(cfg)
    LOG.info("STAR genome load: %s", star.genome_load())
    try:
        ok, failed = _pool(cfg, qc_sample, [(a, cfg, src[a], meta) for a in accs], "qc")
        all_failed += failed
        reports.qc_tsv(cfg, accs)
        summary["qc"] = _tally(cfg, ok, "qc_json", "verdict")
        LOG.info("QC verdicts: %s", ", ".join(f"{k} {v}" for k, v in summary["qc"].items()))

        if cfg["pipeline.end"] == "qc":
            return _finish(cfg, accs, summary, all_failed)

        # architecture, on the runs QC kept
        keep = set(cfg["pipeline.keep_verdicts"])
        passed = [a for a in ok if (read_json(Sample(a, cfg.workdir).qc_json, {})
                                    .get("verdict") in keep)]
        LOG.info("carrying %d/%d dataset(s) past QC (verdicts kept: %s)",
                 len(passed), len(ok), ", ".join(sorted(keep)))

        ok_a, failed = _pool(cfg, arch_sample, [(a, cfg) for a in passed], "architecture")
        all_failed += failed
        reports.arch_tsv(cfg, passed)
        summary["architecture"] = _tally(cfg, ok_a, "arch_json", "status")

        if cfg["pipeline.end"] == "architecture":
            return _finish(cfg, accs, summary, all_failed)

        # the full datasets
        if not cfg["architecture.process_undetermined"]:
            ok_a = [a for a in ok_a
                    if read_json(Sample(a, cfg.workdir).arch_json, {}).get("status") == "ok"]
        ok_b, failed = _pool(cfg, process_sample, [(a, cfg, src[a]) for a in ok_a], "bam")
        all_failed += failed
        reports.process_tsv(cfg, ok_a)
        reports.counts_tsv(cfg, ok_b)
        summary["bam"] = len(ok_b)
    finally:
        star.unload_genome(cfg)

    return _finish(cfg, accs, summary, all_failed)


def _finish(cfg: Config, accs: list[str], summary: dict, failed: list[dict]) -> dict:
    if failed:
        p = write_tsv(os.path.join(cfg.workdir, "failed.tsv"), failed,
                      ["run_accession", "stage", "error"])
        LOG.warning("%d dataset(s) failed -- see %s", len(failed), p)
    summary["n_failed"] = len(failed)
    print("\n" + reports.summary_line(cfg, accs))
    return summary


def _tally(cfg: Config, accs: list[str], json_attr: str, key: str) -> dict:
    """How many samples carry each value of `key` in one of their stage JSONs."""
    out: dict[str, int] = {}
    for a in accs:
        v = read_json(getattr(Sample(a, cfg.workdir), json_attr), {}).get(key, "?")
        out[v] = out.get(v, 0) + 1
    return out
