"""Stage orchestration: what runs, in what order, and what a resume skips.

The pipeline has four stages and three entry points into them:

    query ──► qc ──► architecture ──► bam
      ▲        ▲
      │        └── start="accessions" (a list) or start="fastq" (a local directory)
      └── start="query"

and it stops after the stage named by `pipeline.end` ("qc", "architecture" or
"bam"). The two early end points are real deliverables, not debug hooks: each
writes a TSV and a per-dataset figure.

Why the stages split where they do: the QC stage already samples, filters, maps
(locally) and profiles the reads, and the architecture stage reads *only* the
profile JSON that produced. So calling the architecture costs milliseconds once
QC has run, and the architecture thresholds can be re-tuned over a whole cohort
without touching a BAM.

Per-sample work runs in a process pool. A sample that fails is recorded in
failed.tsv and does not take the run down -- at 500 datasets, something always
fails.
"""
from __future__ import annotations

import os
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool

from . import reports
from .arch import infer as arch_infer
from .arch import plot as arch_plot
from .arch import trim as arch_trim
from .config import Config, stages_to_run
from .process import counts
from .process import dedup as proc_dedup
from .process import star
from .qc import annotation, contaminants, pileups, plot as qc_plot, profile, verdict
from .sra import download, query
from .utils import (LOG, Sample, human, nonempty, read_json, read_lines, read_tsv,
                    rm, write_json, write_tsv)


# ---------------------------------------------------------------------------
# sample list: the three start points
# ---------------------------------------------------------------------------
def resolve_inputs(cfg: Config) -> tuple[list[str], dict[str, str]]:
    """(accessions, source_map). `source_map[acc]` is what the reads come from:
    an accession (streamed from ENA) or a local FASTQ path."""
    start = cfg["pipeline.start"]

    if start == "fastq":
        d = cfg._abs(cfg["pipeline.fastq_dir"])
        exts = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
        files = sorted(f for f in os.listdir(d) if f.endswith(exts))
        if not files:
            raise RuntimeError(f"no FASTQ files in {d}")
        src = {}
        for f in files:
            label = f
            for e in exts:            # strip the longest matching extension
                if label.endswith(e):
                    label = label[: -len(e)]
                    break
            src[label] = os.path.join(d, f)
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
    rows = read_tsv(tsv)
    accs = [r["run_accession"] for r in rows]
    LOG.info("start=query: %d candidate run(s) from %s", len(accs), tsv)
    return accs, {a: a for a in accs}


def _ensure_fasta_index(cfg: Config) -> None:
    """Build the genome `.fai` up front.

    `pysam.FastaFile` builds a missing `.fai` *in place*, silently. In a process
    pool that means every worker starts building the same index for the same 3 GB
    FASTA at the same time, and their interleaved writes leave a truncated one --
    after which `fetch()` returns the wrong sequence and every architecture call
    downstream is quietly fabricated. Build it once, here, where only one process
    is running.
    """
    fasta = cfg.ref("genome_fasta")
    fai = fasta + ".fai"
    if nonempty(fai):
        return
    import pysam

    LOG.info("indexing the genome FASTA (once): %s", fasta)
    t0 = time.time()
    pysam.faidx(fasta)
    LOG.info("  -> %s (%.0fs)", fai, time.time() - t0)


def _run_meta(cfg: Config) -> dict[str, dict]:
    """run_accession -> its candidates.tsv row (read_count, study, ...), if any."""
    tsv = os.path.join(cfg.workdir, "meta", "candidates.tsv")
    if not nonempty(tsv):
        return {}
    return {r["run_accession"]: r for r in read_tsv(tsv)}


# ---------------------------------------------------------------------------
# stage 2: read sample -> contaminant filter -> local map -> profile -> QC
# ---------------------------------------------------------------------------
def qc_sample(acc: str, cfg: Config, src: str, meta: dict) -> dict:
    s = Sample(acc, cfg.workdir)
    resume = cfg["pipeline.resume"]
    threads = cfg["project.threads"]

    if resume and nonempty(s.qc_json) and nonempty(s.profile_json):
        LOG.info("[%s] qc: cached", acc)
        return read_json(s.qc_json)

    t0 = time.time()

    # 1. a uniform read sample -- nothing is downloaded in full at this stage
    if not (resume and nonempty(s.sample_fastq)):
        st = download.sample_reads(src, s.sample_fastq, n=cfg["qc.sample_reads"],
                                   scan=cfg["qc.scan_reads"], seed=cfg["project.seed"])
        LOG.info("[%s] sampled %s reads of %s scanned", acc,
                 human(st["n_sampled"]), human(st["n_scanned"]))
        if st.get("sorted_warning"):
            LOG.warning("[%s] %s", acc, st["sorted_warning"])

    # 2. contaminants: rRNA/tRNA/snRNA/Mt (bowtie2) + low-complexity. Reads stay
    #    UNTRIMMED -- the architecture is read out of the soft clips, so the 5'
    #    construct and the 3' adapter must still be on the read.
    if not (resume and nonempty(s.filtered_fastq) and nonempty(s.contam_json)):
        cst = contaminants.filter_fastq(s.sample_fastq, s.filtered_fastq, cfg,
                                        label=acc, threads=threads, log=s.log)
        write_json(s.contam_json, cst)
    cst = read_json(s.contam_json, {})

    # 3. LOCAL alignment: everything non-genomic lands in the soft clips
    if not (resume and nonempty(s.local_bam)):
        star.align_local(s.filtered_fastq, s.star_local, cfg, threads=threads, log=s.log)

    # 4. position pile-ups (adapter dimers, fixed contaminants) -- data-driven
    if not (resume and nonempty(s.pileup_bam) and nonempty(s.pileup_json)):
        pst = pileups.filter_bam(s.local_bam, s.pileup_bam, cfg, label=acc)
        write_json(s.pileup_json, pst)
    pst = read_json(s.pileup_json, {})

    # 5. positional profile (feeds the architecture caller). The anchor gate is the
    #    architecture's own, so the two cannot drift apart.
    if not (resume and nonempty(s.profile_json)):
        prof = profile.profile_bam(s.pileup_bam, cfg.ref("genome_fasta"), label=acc,
                                   max_reads=cfg["qc.sample_reads"],
                                   min_panel_frac=cfg["architecture.min_anchor_frac"])
        write_json(s.profile_json, prof)

    # 6. the verdict
    total = _total_reads(acc, meta)
    q = verdict.qc(s.pileup_bam, cfg.annotation_index, cfg, star_log=s.local_log,
                   label=acc, contam=cst, pileup=pst, total_reads=total)
    write_json(s.qc_json, q)

    if cfg["plots.enabled"]:
        qc_plot.plot_qc(q, s.qc_plot(cfg["plots.format"]), dpi=cfg["plots.dpi"])

    # 7. the working data behind the verdict and the architecture call. Every number
    #    they produced is already in qc.json / profile.json (and the architecture
    #    stage reads ONLY the profile), so by default the reads and their alignments
    #    go -- at several hundred runs they are the bulk of the workdir. Kept, they
    #    are sorted and indexed: the BAM a call was computed on is the one to open in
    #    a browser when the call looks wrong. STAR writes them unsorted, and the steps
    #    above read them in file order, so that has to happen here, at the end.
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
    row = meta.get(acc) or {}
    try:
        return int(row["read_count"])
    except (KeyError, TypeError, ValueError):
        pass
    from .sra import metadata
    try:
        return metadata.read_count(acc)
    except Exception:  # noqa: BLE001  -- a missing projection must not fail the sample
        return None


# ---------------------------------------------------------------------------
# stage 3: the architecture call (reads only the profile JSON)
# ---------------------------------------------------------------------------
def arch_sample(acc: str, cfg: Config) -> dict:
    s = Sample(acc, cfg.workdir)
    prof = read_json(s.profile_json)
    if prof is None:
        raise RuntimeError(f"{acc}: no profile ({s.profile_json}); run the qc stage first")

    call = arch_infer.infer(prof, arch_infer.Thresholds.from_config(cfg))
    write_json(s.arch_json, call)
    if cfg["plots.enabled"]:
        arch_plot.plot_arch(prof, call, s.arch_plot(cfg["plots.format"]), dpi=cfg["plots.dpi"])

    if call["status"] == "ok":
        LOG.info("[%s] %s", acc, arch_infer.architecture_string(call))
    else:
        LOG.info("[%s] architecture %s: %s", acc, call["status"], call.get("reason", "")[:80])
    return call


# ---------------------------------------------------------------------------
# stage 4: full download -> trim -> filter -> map (-> dedup)
# ---------------------------------------------------------------------------
def is_processed(cfg: Config, s: Sample) -> bool:
    """Has stage 4 already finished for this sample? (What `resume` skips.)

    What "finished" looks like on disk depends on what the config KEEPS. With
    `keep.bam` off the deliverable is deleted on purpose, so its absence is not
    evidence that the sample needs re-processing -- the process JSON is. Testing for
    the BAM would silently re-download and re-map an entire cohort on every resume of
    a counts-only run.
    """
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

    # The process record keeps each step's own stats verbatim, in its own block, so
    # the record stays a faithful account of what ran rather than a lossy summary --
    # and `reports.process_tsv` reads those blocks by name.
    info: dict = {"run_accession": acc}

    # 1. the full dataset. A local FASTQ start point skips the download entirely.
    if os.path.exists(src):
        full = src
        info["download"] = {"route": "local", "bytes": os.path.getsize(src)}
    else:
        dl = download.download_full(acc, s.full_fastq, cfg, log=s.log)
        full = s.full_fastq
        info["download"] = dl
        LOG.info("[%s] downloaded %s via %s at %.1f MB/s", acc,
                 human(dl.get("bytes")), dl["route"], dl.get("mb_per_s") or 0)

    # 2. trim with the inferred architecture: 5' construct off (RT base kept),
    #    3' construct + adapter off, UMI content into the read name.
    if not (resume and nonempty(s.trimmed_fastq) and nonempty(s.trim_json)):
        tst = arch_trim.trim_fastq(full, call, s.trimmed_fastq,
                                   min_len=cfg["process.min_len"], label=acc,
                                   discard_untrimmed=cfg["process.discard_untrimmed"],
                                   min_overlap=cfg["process.adapter_min_overlap"])
        write_json(s.trim_json, tst)
    info["trim"] = read_json(s.trim_json, {})

    # 3. contaminants again -- on the whole run this time, not on the 200k sample
    to_map = s.trimmed_fastq
    contam_full = os.path.join(s.dir, f"{acc}.contam_full.json")
    if cfg["process.filter_contaminants"]:
        if not (resume and nonempty(s.clean_fastq) and nonempty(contam_full)):
            cst = contaminants.filter_fastq(s.trimmed_fastq, s.clean_fastq, cfg,
                                            label=acc, threads=threads, log=s.log)
            write_json(contam_full, cst)
        info["contaminants"] = read_json(contam_full, {})
        to_map = s.clean_fastq

    # 4. map end-to-end -- the deliverable alignment. STAR counts the reads into genes
    #    as it aligns them (--quantMode GeneCounts), so the count matrix costs nothing
    #    extra; `reports.counts_tsv` joins the per-run tables at the end of the batch.
    bam = star.align_final(to_map, s.star_final, cfg, threads=threads, log=s.log)
    info["mapping"] = star.parse_log(os.path.join(s.star_final, "Log.final.out"))

    gene_counts = counts.path(s.star_final)
    if nonempty(gene_counts):
        info["counts"] = counts.read_counts(gene_counts, label=acc)[1]
        LOG.info("[%s] gene counts: %s reads in %s genes (%.0f%% of counted reads)", acc,
                 human(info["counts"]["n_in_genes"]),
                 human(info["counts"]["n_genes_detected"]),
                 100 * (info["counts"].get("frac_in_genes") or 0))

    # 5. pile-up removal (adapter dimers / fixed contaminants that survived the
    #    sequence filter). Cheapest on the unsorted BAM STAR just wrote.
    if cfg["process.filter_pileups"]:
        pb = os.path.join(s.star_final, "aligned.pileup.bam")
        pst = pileups.filter_bam(bam, pb, cfg, label=acc)
        write_json(os.path.join(s.dir, f"{acc}.pileup_full.json"), pst)
        info["pileup"] = pst
        rm(bam)
        bam = pb

    # 6. sort + index -- the BAM users get
    star.sort_index(bam, s.bam, threads=min(threads, 8))
    rm(bam)

    # 7. optional UMI deduplication. OFF by default, and deliberately so: it is only
    #    correct when the library actually carries a UMI, and on a library that does
    #    not it silently collapses genuine duplicate footprints -- which in ribo-seq
    #    are real signal (a highly translated codon IS covered many times).
    info["umi_dedup"] = bool(cfg["process.umi_dedup"])
    if cfg["process.umi_dedup"]:
        tmp = s.bam + ".dedup.bam"
        info["dedup"] = proc_dedup.dedup(s.bam, tmp, cfg, log=s.log)
        os.replace(tmp, s.bam)
        rm(s.bam + ".bai", tmp + ".bai")
        # Both backends hand back the coordinate order they were given, so this is an
        # index and not a re-sort -- but it is `ensure_sorted_indexed` rather than
        # `index` so that the "every BAM RiboMine leaves behind is sorted and indexed"
        # invariant holds because it is CHECKED, not because a third-party tool is
        # assumed to have been well behaved.
        star.ensure_sorted_indexed(s.bam, threads=min(threads, 8))
        LOG.info("[%s] UMI dedup: %s -> %s reads (%.0f%% duplicates)", acc,
                 human(info["dedup"]["n_in"]), human(info["dedup"]["n_out"]),
                 100 * (1 - info["dedup"]["frac_kept"]))
        if "counts" in info:
            # STAR counted during the alignment, which is before this step ran. Say so
            # once, per sample, rather than let someone discover it in a volcano plot.
            LOG.warning("[%s] the gene counts are NOT deduplicated: STAR counts while it "
                        "aligns, and UMI dedup happens after. The BAM is deduplicated; "
                        "the count matrix is of the reads that went into it.", acc)

    # 8. periodicity of the reads we are actually handing over. The QC verdict was
    #    decided on a 200k-read sample of the UNTRIMMED reads, locally aligned; this
    #    is the finished article -- trimmed, filtered, deduplicated, end-to-end. A
    #    trim that cut the footprint boundary wrong shows up here and nowhere else.
    info["periodicity"] = verdict.periodicity(
        s.bam, cfg.annotation_index, max_reads=cfg["qc.max_reads_scored"], label=acc)

    info["bam"] = s.bam
    info["bam_bytes"] = os.path.getsize(s.bam)
    info["keep"] = {k: bool(cfg[f"keep.{k}"])
                    for k in ("bam", "fastq", "trimmed_fastq", "clean_fastq")}
    write_json(s.process_json, info)

    # 9. housekeeping. A ribo-seq run is 1-10 GB and a mining run holds hundreds of
    #    them, so nothing is kept unless `keep` says so. This is the LAST thing the
    #    stage does: every number above was measured before its input was deleted.
    if not cfg["keep.fastq"] and not os.path.exists(src):
        rm(s.full_fastq)                       # never the user's own input FASTQ
    if not cfg["keep.trimmed_fastq"]:
        rm(s.trimmed_fastq)
    if not cfg["keep.clean_fastq"]:
        rm(s.clean_fastq)
    if not cfg["keep.bam"]:
        # the counts and the periodicity are already measured and written; the gene
        # count table STAR wrote stays, so the matrix can still be built
        rm(s.bam, s.bam + ".bai")

    LOG.info("[%s] BAM %s (%.0f%% uniquely mapped)", acc, human(info["bam_bytes"]),
             100 * (info.get("mapping", {}).get("frac_unique") or 0))
    return info


# ---------------------------------------------------------------------------
# the driver
# ---------------------------------------------------------------------------
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
        # A worker that raises is recorded and the pool carries on -- but a worker that
        # is *killed* (the OOM killer picking one process under the startup memory spike,
        # a segfault in STAR/samtools) is different: it puts the whole ProcessPoolExecutor
        # into a broken state, and every future still pending then fails at once with
        # BrokenProcessPool. Left unhandled that discards the entire rest of the shard --
        # one dead worker throwing away ~1000 samples. So we treat a broken pool as a
        # transient event: the futures that had not run yet are retried in a fresh pool at
        # half the width (a narrower wave is less likely to trip the same OOM), down to 1,
        # for a bounded number of rounds. Only then is the remainder recorded as failed.
        remaining = list(items)
        workers = jobs
        rounds_at_one = 0
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
                        # This future never ran (or was the one killed); leave it out of
                        # `done` so it is retried below. Do not record it as a failure.
                        broke = True
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
            # Already down to a single worker and it still died: this is not the memory
            # wave, it is one dataset that kills its process. Give it one more lone pass to
            # get past it, then stop looping and record whatever is left as failed rather
            # than churn forever.
            rounds_at_one += 1
            if rounds_at_one >= 2:
                LOG.error("stage %s: %d dataset(s) still break a single-worker pool; "
                          "recording them as failed", stage, len(remaining))
                for args in remaining:
                    failed.append({"run_accession": args[0], "stage": stage,
                                   "error": "worker killed repeatedly (pool broken); "
                                            "likely a crash or OOM on this dataset"})
                remaining = []
            else:
                LOG.error("stage %s: a single-worker pool broke; one more lone pass over "
                          "%d remaining dataset(s)", stage, len(remaining))

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

    # ---- everything the workers SHARE gets built here, in the parent, before the
    # pool forks. Anything built lazily inside a worker is built by all of them at
    # once, and two processes writing the same file is corruption, not a slowdown.
    annotation.ensure_index(cfg)     # the 22 MB GTF pickle
    contaminants.ensure_index(cfg)   # the bowtie2 contaminant index (bundled FASTA by default)
    _ensure_fasta_index(cfg)         # the genome .fai -- pysam builds it silently otherwise
    meta = _run_meta(cfg)

    # STAR's ~28 GB genome index goes into shared memory once for the whole batch
    # rather than once per worker. At jobs=8 the difference is 28 GB against 224 GB;
    # this is what makes a several-hundred-dataset run possible at all.
    star.load_genome(cfg)
    LOG.info("STAR genome load: %s", star.genome_load())
    try:
        ok, failed = _pool(cfg, qc_sample, [(a, cfg, src[a], meta) for a in accs], "qc")
        all_failed += failed
        reports.qc_tsv(cfg, accs)
        summary["qc"] = _verdict_counts(cfg, ok)
        LOG.info("QC verdicts: %s", ", ".join(f"{k} {v}" for k, v in summary["qc"].items()))

        if cfg["pipeline.end"] == "qc":
            return _finish(cfg, accs, summary, all_failed)

        # ---- stage 3: architecture, on the runs QC kept
        keep = set(cfg["pipeline.keep_verdicts"])
        passed = [a for a in ok if (read_json(Sample(a, cfg.workdir).qc_json, {})
                                    .get("verdict") in keep)]
        LOG.info("carrying %d/%d dataset(s) past QC (verdicts kept: %s)",
                 len(passed), len(ok), ", ".join(sorted(keep)))

        ok_a, failed = _pool(cfg, arch_sample, [(a, cfg) for a in passed], "architecture")
        all_failed += failed
        reports.arch_tsv(cfg, passed)
        summary["architecture"] = _arch_counts(cfg, ok_a)

        if cfg["pipeline.end"] == "architecture":
            return _finish(cfg, accs, summary, all_failed)

        # ---- stage 4: the full datasets
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


def _verdict_counts(cfg: Config, accs: list[str]) -> dict:
    out: dict[str, int] = {}
    for a in accs:
        v = read_json(Sample(a, cfg.workdir).qc_json, {}).get("verdict", "?")
        out[v] = out.get(v, 0) + 1
    return out


def _arch_counts(cfg: Config, accs: list[str]) -> dict:
    out: dict[str, int] = {}
    for a in accs:
        st = read_json(Sample(a, cfg.workdir).arch_json, {}).get("status", "?")
        out[st] = out.get(st, 0) + 1
    return out
