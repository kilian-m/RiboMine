"""STAR alignment: the permissive LOCAL pass, and the end-to-end pass.

Two alignments, for two different jobs:

* `align_local` -- the QC / architecture alignment. The reads are **untrimmed**:
  the adapter is an *output* of this pipeline, not an input, so nothing may be
  removed before we look. STAR runs with `--alignEndsType Local` and permissive
  length/score filters so that everything non-genomic (5' UMI, RT additions,
  3' UMI, sample barcode, adapter) is pushed into the **soft clips** instead of
  preventing the alignment. The architecture is then read back out of those
  clips. These parameters are what make the architecture readable at all, so
  they are hard-coded here rather than exposed in the config.

* `align_final` -- the deliverable BAM. The trimmed reads, aligned end-to-end
  with the `mapping` section of the config. Mapping metrics are taken from this
  alignment, never from the local one (a permissive local alignment inflates
  multimapping).

Both write an unsorted BAM; `sort_index` produces the coordinate-sorted+indexed
file. That split is not cosmetic: STAR's own BAM sorting is incompatible with
the shared-memory genome (`--genomeLoad LoadAndKeep`), and the shared genome is
what makes a 500-sample run affordable -- the 28 GB index is loaded once for the
whole batch instead of once per sample. `pipeline.run` calls `load_genome` before
the process pool and `unload_genome` in a `finally`; see below.
"""
from __future__ import annotations

import os
import re
import shutil

from ..config import Config
from ..utils import LOG, nonempty, require_tools, rm, run

# --- the shared-memory genome: how the mode gets from the parent to the workers
#
# The per-sample aligners run in a `project.jobs`-wide PROCESS POOL. Without a
# shared genome each worker mmaps its own copy of the ~28 GB index -- at jobs=8
# that is 224 GB and the box dies. So the batch loads the index into shared
# memory ONCE (`--genomeLoad LoadAndExit`), every sample then attaches to that
# one segment (`--genomeLoad LoadAndKeep`), and the batch drops it at the end
# (`--genomeLoad Remove`). This is what the reference batch driver
# (read_architecture/bin/run_batch.sh) does, and it is what makes a
# several-hundred-dataset run affordable.
#
# The "are we sharing?" decision is made in the parent but has to be visible to
# the workers, so it is kept in the ENVIRONMENT rather than in a module global:
# a global mutated in the parent after the pool has forked would not reach the
# children (and would not survive a `spawn` start method at all), whereas
# `os.environ` is inherited by both fork and spawn children. `load_genome`
# writes the effective mode into RIBOMINE_STAR_GENOME_LOAD before the pool is
# created; `align_local`/`align_final` read it back through `genome_load()`.
#
# Resolution order, most specific first:
#   RIBOMINE_STAR_GENOME_LOAD -- the decision this batch actually made (it already
#                                accounts for the user's wish AND for whether the
#                                shared load succeeded)
#   STAR_GENOME_LOAD          -- the user's override, honoured when no batch
#                                decision has been made (e.g. a single align_local
#                                call outside `pipeline.run`)
#   NoSharedMemory            -- safe default: every process loads its own index
STATE_ENV = "RIBOMINE_STAR_GENOME_LOAD"
GENOME_LOAD_ENV = "STAR_GENOME_LOAD"
DEFAULT_GENOME_LOAD = "NoSharedMemory"
SHARED = "LoadAndKeep"
# the modes an *alignment* may legally run with. LoadAndExit/Remove are genome
# management ops, not alignment modes; LoadAndRemove has the first sample to
# finish unload the genome out from under the rest of the pool, so the pipeline
# never selects it -- but a user who asks for it explicitly gets it.
ALIGN_LOADS = ("NoSharedMemory", "LoadAndKeep", "LoadAndRemove")

# The permissive LOCAL alignment (identical to the reference map_reads.sh).
#   MD                            -- cheap, and consumers (profiler, pile-ups) use it
#   outSAMunmapped None           -- the profile only ever looks at aligned reads
#   alignEndsType Local           -- the whole point: non-genomic ends -> soft clips
#   outFilterMultimapNmax 1       -- unique only; a multimapper has no single genomic
#                                    context to compare the read bases against
#   outFilterMatchNmin 20         -- 20 genomic bases is enough to place a footprint...
#   outFilterMatchNminOverLread 0 -- ...and the *fraction* filters must be off, or a
#   outFilterScoreMinOverLread 0     read that is half construct is thrown away for
#                                    being "too short" -- exactly the reads we need
#   outFilterMismatchNmax 3 / NoverLmax 0.12 -- tolerate real mismatches inside the
#                                    footprint without letting construct bases in
#   seedSearchStartLmax 20        -- seed inside a short footprint that is flanked by
#                                    non-genomic sequence
#   outSJtype None                -- no splice-junction output; nothing reads it
LOCAL_ARGS: list[str] = [
    "--outSAMtype", "BAM", "Unsorted",
    "--outSAMattributes", "NH", "HI", "AS", "nM", "MD",
    "--outSAMunmapped", "None",
    "--alignEndsType", "Local",
    "--outFilterMultimapNmax", "1",
    "--outFilterMatchNmin", "20",
    "--outFilterMatchNminOverLread", "0",
    "--outFilterScoreMinOverLread", "0",
    "--outFilterMismatchNmax", "3",
    "--outFilterMismatchNoverLmax", "0.12",
    "--seedSearchStartLmax", "20",
    "--alignSJoverhangMin", "8",
    "--alignSJDBoverhangMin", "2",
    "--outSJtype", "None",
]

_STAR_BAM = "Aligned.out.bam"     # what STAR writes
_OUT_BAM = "aligned.bam"          # what we hand on (Sample.local_bam / star_final)
_LOG_FINAL = "Log.final.out"


# --- alignment -------------------------------------------------------------
def align_local(fastq: str, outdir: str, cfg: Config, *, threads: int = 8,
                log: str = "") -> str:
    """Permissive LOCAL alignment of the untrimmed sampled reads. Returns the BAM.

    Nothing is trimmed beforehand -- the 5' UMI, RT nt, 3' UMI, barcode and
    adapter are supposed to land in the soft clips, which is where the profiler
    reads the architecture from.
    """
    return _align(fastq, outdir, cfg, args=LOCAL_ARGS, threads=threads, log=log,
                  what="local")


def align_final(fastq: str, outdir: str, cfg: Config, *, threads: int = 8,
                log: str = "") -> str:
    """End-to-end alignment of the TRIMMED reads: the deliverable BAM.

    Parameters come from `cfg['mapping']` -- this alignment is the user's, and
    what "a mapped read" means for their downstream analysis is theirs to set.
    """
    args = [
        "--outSAMtype", "BAM", "Unsorted",
        "--outSAMattributes", "NH", "HI", "AS", "nM", "MD",
        "--alignEndsType", str(cfg.get("mapping.align_ends_type", "EndToEnd")),
        "--outFilterMultimapNmax", str(int(cfg.get("mapping.multimap_nmax", 1))),
        "--outFilterMismatchNmax", str(int(cfg.get("mapping.mismatch_nmax", 3))),
        "--outFilterMismatchNoverLmax", str(float(cfg.get("mapping.mismatch_noverlmax", 0.1))),
        "--outFilterMatchNminOverLread", str(float(cfg.get("mapping.match_nmin_over_lread", 0.9))),
        "--outSJtype", "None",
    ]
    intron_max = int(cfg.get("mapping.align_intron_max", 0) or 0)
    if intron_max:      # 0 == STAR's own default (spliced); only pass it when set
        args += ["--alignIntronMax", str(intron_max)]
    args += [str(a) for a in (cfg.get("mapping.extra_args") or [])]
    return _align(fastq, outdir, cfg, args=args, threads=threads, log=log, what="final")


def _align(fastq: str, outdir: str, cfg: Config, *, args: list[str], threads: int,
           log: str, what: str) -> str:
    require_tools("STAR")
    if not nonempty(fastq):
        raise ValueError(f"no reads to align: {fastq}")
    index = _index(cfg)
    os.makedirs(outdir, exist_ok=True)
    # a killed STAR leaves _STARtmp behind and the next run refuses to start
    rm(os.path.join(outdir, "_STARtmp"))

    # genome_load() reads the mode the batch settled on out of the environment, so a
    # pool worker sees it across the fork. NOTE: the shared genome (LoadAndKeep) is
    # incompatible with STAR's own BAM sorting -- which is why both arg sets say
    # `--outSAMtype BAM Unsorted` and sort_index() sorts with samtools afterwards.
    cmd = [
        "STAR",
        "--runMode", "alignReads",
        "--runThreadN", str(int(threads)),
        "--genomeDir", index,
        "--genomeLoad", genome_load(),
        "--readFilesIn", fastq,
        *_read_files_command(fastq),
        "--outFileNamePrefix", os.path.join(outdir, ""),
        *args,
    ]
    LOG.info("STAR %s: %s -> %s", what, os.path.basename(fastq), outdir)
    run(cmd, log_to=log or None)

    raw = os.path.join(outdir, _STAR_BAM)
    if not nonempty(raw):
        raise RuntimeError(f"STAR produced no alignments: {raw} (see {os.path.join(outdir, _LOG_FINAL)})")
    bam = os.path.join(outdir, _OUT_BAM)
    os.replace(raw, bam)
    return bam


def _index(cfg: Config) -> str:
    index = cfg.ref("star_index")
    if not index or not os.path.isdir(index):
        raise ValueError(f"reference.star_index is not a directory: {index!r}")
    return index


def _read_files_command(fastq: str) -> list[str]:
    return ["--readFilesCommand", "zcat"] if fastq.endswith(".gz") else []


# --- the shared-memory genome ----------------------------------------------
def genome_load() -> str:
    """The effective `--genomeLoad` mode for an alignment (see the note above).

    Public so the pipeline can log which mode a run actually got -- "is the genome
    shared?" is the difference between 28 GB and jobs x 28 GB of resident memory.
    """
    mode = (os.environ.get(STATE_ENV) or os.environ.get(GENOME_LOAD_ENV)
            or DEFAULT_GENOME_LOAD)
    if mode not in ALIGN_LOADS:
        # a typo (or a management op) in the override would make every STAR call
        # die; fall back rather than take the run down
        LOG.warning("ignoring --genomeLoad %r (not one of %s); using %s",
                    mode, ", ".join(ALIGN_LOADS), DEFAULT_GENOME_LOAD)
        return DEFAULT_GENOME_LOAD
    return mode


def load_genome(cfg: Config) -> bool:
    """Load the STAR index into shared memory once, for the whole batch.

    Call before the process pool is created: on success the per-sample aligners
    attach to the one segment with `--genomeLoad LoadAndKeep` instead of each
    loading their own ~28 GB copy.

    Idempotent, and it never raises. A shared segment can be legitimately
    unavailable (SHMALL/SHMMAX too small, another user's stale segment, no STAR on
    PATH); that is a performance problem, not a correctness one, so we warn and
    fall back to NoSharedMemory -- the run still produces the same BAMs, it just
    pays for the index per worker.

    Returns True if the batch is sharing the genome.
    """
    # already decided (a second call, e.g. a nested driver) -- do not re-load
    if STATE_ENV in os.environ:
        return os.environ[STATE_ENV] == SHARED

    want = os.environ.get(GENOME_LOAD_ENV)
    if want and want != SHARED:
        # the user asked for something else on purpose (small SHMALL, shared box)
        LOG.info("STAR shared genome not used ($%s=%s)", GENOME_LOAD_ENV, want)
        os.environ[STATE_ENV] = want if want in ALIGN_LOADS else DEFAULT_GENOME_LOAD
        return False

    if _genome_op(cfg, "LoadAndExit", "_star_load."):
        os.environ[STATE_ENV] = SHARED
        LOG.info("STAR genome loaded into shared memory (one copy for the whole batch)")
        return True

    # fall back: correctness is untouched, memory is not
    os.environ[STATE_ENV] = DEFAULT_GENOME_LOAD
    LOG.warning("STAR shared genome unavailable -- falling back to %s: each of the "
                "%d concurrent job(s) will load its OWN copy of the index. Lower "
                "project.jobs if the machine cannot hold that.",
                DEFAULT_GENOME_LOAD, int(cfg.get("project.jobs", 1) or 1))
    return False


def unload_genome(cfg: Config) -> bool:
    """Drop the shared-memory genome. Call it in a `finally` -- a leaked segment
    holds 28 GB of RAM until someone runs `STAR --genomeLoad Remove` by hand.

    Idempotent (a second call is a no-op) and never raises: a run that finished
    must not fail in its cleanup.
    """
    shared = os.environ.get(STATE_ENV) == SHARED
    # forget the decision, so this is a no-op the second time and a later
    # load_genome() may start over
    os.environ.pop(STATE_ENV, None)
    if not shared:
        return False
    if _genome_op(cfg, "Remove", "_star_remove."):
        LOG.info("STAR shared-memory genome removed")
        return True
    LOG.warning("could not remove the STAR shared-memory genome -- it may still hold "
                "RAM. Free it by hand: STAR --genomeLoad Remove --genomeDir %s",
                cfg.get("reference.star_index"))
    return False


def _genome_op(cfg: Config, mode: str, prefix: str) -> bool:
    """`STAR --genomeLoad <mode>` on the index. True on success; never raises."""
    try:
        if not shutil.which("STAR"):
            LOG.warning("STAR is not on PATH -- cannot share the genome")
            return False
        logdir = cfg.dir("logs")
        out_prefix = os.path.join(logdir, prefix)
        rm(out_prefix + "_STARtmp")        # a killed STAR leaves this and blocks the next
        run(["STAR", "--genomeLoad", mode, "--genomeDir", _index(cfg),
             "--outFileNamePrefix", out_prefix],
            log_to=os.path.join(logdir, "star_genome.log"))
        return True
    except Exception as exc:               # noqa: BLE001 -- advisory, never fatal
        LOG.warning("STAR --genomeLoad %s failed: %s", mode, exc)
        return False


# --- post-processing / stats -----------------------------------------------
def sort_index(bam: str, out_bam: str, *, threads: int = 4) -> str:
    """Coordinate-sort + index a BAM (STAR wrote it unsorted). Returns out_bam."""
    require_tools("samtools")
    if not nonempty(bam):
        raise ValueError(f"cannot sort an empty/missing BAM: {bam}")
    os.makedirs(os.path.dirname(os.path.abspath(out_bam)) or ".", exist_ok=True)
    n = str(int(threads))
    run(["samtools", "sort", "-@", n, "-o", out_bam, bam])
    run(["samtools", "index", "-@", n, out_bam])
    return out_bam


def index(bam: str, *, threads: int = 4) -> str:
    """Index an already coordinate-sorted BAM.

    Separate from `sort_index` because `umi_tools dedup` preserves the coordinate
    order of its input: re-sorting its output would be pure waste, but the old
    index no longer matches the new file and must be rebuilt.
    """
    require_tools("samtools")
    if not nonempty(bam):
        raise ValueError(f"cannot index an empty/missing BAM: {bam}")
    run(["samtools", "index", "-@", str(int(threads)), bam])
    return bam


def parse_log(star_log: str) -> dict:
    """Mapping stats from a STAR `Log.final.out`. Empty dict if there is no log.

    `unique + multi + unmapped` sums to 1 by construction: STAR reports the
    mapped classes only, and "unmapped" is whatever is left of the input.
    """
    if not star_log or not os.path.exists(star_log):
        return {}
    with open(star_log) as fh:
        txt = fh.read()

    def grab(pat: str, cast=int):
        m = re.search(pat + r"\s*\|\s*([\d.]+)", txt)
        return cast(m.group(1)) if m else None

    inp = grab("Number of input reads")
    uniq = grab("Uniquely mapped reads number")
    multi = grab("Number of reads mapped to multiple loci") or 0
    toomany = grab("Number of reads mapped to too many loci") or 0
    out: dict = {"n_input": inp, "n_unique": uniq, "n_multi": multi + toomany}
    if inp:
        mapped = (uniq + multi + toomany) if uniq is not None else None
        out["frac_mapped"] = mapped / inp if mapped is not None else None
        out["frac_multimapping"] = (multi + toomany) / inp
        out["frac_unique"] = uniq / inp if uniq is not None else None
        out["frac_unmapped"] = (1 - mapped / inp) if mapped is not None else None
    return out
