"""STAR alignment: the permissive local pass and the final pass.

* `align_local` -- QC / architecture alignment of the untrimmed reads, with the
  fixed permissive settings `LOCAL_ARGS` (from `fqdissect.align`): non-genomic
  sequence (UMIs, RT additions, barcode, adapter) ends up in the soft clips,
  where fqdissect reads the architecture.
* `align_final` -- the deliverable BAM: the trimmed reads, aligned with the
  `mapping` section of the config. Mapping metrics and gene counts
  (`--quantMode GeneCounts`) come from this alignment; the permissive local one
  inflates multimapping.

Both write an unsorted BAM that `sort_index` sorts with samtools, because
STAR's own sorting is incompatible with the shared-memory genome.
"""
from __future__ import annotations

import os
import re
import shutil

from fqdissect.align import LOCAL_ARGS

from ..config import Config
from ..utils import LOG, nonempty, require_tools, rm, run

# --- genome-load mode ------------------------------------------------------
# Without a shared genome, each worker of the `project.jobs` process pool loads
# its own copy of the index (~28 GB for human). `load_genome` loads it once
# (`--genomeLoad LoadAndExit`), the alignments attach to it (`LoadAndKeep`), and
# `unload_genome` removes it (`Remove`).
#
# The effective mode reaches the workers through the environment, which both
# fork and spawn children inherit. Resolution order:
#   RIBOMINE_STAR_GENOME_LOAD -- set by `load_genome`: the mode this batch uses
#   STAR_GENOME_LOAD          -- the user's override
#   NoSharedMemory            -- default: every process loads its own index
STATE_ENV = "RIBOMINE_STAR_GENOME_LOAD"
GENOME_LOAD_ENV = "STAR_GENOME_LOAD"
DEFAULT_GENOME_LOAD = "NoSharedMemory"
SHARED = "LoadAndKeep"
# Modes an alignment may run with. The pipeline never selects LoadAndRemove: the
# first finished sample would unload the genome for the rest of the pool.
ALIGN_LOADS = ("NoSharedMemory", "LoadAndKeep", "LoadAndRemove")

_STAR_BAM = "Aligned.out.bam"     # what STAR writes
_OUT_BAM = "aligned.bam"          # what we hand on (Sample.local_bam / star_final)
_LOG_FINAL = "Log.final.out"
# Written by --quantMode GeneCounts next to the BAM, and the index's gene table
# (id, name, biotype). Both are read by ribomine.process.counts.
GENE_COUNTS = "ReadsPerGene.out.tab"
GENE_INFO = "geneInfo.tab"


# --- alignment -------------------------------------------------------------
def align_local(fastq: str, outdir: str, cfg: Config, *, threads: int = 8,
                log: str = "") -> str:
    """Permissive local alignment of the untrimmed sampled reads. Returns the BAM.

    UMIs, RT additions, barcode and adapter land in the soft clips, from which
    fqdissect reads the architecture.
    """
    return _align(fastq, outdir, cfg, args=LOCAL_ARGS, threads=threads, log=log,
                  what="local")


def align_final(fastq: str, outdir: str, cfg: Config, *, threads: int = 8,
                log: str = "") -> str:
    """Align the trimmed reads with the `mapping` config section. Returns the BAM.

    Also writes `ReadsPerGene.out.tab` (`--quantMode GeneCounts`) when the index
    carries an annotation; see `ribomine.process.counts` for what is counted.
    """
    args = [
        # Unsorted: STAR's sorter is incompatible with the shared-memory genome;
        # sort_index() sorts and indexes with samtools afterwards.
        "--outSAMtype", "BAM", "Unsorted",
        "--outSAMattributes", "nM", "MD", "NH",
        "--alignEndsType", str(cfg.get("mapping.align_ends_type", "Local")),
        "--outFilterMultimapNmax", str(int(cfg.get("mapping.multimap_nmax", 10))),
        "--outFilterMismatchNmax", str(int(cfg.get("mapping.mismatch_nmax", 3))),
        "--outFilterMismatchNoverLmax", str(float(cfg.get("mapping.mismatch_noverlmax", 0.1))),
        "--outFilterMatchNminOverLread", str(float(cfg.get("mapping.match_nmin_over_lread", 0.9))),
        "--outSJtype", "None",
    ]
    if has_annotation(cfg):
        args += ["--quantMode", "GeneCounts"]
    else:
        LOG.warning(
            "the STAR index was built without a GTF (--sjdbGTFfile), so STAR cannot "
            "count reads into genes: no read-count matrix will be written. Rebuild the "
            "index with --sjdbGTFfile %s to get one.", cfg.get("reference.gtf"))
    intron_max = int(cfg.get("mapping.align_intron_max", 0) or 0)
    if intron_max:      # 0 == STAR's own default (spliced); only pass it when set
        args += ["--alignIntronMax", str(intron_max)]
    args += [str(a) for a in (cfg.get("mapping.extra_args") or [])]
    return _align(fastq, outdir, cfg, args=args, threads=threads, log=log, what="final")


def has_annotation(cfg: Config) -> bool:
    """True if the STAR index was generated with a GTF.

    `--quantMode GeneCounts` needs the annotation in the index: supplying it at
    mapping time (`--sjdbGTFfile`) inserts junctions on the fly, which STAR
    refuses to do against a shared-memory genome. `geneInfo.tab` holds the gene
    count on its first line; an index built without a GTF has a 0 there.
    """
    gene_info = os.path.join(_index(cfg), GENE_INFO)   # a missing index raises here
    try:
        with open(gene_info) as fh:
            return int((fh.readline() or "0").strip()) > 0
    except (OSError, ValueError):
        return False


def _align(fastq: str, outdir: str, cfg: Config, *, args: list[str], threads: int,
           log: str, what: str) -> str:
    require_tools("STAR")
    if not nonempty(fastq):
        raise ValueError(f"no reads to align: {fastq}")
    index = _index(cfg)
    os.makedirs(outdir, exist_ok=True)
    # a killed STAR leaves _STARtmp behind and the next run refuses to start
    rm(os.path.join(outdir, "_STARtmp"))

    # genome_load() returns the mode the batch settled on (from the environment)
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


# --- shared-memory genome --------------------------------------------------
def genome_load() -> str:
    """The effective `--genomeLoad` mode for an alignment (resolution order above)."""
    mode = (os.environ.get(STATE_ENV) or os.environ.get(GENOME_LOAD_ENV)
            or DEFAULT_GENOME_LOAD)
    if mode not in ALIGN_LOADS:
        # an invalid override would make every STAR call fail; fall back instead
        LOG.warning("ignoring --genomeLoad %r (not one of %s); using %s",
                    mode, ", ".join(ALIGN_LOADS), DEFAULT_GENOME_LOAD)
        return DEFAULT_GENOME_LOAD
    return mode


def load_genome(cfg: Config) -> bool:
    """Load the STAR index into shared memory once, for the whole batch.

    Call before the process pool is created. Idempotent; never raises: if shared
    memory is unavailable (e.g. SHMALL/SHMMAX too small) it warns and falls back
    to NoSharedMemory, which gives the same output with one index copy per
    worker. Returns True if the batch shares the genome.
    """
    # already decided by an earlier call
    if STATE_ENV in os.environ:
        return os.environ[STATE_ENV] == SHARED

    want = os.environ.get(GENOME_LOAD_ENV)
    if want and want != SHARED:
        # the user asked for another mode
        LOG.info("STAR shared genome not used ($%s=%s)", GENOME_LOAD_ENV, want)
        os.environ[STATE_ENV] = want if want in ALIGN_LOADS else DEFAULT_GENOME_LOAD
        return False

    if _genome_op(cfg, "LoadAndExit", "_star_load."):
        os.environ[STATE_ENV] = SHARED
        LOG.info("STAR genome loaded into shared memory (one copy for the whole batch)")
        return True

    # fall back: one copy of the index per worker
    os.environ[STATE_ENV] = DEFAULT_GENOME_LOAD
    LOG.warning("STAR shared genome unavailable -- falling back to %s: each of the "
                "%d concurrent job(s) will load its OWN copy of the index. Lower "
                "project.jobs if the machine cannot hold that.",
                DEFAULT_GENOME_LOAD, int(cfg.get("project.jobs", 1) or 1))
    return False


def unload_genome(cfg: Config) -> bool:
    """Remove the shared-memory genome. Call it in a `finally`: a leaked segment
    holds its RAM until `STAR --genomeLoad Remove` is run by hand.

    Idempotent; never raises.
    """
    shared = os.environ.get(STATE_ENV) == SHARED
    # clear the decision: a second call is a no-op, and load_genome() can run again
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


def is_coordinate_sorted(bam: str) -> bool:
    import pysam

    try:
        with pysam.AlignmentFile(bam, "rb", check_sq=False) as fh:
            return fh.header.get("HD", {}).get("SO") == "coordinate"
    except Exception:  # noqa: BLE001 -- an unreadable BAM is "not sorted"
        return False


def ensure_sorted_indexed(bam: str, *, threads: int = 4) -> str:
    """Leave `bam` coordinate-sorted and indexed, in place. Idempotent.

    Applied to every BAM that is kept. Call it only after the analysis steps
    have read the BAM: the pile-up filter and the architecture profile read
    their input in file order.
    """
    if is_coordinate_sorted(bam) and nonempty(bam + ".bai"):
        return bam
    if not is_coordinate_sorted(bam):
        tmp = bam + ".sorted.tmp.bam"
        run(["samtools", "sort", "-@", str(int(threads)), "-o", tmp, bam])
        os.replace(tmp, bam)
    return index(bam, threads=threads)


def index(bam: str, *, threads: int = 4) -> str:
    """Index an already coordinate-sorted BAM.

    For a BAM that is sorted but has no valid index, such as a deduplicator's
    output, which keeps the coordinate order of its input.
    """
    require_tools("samtools")
    if not nonempty(bam):
        raise ValueError(f"cannot index an empty/missing BAM: {bam}")
    run(["samtools", "index", "-@", str(int(threads)), bam])
    return bam


def parse_log(star_log: str) -> dict:
    """Mapping stats from a STAR `Log.final.out`. Empty dict if there is no log.

    `frac_unmapped` is the remainder of the input, so unique + multi + unmapped
    sums to 1.
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
    out: dict = {
        "n_input": inp, "n_unique": uniq, "n_multi": multi + toomany,
        # STAR's input is the trimmed, contaminant-filtered FASTQ, so this is the mean
        # length of the candidate footprints (`mean_len_after_trim` still includes
        # the contaminants).
        "avg_input_len": grab("Average input read length", float),
        "avg_mapped_len": grab("Average mapped length", float),
    }
    if inp:
        mapped = (uniq + multi + toomany) if uniq is not None else None
        out["frac_mapped"] = mapped / inp if mapped is not None else None
        out["frac_multimapping"] = (multi + toomany) / inp
        out["frac_unique"] = uniq / inp if uniq is not None else None
        out["frac_unmapped"] = (1 - mapped / inp) if mapped is not None else None
    return out
