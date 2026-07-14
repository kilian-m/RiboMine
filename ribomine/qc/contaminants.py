"""Remove sequence-based contaminant reads before alignment.

Two adapter-agnostic filters (adapter- / primer-dimers and any other single-locus
pile-up are handled *after* mapping, by position, in `ribomine.qc.pileups` -- no
adapter sequences are hard-coded anywhere):

  1. STRUCTURED-RNA contaminants (rRNA, tRNA, snRNA, snoRNA, Mt) -- the dominant
     junk in ribo-seq and the reason multimapping looks huge (each rRNA has
     hundreds of genomic copies, so every rRNA read is a multimapper). Removed by
     aligning to a contaminant FASTA with `bowtie2 --very-sensitive-local` (local
     so the still-present 5' construct / 3' adapter are soft-clipped and the
     footprint core is matched). This uses a reference of contaminant *sequences*,
     not any hard-coded adapter.

  2. LOW-COMPLEXITY / homopolymer reads (poly-A / poly-G, simple repeats) --
     dropped by a Shannon-entropy / single-base-dominance test on the read.

Keeps the FULL (untrimmed) surviving reads -- architecture detection still needs
the soft-clippable 5' construct and 3' adapter.

`build_index()` is the other half: `ribomine setup` turns a contaminant FASTA
into the bowtie2 index this module needs.
"""
from __future__ import annotations

import contextlib
import functools
import glob
import gzip
import math
import os
import shutil
import subprocess
import tempfile
import time
from collections.abc import Iterator
from typing import TextIO

from ..config import Config
from ..utils import LOG, ToolError, count_fastq_reads, have, open_fastq, run

BASES = "ACGT"

# The gzipped FASTQs this module writes are throwaway intermediates (stage 4's
# clean FASTQ is deleted once the BAM exists) that run to multiple GB. Python's
# gzip defaults to compresslevel=9, which costs ~3-5x the CPU of level 6 for a
# couple of percent of size -- and is single-threaded on top. So: pigz when it is
# on PATH, else gzip at 6.
GZIP_LEVEL = 6
FLUSH_RECORDS = 8192          # records buffered per write() -- bounded memory

# bowtie2's `--un` file is written by the worker threads, and past a handful of them
# their records INTERLEAVE: a FASTQ whose 4-line records are spliced into each other.
# It is not a crash and not an error code -- it is a corrupt FASTQ that the next step
# maps without complaint, so the damage surfaces (if at all) as an inexplicable read
# count or a mangled alignment much later. The reads that survive this filter are the
# reads we map, so this is the one place in the pipeline where extra threads can
# silently change the DATA rather than just the speed.
#
# 8 is the cap, whatever `project.threads` says. bowtie2 is not the bottleneck of a
# sample (STAR is), and correctness is not negotiable for throughput.
BOWTIE2_MAX_THREADS = 8


def filter_fastq(in_fq: str, out_fq: str, cfg: Config, *, label: str = "",
                 threads: int = 8, log: str = "") -> dict:
    """Filter structured-RNA and low-complexity reads out of `in_fq` -> `out_fq`.

    Either path may be gzipped, independently of the other. bowtie2 reads .gz
    natively but `--un` always writes PLAIN FASTQ, so the intermediate is
    uncompressed whatever the input was, and the low-complexity screen decides the
    output's gz-ness from `out_fq`'s suffix alone.

    Returns the stats dict for `Sample.contam_json`.
    """
    label = label or os.path.basename(in_fq)
    if not os.path.isfile(in_fq):
        raise ValueError(f"input FASTQ not found: {in_fq}")

    if not cfg.get("contaminants.enabled", True):
        n_total = count_fastq_reads(in_fq)
        LOG.info("%s: contaminant filtering disabled; passing %s reads through",
                 label, f"{n_total:,}")
        _passthrough(in_fq, out_fq, threads=threads)
        return _stats(label, n_total, n_contaminant=0, n_low_complexity=0, n_kept=n_total)

    min_entropy = float(cfg.get("contaminants.min_entropy", 1.1))
    max_base_frac = float(cfg.get("contaminants.max_base_frac", 0.85))
    index = index_path(cfg)

    os.makedirs(os.path.dirname(os.path.abspath(out_fq)) or ".", exist_ok=True)

    n_bowtie_total = 0
    with _noncontam_tmp(out_fq) as tmp:
        if _index_exists(index):
            # bowtie2 --un writes plain FASTQ, so this temp is uncompressed. It
            # lives next to the output rather than in /tmp: it is as large as the
            # input, and a 500-sample run would otherwise fill the system tmpfs.
            n_bowtie_total, n_contaminant = _run_bowtie2(in_fq, index, tmp,
                                                         threads=threads, log=log)
            src = tmp                       # plain FASTQ, whatever in_fq was
        else:
            _warn_missing_index(index)
            n_contaminant = 0
            src = in_fq                     # may be .gz -- _screen handles both
        s = _screen(src, out_fq, min_entropy=min_entropy, max_base_frac=max_base_frac,
                    threads=threads)

    # Every input read is either aligned (contaminant) or written to --un, so the
    # total falls out of the two counters -- no need to decompress the whole input
    # a second time just to count its reads (a multi-GB pass in stage 4). bowtie2's
    # own "N reads; of these:" is the cross-check, and the fallback if it moves.
    n_total = n_bowtie_total or (n_contaminant + s["n_in"])

    # ... and that identity is also an exact integrity check on the `--un` file, for
    # free. bowtie2 says how many reads it read and how many it aligned; every other
    # read must be in --un. If the file is short, long, or spliced (see
    # BOWTIE2_MAX_THREADS), the arithmetic stops working -- and these are the reads we
    # are about to map, so a corrupt one must stop the sample, not quietly become a BAM.
    if n_bowtie_total:
        expected = n_bowtie_total - n_contaminant
        if s["n_in"] != expected:
            raise RuntimeError(
                f"{label}: bowtie2's --un output is corrupt -- it read {n_bowtie_total:,} "
                f"reads and aligned {n_contaminant:,}, so --un must hold {expected:,} "
                f"reads, but it holds {s['n_in']:,}. These are the reads that would be "
                f"mapped, so the sample is failed rather than filtered wrongly."
            )

    stats = _stats(label, n_total, n_contaminant=n_contaminant,
                   n_low_complexity=s["n_low_complexity"], n_kept=s["n_kept"])
    LOG.info("%s: %s/%s reads kept (%.0f%%) | removed: rRNA/tRNA/etc %.0f%%, "
             "low-complexity %.0f%%", label, f"{stats['n_kept']:,}", f"{n_total:,}",
             100 * stats["frac_kept"], 100 * stats["frac_contaminant_structured_rna"],
             100 * stats["frac_low_complexity"])
    return stats


# --- the index --------------------------------------------------------------
def build_index(fasta: str, cfg: Config) -> str:
    """bowtie2-build the contaminant FASTA into an index prefix; return the prefix.

    The prefix is the conventional one `filter_fastq` falls back to when
    `reference.contaminant_index` is unset: <workdir>/refs/<fasta-stem>.
    """
    fasta = os.path.abspath(fasta)
    if not os.path.isfile(fasta):
        raise ValueError(f"contaminant FASTA not found: {fasta}")

    prefix = _index_prefix(cfg, fasta)
    if _index_is_current(prefix, fasta):
        LOG.info("contaminant index already up to date: %s", prefix)
        return prefix

    exe = _bowtie2_build_bin()
    threads = max(1, int(cfg.get("project.threads", 8) or 1))
    LOG.info("building bowtie2 contaminant index: %s -> %s", fasta, prefix)
    t0 = time.time()

    # Build into a staging dir in the same directory and move the files in
    # afterwards: a killed bowtie2-build otherwise leaves a half-written
    # <prefix>.1.bt2 that the next run would find, trust, and align against.
    stage = tempfile.mkdtemp(dir=os.path.dirname(prefix), prefix=".bt2-build-")
    try:
        run([exe, "--threads", threads, fasta, os.path.join(stage, os.path.basename(prefix))])
        _install_index(stage, prefix)
    finally:
        shutil.rmtree(stage, ignore_errors=True)

    LOG.info("contaminant index built in %.1fs: %s", time.time() - t0, prefix)
    return prefix


def resolve_fasta(cfg: Config) -> str:
    """The contaminant FASTA to filter against.

    `reference.contaminant_fasta` wins; otherwise the human rRNA/tRNA/snRNA/snoRNA/Mt
    reference bundled with RiboMine. Bundling it means the common case needs no
    configuration at all -- and, more to the point, that the QC stage cannot silently
    run without a contaminant filter, which is the failure mode that makes a library
    look like 78% multimapping junk (each rRNA has hundreds of genomic copies, so
    every rRNA read maps to "too many loci").
    """
    fa = cfg.ref("contaminant_fasta")
    if fa:
        if not os.path.isfile(fa):
            raise ValueError(f"reference.contaminant_fasta not found: {fa}")
        return fa
    from ribomine import data

    return data.human_contaminants()


def index_path(cfg: Config) -> str | None:
    """The bowtie2 contaminant index prefix to filter against, or None if unbuilt.

    `reference.contaminant_index` wins. Otherwise the index built from
    `resolve_fasta` -- looked up where `ensure_index` / `ribomine setup` put it, so a
    user who ran `setup` (which does not rewrite their config) does not silently fall
    back to the low-complexity screen alone.
    """
    idx = cfg.ref("contaminant_index")
    if idx:
        return idx

    p = _index_prefix(cfg, resolve_fasta(cfg))
    if _index_exists(p):
        return p

    # Last resort: `ribomine setup --contaminant-fasta X` leaves its index in refs/.
    # Adopt it only if there is exactly one -- guessing between several would
    # silently filter against the wrong reference.
    refs = cfg.dir("refs")
    hits = [q for suf in (".1.bt2", ".1.bt2l")
            for q in glob.glob(os.path.join(refs, "*" + suf))
            # bowtie2-build also writes <prefix>.rev.1.bt2; that is not a prefix
            if not q.endswith(".rev" + suf)]
    if len(hits) == 1:
        return hits[0].rsplit(".1.bt2", 1)[0]
    return None


def ensure_index(cfg: Config) -> str | None:
    """Index prefix, building it from the resolved FASTA if it does not exist yet.

    MUST be called from the parent process, before the sample pool forks: N workers
    all discovering a missing index and all running bowtie2-build into the same
    prefix is corruption, not a slowdown. `pipeline.run` does this alongside the
    annotation index.

    It also LOGS which reference is in use. `reference.contaminant_fasta: null` means
    "the bundled human one", but in a JSON file a null reads as "nothing", so the run
    log has to say which it is -- otherwise the one thing a reader most needs to check
    (am I filtering against the right organism?) is the one thing they cannot see.
    """
    if not cfg["contaminants.enabled"]:
        LOG.warning("contaminants.enabled=false: rRNA/tRNA reads will NOT be removed. "
                    "Expect inflated multimapping -- each rRNA has hundreds of genomic copies.")
        return None
    idx = cfg.ref("contaminant_index")
    if idx:
        LOG.info("contaminant index: %s (configured)", idx)
        return idx
    fasta = resolve_fasta(cfg)
    bundled = not cfg.ref("contaminant_fasta")
    LOG.info("contaminant reference: %s%s", fasta,
             "  (bundled with RiboMine -- human)" if bundled else "  (from reference.contaminant_fasta)")
    return build_index(fasta, cfg)


# --- the two filters --------------------------------------------------------
def low_complexity(seq: str, *, min_entropy: float, max_base_frac: float) -> bool:
    """A homopolymer / simple repeat (poly-A, poly-G, ...), not a footprint."""
    if len(seq) < 5:
        return True
    counts = [seq.count(b) for b in BASES]
    n = sum(counts) or 1
    if max(counts) / n > max_base_frac:
        return True
    p = [c / n for c in counts if c]
    ent = -sum(pi * math.log2(pi) for pi in p)
    return ent < min_entropy


def _screen(fastq: str, out_fq: str, *, min_entropy: float, max_base_frac: float,
            threads: int = 1) -> dict:
    """Drop low-complexity / homopolymer reads; keep the rest at FULL length.

    Streams: one record in, one record out, never the file in memory. Input gz-ness
    comes from `fastq`'s suffix, output gz-ness from `out_fq`'s -- the two are
    independent (bowtie2's `--un` temp is plain even when the input was gzipped).
    """
    n_in = n_lowc = n_kept = 0
    # ... and write via a temp so a killed run cannot leave a truncated FASTQ that
    # `resume` (which only asks "is it non-empty?") would then map.
    partial = out_fq + ".partial"
    buf: list[str] = []
    try:
        with open_fastq(fastq, "rt") as fh, \
                _fastq_writer(partial, gz=out_fq.endswith(".gz"), threads=threads) as out:
            while True:
                h = fh.readline()
                if not h:
                    break
                seq = fh.readline().rstrip("\n")
                plus = fh.readline()
                qual = fh.readline().rstrip("\n")
                n_in += 1
                if low_complexity(seq.upper(), min_entropy=min_entropy,
                                  max_base_frac=max_base_frac):
                    n_lowc += 1
                    continue
                buf.append(h + seq + "\n" + plus + qual + "\n")
                n_kept += 1
                if len(buf) >= FLUSH_RECORDS:
                    out.write("".join(buf))
                    buf.clear()
            if buf:
                out.write("".join(buf))
        os.replace(partial, out_fq)
    finally:
        if os.path.exists(partial):
            os.remove(partial)
    return {"n_in": n_in, "n_low_complexity": n_lowc, "n_kept": n_kept}


def _run_bowtie2(fastq: str, index: str, noncontam_fq: str, *, threads: int,
                 log: str = "") -> tuple[int, int]:
    """bowtie2 --local; `--un` writes the non-contaminant reads (plain FASTQ, even
    for a .gz input, which bowtie2 reads natively). Returns (n_total, n_aligned).

    --very-sensitive-local finds the rRNA/tRNA core despite the appended adapter
    and RT modifications, without the over-aggressive short-match acceptance that
    a lowered --score-min would bring (which risks removing genuine footprints
    that chance-match a contaminant over ~20 nt).

    Threads are capped at BOWTIE2_MAX_THREADS: `--un` corrupts above that. See there.
    """
    bt2 = _bowtie2_bin()
    p_threads = max(1, min(int(threads), BOWTIE2_MAX_THREADS))
    if int(threads) > BOWTIE2_MAX_THREADS:
        LOG.debug("bowtie2: %d thread(s) requested, using %d (--un corrupts above it)",
                  int(threads), p_threads)
    cmd = [bt2, "--local", "--very-sensitive-local", "-p", str(p_threads),
           "-x", index, "-U", fastq, "--un", noncontam_fq, "-S", os.devnull]
    p = run(cmd, log_to=log or None)
    return _parse_bowtie2_summary(p.stderr or "")


def _parse_bowtie2_summary(stderr: str) -> tuple[int, int]:
    """"N reads; of these: ... M (P%) aligned ..." -- aligned = contaminant."""
    n_total = n_aligned = 0
    for line in stderr.splitlines():
        s = line.strip()
        if s.endswith("reads; of these:"):
            n_total = int(s.split()[0])
        elif "aligned exactly 1 time" in s or "aligned >1 times" in s:
            n_aligned += int(s.split()[0])
    return n_total, n_aligned


# --- gz plumbing ------------------------------------------------------------
@contextlib.contextmanager
def _fastq_writer(path: str, *, gz: bool, threads: int = 1) -> Iterator[TextIO]:
    """Text-mode FASTQ sink, gzipped iff `gz`.

    `gz` is passed in rather than sniffed off `path`, because `path` here is a
    `.partial` staging name whose suffix says nothing about the final file's
    compression.

    See GZIP_LEVEL: pigz across `threads` when available, else single-threaded
    gzip at level 6 -- never Python's default level 9 on a multi-GB intermediate.
    """
    if not gz:
        with open(path, "wt") as fh:
            yield fh
        return

    if not have("pigz"):
        with gzip.open(path, "wt", compresslevel=GZIP_LEVEL) as fh:
            yield fh
        return

    cmd = ["pigz", "-p", str(max(1, threads)), f"-{GZIP_LEVEL}", "-c"]
    with open(path, "wb") as raw:
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=raw,
                             stderr=subprocess.PIPE, text=True)
        try:
            yield p.stdin
        finally:
            # On the exception path this closes the pipe, lets pigz drain and
            # exit, and reaps it -- no orphan compressor writing into a file the
            # caller is about to delete.
            with contextlib.suppress(BrokenPipeError, OSError):
                p.stdin.close()
            err = p.stderr.read()
            p.stderr.close()
            rc = p.wait()
    if rc != 0:
        raise ToolError(cmd, rc, err)


@contextlib.contextmanager
def _noncontam_tmp(out_fq: str) -> Iterator[str]:
    """Path for bowtie2's `--un` output, removed on EVERY exit path.

    It is as large as the input (multi-GB in stage 4), so a leak here is not
    cosmetic: it fills the disk one failed sample at a time. It is created next to
    the output, not in /tmp, for the same reason.
    """
    d = os.path.dirname(os.path.abspath(out_fq))
    fd, tmp = tempfile.mkstemp(dir=d, suffix=".noncontam.fastq")
    os.close(fd)
    try:
        yield tmp
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def _passthrough(in_fq: str, out_fq: str, *, threads: int = 1) -> None:
    """contaminants.enabled=false: the output is the input, gz-ness respected."""
    if os.path.abspath(in_fq) == os.path.abspath(out_fq):
        return
    if in_fq.endswith(".gz") == out_fq.endswith(".gz"):
        try:
            if os.path.exists(out_fq):
                os.remove(out_fq)
            os.link(in_fq, out_fq)          # same filesystem: no second copy
            return
        except OSError:
            shutil.copyfile(in_fq, out_fq)
            return
    # gz-ness differs: recompress / decompress through the fast writer, streaming.
    partial = out_fq + ".partial"
    try:
        with open_fastq(in_fq, "rt") as fh, \
                _fastq_writer(partial, gz=out_fq.endswith(".gz"), threads=threads) as out:
            shutil.copyfileobj(fh, out)
        os.replace(partial, out_fq)
    finally:
        if os.path.exists(partial):
            os.remove(partial)


# --- plumbing ---------------------------------------------------------------
def _bowtie2_bin() -> str:
    exe = os.environ.get("BOWTIE2", "bowtie2")
    path = shutil.which(exe)
    if path is None:
        raise RuntimeError(
            f"bowtie2 not found (looked for {exe!r} on PATH; $BOWTIE2 overrides it). "
            f"Activate the RiboMine environment, or set contaminants.enabled=false "
            f"to skip the structured-RNA filter."
        )
    return path


def _bowtie2_build_bin() -> str:
    exe = os.environ.get("BOWTIE2_BUILD", "bowtie2-build")
    path = shutil.which(exe)
    if path is None:
        raise RuntimeError(
            f"bowtie2-build not found (looked for {exe!r} on PATH; $BOWTIE2_BUILD "
            f"overrides it), so the contaminant index cannot be built. Activate the "
            f"RiboMine environment: conda env create -f environment.yml && "
            f"conda activate ribomine"
        )
    return path


def _index_prefix(cfg: Config, fasta: str) -> str:
    """<workdir>/refs/<fasta-stem> -- the conventional index location."""
    stem = os.path.basename(fasta)
    if stem.endswith(".gz"):
        stem = stem[: -len(".gz")]
    stem = os.path.splitext(stem)[0]
    return os.path.join(cfg.dir("refs"), stem)


def _index_exists(index: str | None) -> bool:
    """bowtie2 writes .bt2 (or .bt2l for a large index)."""
    if not index:
        return False
    return any(os.path.exists(index + suf) for suf in (".1.bt2", ".1.bt2l"))


def _index_is_current(prefix: str, fasta: str) -> bool:
    """An index is reusable only if it is not older than the FASTA it came from."""
    fa_mtime = os.path.getmtime(fasta)
    for suf in (".1.bt2", ".1.bt2l"):
        p = prefix + suf
        if os.path.exists(p) and os.path.getmtime(p) >= fa_mtime:
            return True
    return False


def _install_index(stage: str, prefix: str) -> None:
    """Move a finished index out of the staging dir, the `.1.bt2` file LAST.

    `.1.bt2` is what _index_exists / _index_is_current test for, so it must be the
    last name to appear: interrupt the move and the next run rebuilds, rather than
    trusting an index whose other files are missing.
    """
    base = os.path.basename(prefix)
    outdir = os.path.dirname(prefix)
    files = sorted(f for f in os.listdir(stage) if f.startswith(base + "."))
    if not files:
        raise RuntimeError(f"bowtie2-build produced no index files for {prefix}")
    # exactly <base>.1.bt2 / .1.bt2l -- NOT <base>.rev.1.bt2, which is a different file
    last = [f for f in files if f in (base + ".1.bt2", base + ".1.bt2l")]
    for f in [f for f in files if f not in last] + last:
        os.replace(os.path.join(stage, f), os.path.join(outdir, f))


@functools.lru_cache(maxsize=None)
def _warn_missing_index(index: str | None) -> None:
    """Once per process: no index => low-complexity screen only, never a crash."""
    if index:
        LOG.warning("no bowtie2 index at %s; skipping the rRNA/tRNA filter "
                    "(low-complexity screen only)", index)
    else:
        LOG.warning("no contaminant index (reference.contaminant_index is unset and "
                    "none was found in <workdir>/refs); skipping the rRNA/tRNA filter "
                    "(low-complexity screen only). `ribomine setup` builds an index "
                    "from reference.contaminant_fasta.")


def _stats(label: str, n_input: int, *, n_contaminant: int, n_low_complexity: int,
           n_kept: int) -> dict:
    denom = max(n_input, 1)
    return {
        "label": label,
        "n_input": n_input,
        "n_contaminant_rRNA_tRNA_etc": n_contaminant,
        "n_low_complexity": n_low_complexity,
        "n_kept": n_kept,
        "frac_contaminant_structured_rna": round(n_contaminant / denom, 4),
        "frac_low_complexity": round(n_low_complexity / denom, 4),
        "frac_kept": round(n_kept / denom, 4),
    }
