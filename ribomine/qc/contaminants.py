"""Remove sequence-based contaminant reads before alignment.

Two filters, neither of which needs an adapter sequence (adapter dimers and other
single-locus pile-ups are removed after mapping, in `ribomine.qc.pileups`):

  1. Structured RNA (rRNA, tRNA, snRNA, snoRNA, Mt): reads that align to a
     contaminant FASTA with `bowtie2 --very-sensitive-local`. Local alignment
     soft-clips the 5' construct and 3' adapter that are still on the read.
  2. Low-complexity reads (homopolymers, simple repeats): a Shannon-entropy and
     single-base-dominance test. The same pass drops reads longer than
     MAX_READ_LEN and malformed records.

Surviving reads are written unchanged (never trimmed here): in the QC stage,
read-architecture detection still needs the 5' construct and 3' adapter.
`build_index()` builds the bowtie2 index (`ribomine setup`).
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

# Compression level for the intermediate FASTQs written here (pigz if on PATH, else
# gzip). Python's default of 9 costs several times the CPU for a few percent of size.
GZIP_LEVEL = 6
FLUSH_RECORDS = 8192          # records buffered per write()

# Thread cap for bowtie2, whatever `project.threads` says. With many threads the
# records in the `--un` file can interleave, which silently corrupts the reads that
# are mapped next.
BOWTIE2_MAX_THREADS = 8

# Reads longer than this (nt) are dropped before STAR, whose short-read parser
# segfaults on multi-kilobase reads (long-read runs caught by the archive query).
# An untrimmed footprint read with construct, adapter and UMI is well under 150 nt.
MAX_READ_LEN = 300


def filter_fastq(in_fq: str, out_fq: str, cfg: Config, *, label: str = "",
                 threads: int = 8, log: str = "") -> dict:
    """Filter structured-RNA and low-complexity reads out of `in_fq` -> `out_fq`.

    Either path may be gzipped, independently of the other (decided by suffix).
    `log` is a file that receives bowtie2's output. Returns the stats dict for
    `Sample.contam_json`.
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
    n_nonbase = 0
    with _noncontam_tmp(out_fq) as tmp:
        if _index_exists(index):
            # `tmp` is bowtie2's --un output: plain FASTQ, as large as the input
            try:
                n_bowtie_total, n_contaminant = _run_bowtie2(in_fq, index, tmp,
                                                             threads=threads, log=log)
            except ToolError as exc:
                # bowtie2 aborts (SIGABRT) on reads that are not base-space, e.g. SOLiD
                # colour-space. Drop the reads it cannot parse and retry once, so the
                # sample gets a verdict on whatever survives.
                if "more quality values than read characters" not in (exc.stderr or ""):
                    raise
                with _base_space_tmp(in_fq, out_fq, threads=threads) as (clean_in, n_nonbase):
                    LOG.warning("%s: bowtie2 aborted -- dropped %s read(s) it could not parse "
                                "as base-space (non-ACGTN sequence; looks like SOLiD/colour-"
                                "space data mis-caught by the query) and retried.",
                                label, f"{n_nonbase:,}")
                    n_bowtie_total, n_contaminant = _run_bowtie2(clean_in, index, tmp,
                                                                 threads=threads, log=log)
            src = tmp                       # plain FASTQ, whatever in_fq was
        else:
            _warn_missing_index(index)
            n_contaminant = 0
            src = in_fq                     # may be .gz -- _screen handles both
        s = _screen(src, out_fq, min_entropy=min_entropy, max_base_frac=max_base_frac,
                    threads=threads)

    # Input total without a second pass over the input: bowtie2's own read count plus
    # any reads dropped before the retry; without an index, the reads _screen saw.
    n_total = (n_bowtie_total + n_nonbase) or (n_contaminant + s["n_in"])

    # Integrity check on --un: it must hold exactly the reads bowtie2 read but did not
    # align. A mismatch means a corrupt file (see BOWTIE2_MAX_THREADS); fail the sample.
    if n_bowtie_total:
        expected = n_bowtie_total - n_contaminant
        if s["n_in"] != expected:
            raise RuntimeError(
                f"{label}: bowtie2's --un output is corrupt -- it read {n_bowtie_total:,} "
                f"reads and aligned {n_contaminant:,}, so --un must hold {expected:,} "
                f"reads, but it holds {s['n_in']:,}. These are the reads that would be "
                f"mapped, so the sample is failed rather than filtered wrongly."
            )

    if s.get("n_overlong"):
        LOG.warning("%s: dropped %s read(s) longer than %d nt -- STAR's short-read parser "
                    "overflows on them (fatal error + segfault, exit 139). This run looks "
                    "like long-read data mis-caught by the query; the short reads are "
                    "filtered and mapped so QC can still judge it.",
                    label, f"{s['n_overlong']:,}", MAX_READ_LEN)
    if s.get("n_malformed"):
        LOG.warning("%s: dropped %s read(s) whose quality string length != the sequence "
                    "length -- a corrupt FASTQ record that makes STAR fatal-error and "
                    "segfault (exit 139). The rest of the run is filtered and mapped "
                    "normally.", label, f"{s['n_malformed']:,}")

    # "malformed" = reads dropped as unparseable, by the bowtie2 retry or by _screen
    stats = _stats(label, n_total, n_contaminant=n_contaminant,
                   n_low_complexity=s["n_low_complexity"], n_kept=s["n_kept"],
                   n_malformed=s.get("n_malformed", 0) + n_nonbase,
                   n_overlong=s.get("n_overlong", 0))
    LOG.info("%s: %s/%s reads kept (%.0f%%) | removed: rRNA/tRNA/etc %.0f%%, "
             "low-complexity %.0f%%", label, f"{stats['n_kept']:,}", f"{n_total:,}",
             100 * stats["frac_kept"], 100 * stats["frac_contaminant_structured_rna"],
             100 * stats["frac_low_complexity"])
    return stats


# --- index ------------------------------------------------------------------
def build_index(fasta: str, cfg: Config) -> str:
    """bowtie2-build the contaminant FASTA; return the index prefix,
    <workdir>/refs/<fasta-stem>. An up-to-date index is reused.
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

    # Build in a staging dir and move the files in afterwards, so a killed build
    # leaves no partial index that a later run would use.
    stage = tempfile.mkdtemp(dir=os.path.dirname(prefix), prefix=".bt2-build-")
    try:
        run([exe, "--threads", threads, fasta, os.path.join(stage, os.path.basename(prefix))])
        _install_index(stage, prefix)
    finally:
        shutil.rmtree(stage, ignore_errors=True)

    LOG.info("contaminant index built in %.1fs: %s", time.time() - t0, prefix)
    return prefix


def resolve_fasta(cfg: Config) -> str:
    """The contaminant FASTA to filter against: `reference.contaminant_fasta` if set,
    otherwise the human rRNA/tRNA/snRNA/snoRNA/Mt reference bundled with RiboMine.
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

    `reference.contaminant_index` if set; otherwise the index built from
    `resolve_fasta` at its conventional location under <workdir>/refs.
    """
    idx = cfg.ref("contaminant_index")
    if idx:
        return idx

    p = _index_prefix(cfg, resolve_fasta(cfg))
    if _index_exists(p):
        return p

    # Last resort: an index left in refs/ by `ribomine setup --contaminant-fasta X`,
    # adopted only if it is the only one there.
    refs = cfg.dir("refs")
    hits = [q for suf in (".1.bt2", ".1.bt2l")
            for q in glob.glob(os.path.join(refs, "*" + suf))
            # <prefix>.rev.1.bt2 is not a prefix
            if not q.endswith(".rev" + suf)]
    if len(hits) == 1:
        return hits[0].rsplit(".1.bt2", 1)[0]
    return None


def ensure_index(cfg: Config) -> str | None:
    """Index prefix, building it from the resolved FASTA if needed; None when the
    filter is disabled. Logs which reference is in use.

    Call it from the parent process before the sample pool starts: concurrent
    bowtie2-build runs into the same prefix would corrupt the index.
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


# --- filters ----------------------------------------------------------------
def low_complexity(seq: str, *, min_entropy: float, max_base_frac: float) -> bool:
    """True for a homopolymer or simple repeat: shorter than 5 nt, one base above
    `max_base_frac`, or base entropy below `min_entropy` bits."""
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
    """Stream `fastq` -> `out_fq`, dropping over-length (> MAX_READ_LEN), malformed
    (len(qual) != len(seq)) and low-complexity reads; the rest are kept unchanged.

    The first two would crash STAR. All drops are counted in `n_in`, so
    filter_fastq's --un integrity check balances. Returns the counts.
    """
    n_in = n_lowc = n_kept = n_malformed = n_overlong = 0
    # write via a temp so a killed run leaves no truncated FASTQ for a resume to map
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
                if len(seq) > MAX_READ_LEN:
                    n_overlong += 1
                    continue
                if len(seq) != len(qual):
                    n_malformed += 1
                    continue
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
    return {"n_in": n_in, "n_low_complexity": n_lowc, "n_kept": n_kept,
            "n_malformed": n_malformed, "n_overlong": n_overlong}


def _run_bowtie2(fastq: str, index: str, noncontam_fq: str, *, threads: int,
                 log: str = "") -> tuple[int, int]:
    """Align to the contaminant index; `--un` writes the unaligned reads to
    `noncontam_fq` as plain FASTQ. Returns (n_total, n_aligned).

    --very-sensitive-local matches the rRNA/tRNA core despite the adapter, without
    lowering --score-min, which would also remove footprints with a short chance
    match. Threads are capped at BOWTIE2_MAX_THREADS.
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
    """(n_total, n_aligned) from bowtie2's stderr summary; aligned = contaminant."""
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
    """Text-mode FASTQ sink, gzipped iff `gz` (`path` may be a `.partial` name, so
    its suffix is not used). pigz across `threads` when available, else gzip, both
    at GZIP_LEVEL.
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
            # also on the exception path: close the pipe and reap pigz
            with contextlib.suppress(BrokenPipeError, OSError):
                p.stdin.close()
            err = p.stderr.read()
            p.stderr.close()
            rc = p.wait()
    if rc != 0:
        raise ToolError(cmd, rc, err)


@contextlib.contextmanager
def _noncontam_tmp(out_fq: str) -> Iterator[str]:
    """Temp path for bowtie2's `--un` output, removed on every exit path. It is as
    large as the input, so it is created next to the output rather than in /tmp.
    """
    d = os.path.dirname(os.path.abspath(out_fq))
    fd, tmp = tempfile.mkstemp(dir=d, suffix=".noncontam.fastq")
    os.close(fd)
    try:
        yield tmp
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


# Translation table that deletes ACGTN; anything left in a sequence is a character
# bowtie2 cannot parse (e.g. the digits of SOLiD colour-space reads).
_NON_BASE = str.maketrans("", "", "ACGTNacgtn")


def _copy_base_space(in_fq: str, out_fq: str) -> int:
    """Stream `in_fq` -> `out_fq` (plain FASTQ), dropping records bowtie2 cannot
    parse: characters outside ACGTN, or len(qual) != len(seq). Returns the number
    dropped. Used only for the retry after bowtie2 has aborted.
    """
    n_dropped = 0
    buf: list[str] = []
    with open_fastq(in_fq, "rt") as fh, open(out_fq, "wt") as out:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            s = seq.rstrip("\n")
            if len(s) != len(qual.rstrip("\n")) or s.translate(_NON_BASE):
                n_dropped += 1
                continue
            buf.append(h + seq + plus + qual)
            if len(buf) >= FLUSH_RECORDS:
                out.write("".join(buf))
                buf.clear()
        if buf:
            out.write("".join(buf))
    return n_dropped


@contextlib.contextmanager
def _base_space_tmp(in_fq: str, out_fq: str, *, threads: int = 1) -> Iterator[tuple[str, int]]:
    """Temp copy of `in_fq` without the records bowtie2 cannot parse (see
    `_copy_base_space`). Yields (path, n_dropped); created next to the output and
    removed on every exit path."""
    d = os.path.dirname(os.path.abspath(out_fq))
    fd, tmp = tempfile.mkstemp(dir=d, suffix=".basespace.fastq")
    os.close(fd)
    try:
        n_dropped = _copy_base_space(in_fq, tmp)
        yield tmp, n_dropped
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def _passthrough(in_fq: str, out_fq: str, *, threads: int = 1) -> None:
    """contaminants.enabled=false: copy (or hard-link) the input to the output,
    recompressing if their gz suffixes differ."""
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
    # gz suffixes differ: recompress or decompress, streaming
    partial = out_fq + ".partial"
    try:
        with open_fastq(in_fq, "rt") as fh, \
                _fastq_writer(partial, gz=out_fq.endswith(".gz"), threads=threads) as out:
            shutil.copyfileobj(fh, out)
        os.replace(partial, out_fq)
    finally:
        if os.path.exists(partial):
            os.remove(partial)


# --- helpers ----------------------------------------------------------------
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
    """The conventional index prefix, <workdir>/refs/<fasta-stem>."""
    stem = os.path.basename(fasta)
    if stem.endswith(".gz"):
        stem = stem[: -len(".gz")]
    stem = os.path.splitext(stem)[0]
    return os.path.join(cfg.dir("refs"), stem)


def _index_exists(index: str | None) -> bool:
    """True if <index>.1.bt2 (or .1.bt2l, for a large index) exists."""
    if not index:
        return False
    return any(os.path.exists(index + suf) for suf in (".1.bt2", ".1.bt2l"))


def _index_is_current(prefix: str, fasta: str) -> bool:
    """True if the index exists and is not older than its FASTA."""
    fa_mtime = os.path.getmtime(fasta)
    for suf in (".1.bt2", ".1.bt2l"):
        p = prefix + suf
        if os.path.exists(p) and os.path.getmtime(p) >= fa_mtime:
            return True
    return False


def _install_index(stage: str, prefix: str) -> None:
    """Move a finished index out of the staging dir, the `.1.bt2` file last.

    `.1.bt2` is the file _index_exists / _index_is_current test for, so an
    interrupted move leads to a rebuild, not to an incomplete index being used.
    """
    base = os.path.basename(prefix)
    outdir = os.path.dirname(prefix)
    files = sorted(f for f in os.listdir(stage) if f.startswith(base + "."))
    if not files:
        raise RuntimeError(f"bowtie2-build produced no index files for {prefix}")
    # exactly <base>.1.bt2 / .1.bt2l, not <base>.rev.1.bt2
    last = [f for f in files if f in (base + ".1.bt2", base + ".1.bt2l")]
    for f in [f for f in files if f not in last] + last:
        os.replace(os.path.join(stage, f), os.path.join(outdir, f))


@functools.lru_cache(maxsize=None)
def _warn_missing_index(index: str | None) -> None:
    """Warn, once per process, that only the low-complexity screen will run."""
    if index:
        LOG.warning("no bowtie2 index at %s; skipping the rRNA/tRNA filter "
                    "(low-complexity screen only)", index)
    else:
        LOG.warning("no contaminant index (reference.contaminant_index is unset and "
                    "none was found in <workdir>/refs); skipping the rRNA/tRNA filter "
                    "(low-complexity screen only). `ribomine setup` builds an index "
                    "from reference.contaminant_fasta.")


def _stats(label: str, n_input: int, *, n_contaminant: int, n_low_complexity: int,
           n_kept: int, n_malformed: int = 0, n_overlong: int = 0) -> dict:
    denom = max(n_input, 1)
    return {
        "label": label,
        "n_input": n_input,
        "n_contaminant_rRNA_tRNA_etc": n_contaminant,
        "n_low_complexity": n_low_complexity,
        "n_malformed": n_malformed,
        "n_overlong": n_overlong,
        "n_kept": n_kept,
        "frac_contaminant_structured_rna": round(n_contaminant / denom, 4),
        "frac_low_complexity": round(n_low_complexity / denom, 4),
        "frac_malformed": round(n_malformed / denom, 4),
        "frac_overlong": round(n_overlong / denom, 4),
        "frac_kept": round(n_kept / denom, 4),
    }
