"""Getting reads out of the archive -- the slow part, so it is the measured part.

Two very different jobs live here.

**The read sample (QC stage).** Architecture and periodicity are properties of
*every* read, so screening a run does not need the run. We stream the first
`scan` reads straight off ENA's gzipped FASTQ over HTTPS, reservoir-sample `n`
of them, and drop the connection -- nothing is stored, and a 20 GB run costs a
few seconds and ~50 MB of transfer. gzip is a stream format, so a prefix of the
file decompresses to a prefix of the reads; we simply stop reading.

  What the sample *is*: a uniform random sample of the reads it scanned (the
  kept records are shuffled), but only over the run's first `scan` reads. ENA
  stores FASTQ in spot order -- the order reads came off the sequencer, which is
  random with respect to content -- so that prefix is representative. The
  exception is a **sorted or collapsed** deposit, where the prefix is a biased
  slice; that trips a warning, and `qc.scan_reads: 0` reservoir-samples the whole
  run instead.

**The full dataset (processing stage).** ENA over HTTPS with parallel connections
wins, and it is not close -- see `docs/DOWNLOAD.md` for the measurements. So there
is no route to choose and nothing to tune: `ena_https` is simply the route.

  1. `ena_https`  -- ENA serves the submitter's own FASTQ directly. No .sra
                     container to convert, and `aria2c -x16` opens 16 ranged
                     connections, which is what actually beats the throttle: a
                     single HTTPS stream is capped server-side, so parallelism
                     buys far more than bandwidth does.

The two routes below it exist for *availability*, not speed: ENA does not mirror
every run (recent releases, dbGaP), and a run we cannot fetch is a run we cannot
mine. They are tried, in order, only when the one above them cannot serve the
accession at all. Nobody picks between them.

  2. `aws_odp`    -- SRA's Open Data mirror on S3 (`s3://sra-pub-run-odp`),
                     no-sign-request and no egress charge to the downloader.
                     Delivers a .sra that `fasterq-dump` converts locally.
  3. `prefetch`   -- the SRA toolkit itself. Slowest and the most fragile, but it
                     is the only route that always exists.

Every route ends at the same artefact: one gzipped single-end FASTQ at
`out_fastq_gz`, so the caller never has to care which one ran.
"""
from __future__ import annotations

import gzip
import os
import random
import shutil
import subprocess
import sys
import time

from ..config import Config
from ..utils import LOG, ToolError, have, nonempty, rm, run
from . import metadata

# ---------------------------------------------------------------------------
# the read sample
# ---------------------------------------------------------------------------
def _open_stream(src: str):
    """(binary FASTQ text stream, process-or-None). `src` is a URL or a path.

    NOTE the absence of `--retry`. curl is writing into a PIPE, and a curl retry
    restarts the transfer from byte 0 -- into that same pipe. The reader would then
    get the head of the file spliced into the middle of the gzip stream: not a
    recovered download, a corrupt one. On a pipe the retry has to re-open the whole
    stream (which is what `sample_reads` does), so curl is told to fail honestly and
    say why (`--show-error`) instead of trying to fix it here.
    """
    if src.startswith("http"):
        proc = subprocess.Popen(
            ["curl", "-sL", "--fail", "--show-error", "--max-time", "7200", src],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return gzip.GzipFile(fileobj=proc.stdout), proc
    if src.endswith(".gz"):
        return gzip.open(src, "rb"), None
    return open(src, "rb"), None


def _reservoir(stream, scan: float, n: int, rng: random.Random):
    """Uniform sample of n records over the first `scan` reads.

    Also watches whether consecutive reads are monotone in their leading 8-mer or
    in length -- the signature of a sequence-sorted or collapsed deposit, whose
    *prefix* would not represent the run.
    """
    keep: list[list[bytes]] = []
    seen = 0
    prev_lead = prev_len = None
    mono_seq = cmp_n = 0
    len_up = len_down = 0        # strict transitions only; ties (a constant read
                                 # length, the common raw case) are not evidence
    while seen < scan:
        rec = [stream.readline() for _ in range(4)]
        if not rec[0]:
            break
        seen += 1
        seq = rec[1].rstrip()
        lead = seq[:8]
        if prev_lead is not None:
            mono_seq += (lead >= prev_lead)
            len_up += (len(seq) > prev_len)
            len_down += (len(seq) < prev_len)
            cmp_n += 1
        prev_lead, prev_len = lead, len(seq)
        if len(keep) < n:
            keep.append(rec)
        else:
            j = rng.randrange(seen)
            if j < n:
                keep[j] = rec

    tail = bool(stream.readline()) if seen >= scan else False
    nontie = len_up + len_down
    order = {
        "tail_unscanned": tail,
        "cmp_n": cmp_n,
        "seq_monotone": (max(mono_seq, cmp_n - mono_seq) / cmp_n) if cmp_n else 0.0,
        "len_nontie": nontie,
        "len_monotone": (max(len_up, len_down) / nontie) if nontie else 0.0,
    }
    return keep, seen, order


def _stream_once(src: str, limit: float, n: int, rng: random.Random):
    """One attempt: open the stream, reservoir-sample it, and always tear it down.

    Raises OSError if the transfer failed -- which the caller retries on a fresh
    connection.

    A FAILED TRANSFER DOES NOT LOOK LIKE AN ERROR FROM IN HERE, and that is the whole
    point of this function. ENA answers a burst of concurrent requests with 403; curl
    then writes nothing at all, so the gzip stream is simply EMPTY, the reservoir reads
    zero records and returns perfectly normally. Nothing raises. The sample died as "no
    reads obtained" and was never retried: 20 of 100 runs, the first time this pipeline
    asked ENA for 8 samples at once. A mid-transfer SSL error (curl 56) is worse -- the
    reservoir gets SOME reads and returns them, and a short sample is a BIASED sample
    that nothing downstream can tell apart from a good one.

    So the transfer is judged by curl's exit status, not by whether bytes arrived. curl
    still running when we are done means we stopped early, by design; curl already gone
    with a non-zero status means the reads we just took are not the reads we asked for.
    """
    stream, proc = _open_stream(src)
    try:
        got = _reservoir(stream, limit, n, rng)
        rc = proc.poll() if proc is not None else 0
        if rc:
            err = (proc.stderr.read() or b"").decode(errors="replace").strip() \
                if proc.stderr else ""
            raise OSError(f"the read stream failed: curl exited {rc}"
                          + (f" ({err[:120]})" if err else ""))
        return got
    finally:
        try:
            stream.close()
        except Exception:  # noqa: BLE001 -- we are already unwinding; the read result
            pass           # (or the error) is what matters, not the closing of a pipe
        if proc is not None:
            for pipe in (proc.stdout, proc.stderr):
                if pipe:
                    pipe.close()
            proc.terminate()
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                proc.kill()


def sample_reads(source: str, out_fastq: str, *, n: int = 200_000,
                 scan: int = 1_000_000, seed: int = 20260712,
                 tries: int = 4, backoff_s: float = 5.0) -> dict:
    """Reservoir-sample `n` reads from an accession (streamed) or a local FASTQ.

    A dropped connection re-opens the stream and samples again from the start, up to
    `tries` times. It has to be re-opened rather than resumed: the sample is a uniform
    draw over the run's first `scan` reads, and a resume would splice two draws with
    different denominators together. Re-reading a few tens of MB is cheap; getting the
    sample subtly wrong is not.

    The retry is here, and not in curl, on purpose -- see `_open_stream`. Without it a
    single truncated stream permanently fails the run: measured at ~10% of a 20-run
    batch when four samples stream at once, which is a tenth of a mining cohort lost to
    a transient.
    """
    if os.path.exists(source):
        src = source
    else:
        urls = metadata.fastq_urls(source)
        if not urls:
            # No ENA mirror. fastq-dump can stream a bounded slice straight from
            # NCBI without materialising the run, which is the only other way to
            # get a read sample cheaply.
            return _sample_via_sra(source, out_fastq, n=n)
        src = urls[0]
        LOG.debug("streaming %s", src)

    limit = float("inf") if scan <= 0 else scan
    local = os.path.exists(src)
    for attempt in range(1, tries + 1):
        rng = random.Random(seed)          # the same draw every attempt, by construction
        try:
            keep, seen, order = _stream_once(src, limit, n, rng)
            break
        except (EOFError, OSError, gzip.BadGzipFile) as exc:
            # A truncated gzip stream is what a dropped HTTPS connection looks like from
            # in here. A local file that does this is genuinely corrupt, so it is not
            # retried -- re-reading it would fail identically, four times.
            if local or attempt == tries:
                raise RuntimeError(
                    f"{source}: could not stream a read sample ({exc})") from exc
            wait = backoff_s * attempt
            LOG.warning("[%s] read-sample stream broke (%s); retrying in %.0fs (%d/%d)",
                        source, exc, wait, attempt, tries - 1)
            time.sleep(wait)

    if not keep:
        raise RuntimeError(f"no reads obtained from {source}")

    rng.shuffle(keep)
    os.makedirs(os.path.dirname(os.path.abspath(out_fastq)) or ".", exist_ok=True)
    with open(out_fastq, "wb") as out:
        for rec in keep:
            out.write(b"".join(rec))

    warn = ""
    seq_sorted = order["seq_monotone"] > 0.97
    len_sorted = order["len_nontie"] >= 10 and order["len_monotone"] > 0.98
    if order["tail_unscanned"] and order["cmp_n"] > 1000 and (seq_sorted or len_sorted):
        by = "sequence" if seq_sorted else "length"
        warn = (f"the deposit looks sorted by {by} and a tail was left unscanned, so this "
                f"prefix sample may be biased; set qc.scan_reads=0 for a whole-run sample")
    return {"source": src, "n_sampled": len(keep), "n_scanned": seen,
            "sorted_warning": warn}


def _sample_via_sra(acc: str, out_fastq: str, *, n: int) -> dict:
    """No ENA fastq mirror: pull a bounded slice through the SRA toolkit instead.

    `fastq-dump -X n` stops after n spots, so this stays cheap even for a huge run.
    """
    if not have("fastq-dump"):
        raise RuntimeError(
            f"{acc} has no ENA FASTQ mirror and fastq-dump is not installed; "
            f"install sra-tools (it is in environment.yml)")
    LOG.info("[%s] no ENA mirror; sampling %d reads via fastq-dump", acc, n)
    d = os.path.dirname(os.path.abspath(out_fastq))
    run(["fastq-dump", "-X", str(n), "--split-spot", "--skip-technical",
         "-O", d, acc], capture=True)
    produced = os.path.join(d, f"{acc}.fastq")
    if not nonempty(produced):
        raise RuntimeError(f"fastq-dump produced nothing for {acc}")
    if produced != out_fastq:
        os.replace(produced, out_fastq)
    with open(out_fastq, "rb") as fh:
        got = sum(1 for i, _ in enumerate(fh)) // 4
    return {"source": f"sra:{acc}", "n_sampled": got, "n_scanned": got, "sorted_warning": ""}


# ---------------------------------------------------------------------------
# the full dataset: one function per route, all landing on the same artefact
# ---------------------------------------------------------------------------
class RouteUnavailable(RuntimeError):
    """This route cannot serve this accession (no mirror, tool missing). Try the next."""


def _scratch(sra: str, fallback: str) -> str:
    """fasterq-dump wants ~10x the .sra size in temp space. /dev/shm is RAM and by
    far the fastest, but only if the run actually fits -- filling it would take the
    machine down, so fall back to the workdir when it does not."""
    need = os.path.getsize(sra) * 12
    shm = "/dev/shm"
    try:
        if os.path.isdir(shm):
            st = os.statvfs(shm)
            if st.f_bavail * st.f_frsize > need * 2:      # leave the same again free
                d = os.path.join(shm, "ribomine")
                os.makedirs(d, exist_ok=True)
                return d
    except OSError:
        pass
    return fallback


def _finalise(src_fastq: str, out_gz: str, *, threads: int = 4) -> None:
    """gzip a plain FASTQ into place (pigz when available -- it is ~4x faster)."""
    os.makedirs(os.path.dirname(os.path.abspath(out_gz)) or ".", exist_ok=True)
    tmp = out_gz + ".part"
    if have("pigz"):
        with open(tmp, "wb") as out:
            run(["pigz", "-c", "-p", str(threads), src_fastq], stdout=out)
    else:
        with open(src_fastq, "rb") as fin, gzip.open(tmp, "wb") as fout:
            shutil.copyfileobj(fin, fout, length=1 << 22)
    os.replace(tmp, out_gz)
    rm(src_fastq)


def route_ena_https(acc: str, out_gz: str, cfg: Config, log: str = "") -> dict:
    """ENA's own gzipped FASTQ, pulled with N ranged connections.

    ENA hosts what the submitter uploaded, so there is no .sra to convert -- the
    bytes on the wire are the bytes we want. The win comes from `-x/-s`: a single
    HTTPS stream is throttled server-side, so opening 16 of them multiplies the
    throughput almost linearly until the link saturates.
    """
    row = metadata.filereport(acc, ["fastq_ftp", "fastq_md5"])
    paths = [p for p in (row.get("fastq_ftp") or "").split(";") if p]
    md5s = [m for m in (row.get("fastq_md5") or "").split(";") if m]
    if not paths:
        raise RouteUnavailable(f"{acc}: no ENA FASTQ mirror")
    # a paired deposit lists _1 and _2; the footprint read is R1. md5 is positionally
    # aligned with the url list, so pick the checksum by the same index.
    i = 0
    if len(paths) > 1:
        i = next((j for j, p in enumerate(paths) if p.endswith("_1.fastq.gz")), 0)
    url = "https://" + paths[i]
    want_md5 = md5s[i] if i < len(md5s) else ""

    conn = int(cfg["download.connections"])
    tmp = os.path.join(cfg.tmpdir, f"{acc}.ena.fastq.gz")
    rm(tmp)

    if have("aria2c"):
        run(["aria2c", "-x", str(conn), "-s", str(conn), "-k", "8M",
             "--max-tries", "3", "--retry-wait", "5", "--auto-file-renaming=false",
             "--allow-overwrite=true", "--console-log-level=warn", "--summary-interval=0",
             "-d", os.path.dirname(tmp), "-o", os.path.basename(tmp), url], log_to=log)
    else:
        LOG.warning("aria2c not installed -- falling back to a single curl stream. ENA caps "
                    "one stream at ~12 MB/s, so this is several times slower; install aria2.")
        run(["curl", "-sL", "--fail", "--retry", "3", "-o", tmp, url], log_to=log)

    if not nonempty(tmp):
        raise RouteUnavailable(f"{acc}: ENA download produced nothing")
    # ENA hands us the checksum in the same call as the URL, so there is no excuse for
    # not checking it. A truncated download otherwise surfaces as a mangled read count
    # ten steps later, where nobody would suspect the transfer.
    if want_md5:
        got = _md5(tmp)
        if got != want_md5:
            rm(tmp)
            raise RuntimeError(f"{acc}: md5 mismatch (got {got}, ENA says {want_md5})")

    os.makedirs(os.path.dirname(os.path.abspath(out_gz)) or ".", exist_ok=True)
    os.replace(tmp, out_gz)
    return {"url": url, "md5_verified": bool(want_md5)}


def _md5(path: str) -> str:
    import hashlib

    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


AWS_ODP = "https://sra-pub-run-odp.s3.amazonaws.com/sra/{acc}/{acc}"


def route_aws_odp(acc: str, out_gz: str, cfg: Config, log: str = "") -> dict:
    """SRA's Open Data mirror on S3.

    The bucket is public, so it needs no AWS account, no credentials and no `aws`
    CLI -- it is a plain anonymous HTTPS URL, which means the same multi-connection
    `aria2c` that makes the ENA route fast works here too. (`prefetch` resolves to
    this very URL and then fetches it over a single throttled stream, which is why
    it is the slowest route and not the first.)

    Two things to know. The mirror serves a `.sra` container, so `fasterq-dump` has
    to convert it locally -- that cost is CPU and scratch disk, not network. And it
    is the *full-quality* copy: SRA Lite (`.sralite`) substitutes a flat fake
    quality score for every base, which would quietly destroy any quality-aware
    step downstream. This route never touches it.
    """
    if not have("fasterq-dump"):
        raise RouteUnavailable("fasterq-dump not installed (needed to convert the .sra)")
    url = AWS_ODP.format(acc=acc)
    sra = os.path.join(cfg.tmpdir, f"{acc}.sra")
    rm(sra)
    conn = int(cfg["download.connections"])
    try:
        if have("aria2c"):
            run(["aria2c", "-x", str(conn), "-s", str(conn), "-k", "8M",
                 "--max-tries", "2", "--auto-file-renaming=false", "--allow-overwrite=true",
                 "--console-log-level=warn", "--summary-interval=0",
                 "-d", os.path.dirname(sra), "-o", os.path.basename(sra), url], log_to=log)
        else:
            run(["curl", "-sL", "--fail", "-o", sra, url], log_to=log)
    except ToolError as exc:
        raise RouteUnavailable(f"{acc}: not on the S3 Open Data mirror") from exc
    if not nonempty(sra):
        raise RouteUnavailable(f"{acc}: not on the S3 Open Data mirror")

    fq = _fasterq(sra, acc, cfg, log)
    if not cfg["keep.sra"]:
        rm(sra)
    _finalise(fq, out_gz, threads=cfg["project.threads"])
    return {"url": url}


def route_prefetch(acc: str, out_gz: str, cfg: Config, log: str = "") -> dict:
    """The SRA toolkit. Always available, never the fastest."""
    if not have("prefetch") or not have("fasterq-dump"):
        raise RouteUnavailable("sra-tools (prefetch/fasterq-dump) not installed")
    d = cfg.tmpdir
    run(["prefetch", "--max-size", "u", "-O", d, acc], log_to=log)
    sra = os.path.join(d, acc, f"{acc}.sra")
    if not os.path.exists(sra):
        sra = os.path.join(d, f"{acc}.sra")
    if not nonempty(sra):
        raise RouteUnavailable(f"{acc}: prefetch produced no .sra")
    fq = _fasterq(sra, acc, cfg, log)
    if not cfg["keep.sra"]:
        rm(sra, os.path.join(d, acc))
    _finalise(fq, out_gz, threads=cfg["project.threads"])
    return {"url": f"sra:{acc}"}


# fasterq-dump is I/O-bound writing an uncompressed FASTQ, not CPU-bound: measured,
# -e1 8.6s -> -e4 6.7s -> -e16 6.6s. Threads past ~4-6 buy nothing and just take
# cores away from the other samples in the pool.
FASTERQ_THREADS = 6


def _fasterq(sra: str, acc: str, cfg: Config, log: str = "") -> str:
    """.sra -> a single FASTQ of the biological read.

    `--split-3` writes _1/_2 for a paired run and a bare file for a single-end
    one; ribo-seq is single-end, but a paired deposit (an RNA-seq control, or a
    UMI split into R2) must still resolve to the footprint read, which is R1.

    fasterq-dump needs roughly 10x the .sra size in scratch, so `-t` wants fast
    storage; /dev/shm is ideal when the run fits in it.
    """
    d = cfg.tmpdir
    scratch = _scratch(sra, d)
    threads = min(int(cfg["project.threads"]), FASTERQ_THREADS)
    run(["fasterq-dump", "-e", str(threads), "--split-3",
         "--skip-technical", "-t", scratch, "-O", d, sra], log_to=log)
    single = os.path.join(d, f"{acc}.fastq")
    r1 = os.path.join(d, f"{acc}_1.fastq")
    if nonempty(single):
        rm(r1, os.path.join(d, f"{acc}_2.fastq"))
        return single
    if nonempty(r1):
        rm(os.path.join(d, f"{acc}_2.fastq"), os.path.join(d, f"{acc}_3.fastq"))
        return r1
    raise RuntimeError(f"fasterq-dump produced no FASTQ for {acc}")


# ENA first because it is the fastest by a wide margin; the other two are the
# fallback for the runs ENA has not mirrored, tried in order. Not a preference --
# an availability chain. Do not reorder without re-reading docs/DOWNLOAD.md.
ROUTES = (
    ("ena_https", route_ena_https),
    ("aws_odp", route_aws_odp),
    ("prefetch", route_prefetch),
)


def download_full(acc: str, out_fastq_gz: str, cfg: Config, *, log: str = "") -> dict:
    """The whole run, through the first route that works. Idempotent."""
    if nonempty(out_fastq_gz):
        LOG.info("[%s] fastq already present", acc)
        return {"route": "cached", "bytes": os.path.getsize(out_fastq_gz),
                "seconds": 0.0, "mb_per_s": None, "attempts": 0}

    os.makedirs(cfg.tmpdir, exist_ok=True)

    errors = []
    for attempt, (name, fn) in enumerate(ROUTES, 1):
        t0 = time.time()
        try:
            info = _with_retries(fn, acc, out_fastq_gz, cfg, log)
        except RouteUnavailable as exc:
            LOG.info("[%s] route %s unavailable: %s", acc, name, exc)
            errors.append(f"{name}: {exc}")
            continue
        except Exception as exc:  # noqa: BLE001
            LOG.warning("[%s] route %s failed: %s", acc, name, exc)
            errors.append(f"{name}: {exc}")
            rm(out_fastq_gz)
            continue
        dt = time.time() - t0
        size = os.path.getsize(out_fastq_gz)
        return {"route": name, "bytes": size, "seconds": round(dt, 1),
                "mb_per_s": round(size / 1e6 / max(dt, 1e-9), 2),
                "attempts": attempt, "url": info.get("url", "")}

    raise RuntimeError(f"{acc}: every download route failed\n  " + "\n  ".join(errors))


def _with_retries(fn, acc: str, out: str, cfg: Config, log: str):
    """Transient network failures are the norm at this scale, not the exception."""
    tries = int(cfg["download.max_retries"])
    backoff = float(cfg["download.retry_backoff_s"])
    last: Exception | None = None
    for i in range(tries):
        try:
            return fn(acc, out, cfg, log)
        except RouteUnavailable:
            raise                       # a missing mirror will not fix itself
        except Exception as exc:        # noqa: BLE001
            last = exc
            if i == tries - 1:
                break
            wait = backoff * (2 ** i)
            LOG.info("[%s] retry %d/%d in %.0fs (%s)", acc, i + 1, tries, wait, exc)
            time.sleep(wait)
    raise last  # type: ignore[misc]
