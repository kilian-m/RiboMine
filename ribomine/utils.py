"""Shared helpers: logging, subprocess, atomic IO, TSV, and the sample layout.

Every stage writes into one per-sample directory and reads the previous stage's
JSON from it; `Sample` defines those paths in one place.
"""
from __future__ import annotations

import errno
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass

from xopen import xopen

LOG = logging.getLogger("ribomine")


# --- logging ---------------------------------------------------------------
def setup_logging(level: str = "INFO", logfile: str | None = None) -> None:
    fmt = "%(asctime)s %(levelname)-7s %(message)s"
    datefmt = "%H:%M:%S"
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    if logfile:
        os.makedirs(os.path.dirname(os.path.abspath(logfile)), exist_ok=True)
        handlers.append(logging.FileHandler(logfile))
    logging.basicConfig(level=getattr(logging, level.upper(), logging.INFO),
                        format=fmt, datefmt=datefmt, handlers=handlers, force=True)
    # matplotlib and urllib are chatty at DEBUG
    for noisy in ("matplotlib", "urllib3", "PIL"):
        logging.getLogger(noisy).setLevel(logging.WARNING)


# --- subprocess ------------------------------------------------------------
class ToolError(RuntimeError):
    """An external tool exited non-zero."""

    def __init__(self, cmd, returncode, stderr=""):
        self.cmd, self.returncode, self.stderr = cmd, returncode, stderr
        tail = "\n".join((stderr or "").strip().splitlines()[-12:])
        super().__init__(f"{cmd[0]} exited {returncode}\n  cmd: {' '.join(map(str, cmd))}"
                         + (f"\n  stderr:\n{tail}" if tail else ""))


def run(cmd, *, check=True, capture=True, stdout=None, cwd=None, env=None, log_to=None):
    """Run a command. Returns CompletedProcess; raises ToolError on failure.

    `log_to` appends the command's stderr to that (per-sample) log file.
    """
    cmd = [str(c) for c in cmd]
    LOG.debug("run: %s", " ".join(cmd))
    t0 = time.time()
    p = subprocess.run(
        cmd,
        stdout=stdout if stdout is not None else (subprocess.PIPE if capture else None),
        stderr=subprocess.PIPE if capture else None,
        text=True, cwd=cwd, env=env,
    )
    if log_to and p.stderr:
        with open(log_to, "a") as fh:
            fh.write(f"\n$ {' '.join(cmd)}\n{p.stderr}")
    LOG.debug("  -> %d in %.1fs", p.returncode, time.time() - t0)
    if check and p.returncode != 0:
        raise ToolError(cmd, p.returncode, p.stderr or "")
    return p


def require_tools(*tools: str) -> None:
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        raise RuntimeError(
            f"required tool(s) not on PATH: {', '.join(missing)}. "
            f"Activate the RiboMine environment: conda env create -f environment.yml && "
            f"conda activate ribomine"
        )


def have(tool: str) -> bool:
    return shutil.which(tool) is not None


# --- IO --------------------------------------------------------------------
def read_json(path: str, default=None):
    if not path or not os.path.exists(path):
        return default
    with open(path) as fh:
        return json.load(fh)


def write_json(path: str, obj) -> str:
    """Write JSON atomically, so a killed run cannot leave a half-written file
    that `resume` would treat as complete."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(os.path.abspath(path)), suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(obj, fh, indent=1, default=_jsonable)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    return path


def _jsonable(o):
    if hasattr(o, "item"):        # numpy scalars
        return o.item()
    if hasattr(o, "tolist"):      # numpy arrays
        return o.tolist()
    raise TypeError(f"not JSON-serialisable: {type(o)}")


def write_tsv(path: str, rows: list[dict], columns: list[str] | None = None) -> str:
    """Write a list of dicts as a TSV. Missing values render as empty cells."""
    import csv

    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    cols = columns or (list(rows[0].keys()) if rows else [])
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore",
                           lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: _cell(r.get(c)) for c in cols})
    os.replace(tmp, path)
    return path


def _cell(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return f"{v:.4g}"
    if isinstance(v, (list, tuple)):
        return ";".join(str(x) for x in v)
    if isinstance(v, dict):
        return ";".join(f"{k}={v}" for k, v in v.items())
    return str(v)


def read_tsv(path: str) -> list[dict]:
    import csv

    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_lines(path: str) -> list[str]:
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]


def open_fastq(path: str, mode: str = "rt", *, compresslevel: int = 6):
    """Open a FASTQ, transparently gzipped.

    Uses `xopen`, which hands (de)compression to an external pigz/igzip process
    instead of Python's slower in-process `gzip`. Writes default to level 6
    rather than gzip's 9: the files are intermediates, and level 9 costs 3-5x
    the CPU for a few percent of size.
    """
    if not path.endswith(".gz"):
        return open(path, mode)
    if "w" in mode or "a" in mode:
        return xopen(path, mode, compresslevel=compresslevel)
    return xopen(path, mode)


def count_fastq_reads(path: str) -> int:
    n = 0
    with open_fastq(path, "rb") as fh:
        for i, _ in enumerate(fh):
            n = i
    return (n + 1) // 4


def human(n) -> str:
    if n is None:
        return "n/a"
    for unit in ("", "k", "M", "G", "T"):
        if abs(n) < 1000:
            return f"{n:.0f}{unit}" if unit else f"{n:.0f}"
        n /= 1000.0
    return f"{n:.1f}P"


def rm(*paths: str) -> None:
    for p in paths:
        try:
            if os.path.isdir(p) and not os.path.islink(p):
                shutil.rmtree(p)
            elif os.path.exists(p) or os.path.islink(p):
                os.remove(p)
        except OSError as exc:
            if exc.errno != errno.ENOENT:
                LOG.warning("could not remove %s: %s", p, exc)


def nonempty(path: str | None) -> bool:
    return bool(path) and os.path.exists(path) and os.path.getsize(path) > 0


# --- the per-sample layout -------------------------------------------------
@dataclass
class Sample:
    """The paths of every artefact of one run, from the accession and the workdir.

    Stages and `resume` use these properties instead of building paths themselves.
    """

    acc: str
    workdir: str

    # -- directories
    @property
    def dir(self) -> str:
        return _mk(os.path.join(self.workdir, "samples", self.acc))

    @property
    def star_local(self) -> str:
        return _mk(os.path.join(self.dir, "star_local"))

    @property
    def star_final(self) -> str:
        return _mk(os.path.join(self.dir, "star_final"))

    # -- stage 2 (QC) artefacts
    @property
    def sample_fastq(self) -> str:
        return os.path.join(self.dir, f"{self.acc}.sample.fastq")

    @property
    def filtered_fastq(self) -> str:
        return os.path.join(self.dir, f"{self.acc}.filtered.fastq")

    @property
    def local_bam(self) -> str:
        return os.path.join(self.star_local, "aligned.bam")

    @property
    def pileup_bam(self) -> str:
        return os.path.join(self.star_local, "aligned.pileup.bam")

    @property
    def local_log(self) -> str:
        return os.path.join(self.star_local, "Log.final.out")

    @property
    def contam_json(self) -> str:
        return os.path.join(self.dir, f"{self.acc}.contam.json")

    @property
    def pileup_json(self) -> str:
        return os.path.join(self.dir, f"{self.acc}.pileup.json")

    @property
    def profile_json(self) -> str:
        return os.path.join(self.dir, f"{self.acc}.profile.json")

    @property
    def qc_json(self) -> str:
        return os.path.join(self.dir, f"{self.acc}.qc.json")

    # -- stage 3 (architecture)
    @property
    def arch_json(self) -> str:
        return os.path.join(self.dir, f"{self.acc}.architecture.json")

    # -- stage 4 (full processing)
    @property
    def full_fastq(self) -> str:
        return os.path.join(_mk(os.path.join(self.workdir, "fastq")), f"{self.acc}.fastq.gz")

    @property
    def trimmed_fastq(self) -> str:
        return os.path.join(self.dir, f"{self.acc}.trimmed.fastq.gz")

    @property
    def clean_fastq(self) -> str:
        return os.path.join(self.dir, f"{self.acc}.clean.fastq.gz")

    @property
    def trim_json(self) -> str:
        return os.path.join(self.dir, f"{self.acc}.trim.json")

    @property
    def process_json(self) -> str:
        return os.path.join(self.dir, f"{self.acc}.process.json")

    @property
    def bam(self) -> str:
        return os.path.join(_mk(os.path.join(self.workdir, "bams")), f"{self.acc}.bam")

    @property
    def log(self) -> str:
        return os.path.join(_mk(os.path.join(self.workdir, "logs")), f"{self.acc}.log")

    # -- plots
    def qc_plot(self, ext: str = "png") -> str:
        return os.path.join(_mk(os.path.join(self.workdir, "qc", "plots")),
                            f"{self.acc}.qc.{ext}")

    def arch_plot(self, ext: str = "png") -> str:
        return os.path.join(_mk(os.path.join(self.workdir, "architecture", "plots")),
                            f"{self.acc}.arch.{ext}")


def _mk(p: str) -> str:
    os.makedirs(p, exist_ok=True)
    return p
