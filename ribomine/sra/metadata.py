"""ENA portal metadata: the run records everything downstream is keyed on.

ENA and SRA mirror each other's submissions, but ENA's portal API is the one
worth talking to: a single request returns the full run record *including the
direct fastq URLs*, where NCBI would need an esearch, an esummary and then a
separate resolver call to find out where the data actually lives.

Caveat that bites: DRR runs submitted to DDBJ have a browser record but no
portal record, so they come back empty. Those fall back to the SRA toolkit
download route.
"""
from __future__ import annotations

import io
import time
import urllib.error
import urllib.parse
import urllib.request

from ..utils import LOG

PORTAL = "https://www.ebi.ac.uk/ena/portal/api"

# the fields every stage may ask for, in the order they appear in candidates.tsv
FIELDS = [
    "run_accession",
    "experiment_accession",
    "sample_accession",
    "study_accession",
    "secondary_study_accession",
    "scientific_name",
    "tax_id",
    "instrument_platform",
    "instrument_model",
    "library_strategy",
    "library_selection",
    "library_source",
    "library_layout",
    "library_construction_protocol",
    "read_count",
    "base_count",
    "fastq_ftp",
    "fastq_bytes",
    "fastq_md5",
    "sra_ftp",
    "submitted_ftp",
    "study_title",
    "experiment_title",
    "sample_title",
    "description",
    "first_public",
    "center_name",
]


def _post(endpoint: str, params: dict, *, retries: int = 4, timeout: int = 180) -> str:
    data = urllib.parse.urlencode(params).encode()
    url = f"{PORTAL}/{endpoint}"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, data=data, timeout=timeout) as fh:
                return fh.read().decode("utf-8", "replace")
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            if attempt == retries - 1:
                raise
            wait = 5 * (attempt + 1)
            LOG.warning("ENA %s failed (%s); retrying in %ds", endpoint, exc, wait)
            time.sleep(wait)
    return ""


def _parse_tsv(text: str) -> list[dict]:
    import csv

    text = (text or "").strip("\n")
    if not text:
        return []
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    return [r for r in rows if r.get("run_accession")]


def fetch_runs(accessions: list[str], *, batch: int = 100,
               fields: list[str] | None = None) -> list[dict]:
    """ENA read_run records for a list of run accessions.

    Runs ENA does not know about (DDBJ-only DRR runs, very fresh submissions) are
    simply absent from the result -- the caller decides what to do about that.
    """
    out: list[dict] = []
    flds = ",".join(fields or FIELDS)
    for i in range(0, len(accessions), batch):
        chunk = accessions[i: i + batch]
        txt = _post("search", {
            "result": "read_run",
            "includeAccessions": ",".join(chunk),
            "fields": flds,
            "format": "tsv",
            "limit": "0",
        })
        out += _parse_tsv(txt)
        if len(accessions) > batch:
            LOG.info("  ENA metadata %d/%d -> %d rows", min(i + batch, len(accessions)),
                     len(accessions), len(out))
            time.sleep(0.3)
    return out


def search(query: str, *, fields: list[str] | None = None, limit: int = 0) -> list[dict]:
    """ENA portal query-language search over read_run."""
    txt = _post("search", {
        "result": "read_run",
        "query": query,
        "fields": ",".join(fields or FIELDS),
        "format": "tsv",
        "limit": str(limit),
    })
    return _parse_tsv(txt)


def filereport(acc: str, fields: list[str]) -> dict:
    """The lightweight per-run endpoint (used for fastq URLs / read counts)."""
    q = urllib.parse.urlencode({
        "accession": acc, "result": "read_run",
        "fields": ",".join(fields), "format": "tsv",
    })
    url = f"{PORTAL}/filereport?{q}"
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                rows = _parse_tsv(fh.read().decode())
            return rows[0] if rows else {}
        except (urllib.error.URLError, TimeoutError, OSError):
            if attempt == 3:
                return {}
            time.sleep(3 * (attempt + 1))
    return {}


def fastq_urls(acc: str) -> list[str]:
    """Direct https FASTQ URLs for a run, R1 first.

    ENA serves the submitter's own FASTQ, so nothing has to be converted out of
    an .sra container. A run with no ENA fastq mirror returns [] and the caller
    falls back to another route.
    """
    row = filereport(acc, ["fastq_ftp", "fastq_bytes"])
    paths = [p for p in (row.get("fastq_ftp") or "").split(";") if p]
    if not paths:
        return []
    # A paired deposit lists _1 and _2. Ribo-seq is single-end; when a run is
    # paired anyway (an RNA-seq control in the same study, or a UMI read split out
    # into R2), the footprint is in R1.
    if len(paths) > 1:
        r1 = [p for p in paths if p.endswith("_1.fastq.gz")]
        paths = r1 or paths[:1]
    return ["https://" + p for p in paths]


def read_count(acc: str) -> int | None:
    row = filereport(acc, ["read_count"])
    try:
        return int(row["read_count"])
    except (KeyError, TypeError, ValueError):
        return None
