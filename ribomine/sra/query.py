"""Find candidate Ribo-seq runs in the SRA/ENA.

Ribo-seq has no `library_strategy` of its own (it is deposited as `RNA-Seq`,
`OTHER` or `ncRNA-Seq`), so runs are found by their free text, in three steps:

1. Recall. ENA's portal index is token-based and a multi-word wildcard matches
   nothing, so single tokens (`query.terms`) are ORed across the text fields.
2. Precision. Phrases are matched here as regexes: STRONG (unambiguous, e.g.
   "ribosome profiling") and WEAK (also used in plain RNA-seq, e.g.
   "translatome"). EXCLUDE (rRNA-depleted RNA-seq, e.g. "Ribo-Zero") vetoes
   weak hits only, because Ribo-seq protocols deplete rRNA too.
3. Verification is left to the QC stage, so the query is tuned for recall.

`query.source: "entrez" | "both"` adds NCBI's index, which also holds SRA-only
and very recent submissions.
"""
from __future__ import annotations

import io
import re
import time
import urllib.parse
import urllib.request

from ..config import Config
from ..utils import LOG, write_tsv
from . import metadata

# --- tier 1: unambiguous; never vetoed by EXCLUDE ---
STRONG = re.compile(r"""
    (?<![a-z]) ribo [\s._-]* seq (?![a-z])
  | ribosom\w* [\s._-]* (?: profiling | profile | footprint\w* | protected )
  | ribosome [\s._-]* protected [\s._-]* (?: fragment | footprint )\w*
  | \b RPFs? \b
  | (?<![a-z]) ribo [\s._-]* profil\w*
  | \b [QG]? TI [\s._-]* seq \b
  | \b ARTseq \b | \b RiboLace \b
  | monosome [\s._-]* (?: seq | profil\w* )
  | translatome [\s._-]* seq
""", re.I | re.X)

# --- tier 2: suggestive, but also used by plain RNA-seq studies ---
WEAK = re.compile(r"""
    translatom\w*
  | harringtonine | lactimidomycin
  | \b polysome [\s._-]* profil\w*
  | \b 40S [\s._-]* footprint\w*
""", re.I | re.X)

# --- rRNA-depleted RNA-seq: vetoes a WEAK hit, never a STRONG one ---
EXCLUDE = re.compile(r"""
    ribo [\s_-]* zero | ribozero
  | ribosomal \s+ RNA \s+ deplet\w*
  | rRNA [\s_-]* deplet\w*
  | ribo [\s_-]* minus | ribodeplet\w*
  | riboswitch | ribozyme
""", re.I | re.X)

# fields the tiers are matched against
TEXT_FIELDS = ("study_title", "experiment_title", "sample_title",
               "library_construction_protocol", "description")

# fields ENA accepts a wildcard token search on
SEARCHABLE = ("study_title", "experiment_title", "sample_title",
              "library_construction_protocol")


def _blob(row: dict) -> str:
    return " | ".join((row.get(f) or "") for f in TEXT_FIELDS)


def classify(row: dict) -> str:
    """'strong' | 'weak' | 'excluded' | 'none' -- why this run was kept or dropped."""
    b = _blob(row)
    if STRONG.search(b):
        return "strong"
    if WEAK.search(b):
        return "excluded" if EXCLUDE.search(b) else "weak"
    return "none"


# --- route 1: ENA portal ---
def _ena_query(cfg: Config) -> str:
    """OR the wildcard tokens (single words only) across the searchable text fields."""
    tokens = [t.strip() for t in cfg["query.terms"] if t.strip()]
    if not tokens:
        raise ValueError("query.terms is empty")
    if any(" " in t for t in tokens):
        raise ValueError(
            "query.terms must be single words -- ENA's text index is tokenised and a "
            "multi-word wildcard matches nothing. Offending: "
            + ", ".join(repr(t) for t in tokens if " " in t))
    clauses = [f'{f}="*{t}*"' for f in SEARCHABLE for t in tokens]
    taxon = cfg.get("query.taxon_id") or cfg.get("reference.taxon_id")
    parts = [f"tax_eq({int(taxon)})"]
    plat = cfg["query.instrument_platform"]
    if plat and plat.lower() != "any":
        parts.append(f'instrument_platform="{plat}"')
    parts.append("(" + " OR ".join(clauses) + ")")
    return " AND ".join(parts)


def search_ena(cfg: Config) -> list[dict]:
    q = _ena_query(cfg)
    LOG.info("ENA portal search (%d token(s) x %d field(s))",
             len(cfg["query.terms"]), len(SEARCHABLE))
    LOG.debug("  query: %s", q)
    t0 = time.time()
    rows = metadata.search(q)
    LOG.info("  ENA returned %d candidate run(s) in %.0fs", len(rows), time.time() - t0)
    for r in rows:
        r["source"] = "ena"
    return rows


# --- route 2: NCBI Entrez ---
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _entrez_term(cfg: Config) -> str:
    """Entrez supports phrase search, so it gets the precise phrases directly."""
    org = cfg["reference.species"]
    phrases = [
        "ribo-seq", "riboseq", "ribosome profiling", "ribosome footprint",
        "ribosome footprinting", "ribosome protected fragments",
        "ribosome-protected fragments", "ribosome profiling sequencing",
        "RPF", "translatome", "TI-seq", "QTI-seq", "GTI-seq", "ARTseq",
        "harringtonine", "lactimidomycin",
    ]
    ors = " OR ".join(f'"{p}"[All Fields]' for p in phrases)
    return f'"{org}"[Organism] AND ({ors})'


def _eget(endpoint: str, params: dict, *, retries: int = 4) -> str:
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=180) as fh:
                return fh.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            if attempt == retries - 1:
                raise
            LOG.warning("Entrez %s failed (%s); retrying", endpoint, exc)
            time.sleep(4 * (attempt + 1))
    return ""


def search_entrez(cfg: Config) -> list[dict]:
    """esearch (history) -> efetch runinfo CSV, mapped onto the ENA field names."""
    import csv

    key = cfg["query.ncbi_api_key"]
    base = {"db": "sra", "retmode": "json"}
    if key:
        base["api_key"] = key
    if cfg["query.ncbi_email"]:
        base["email"] = cfg["query.ncbi_email"]

    LOG.info("NCBI Entrez search (db=sra)")
    js = _eget("esearch.fcgi", {**base, "term": _entrez_term(cfg),
                                "usehistory": "y", "retmax": "0"})
    import json as _json
    res = _json.loads(js)["esearchresult"]
    total = int(res["count"])
    webenv, qkey = res["webenv"], res["querykey"]
    LOG.info("  Entrez matched %d SRA record(s)", total)

    rows: list[dict] = []
    step = 500
    pause = 0.11 if key else 0.35        # 10/s with a key, 3/s without
    for start in range(0, total, step):
        csv_txt = _eget("efetch.fcgi", {
            "db": "sra", "WebEnv": webenv, "query_key": qkey,
            "rettype": "runinfo", "retmode": "text",
            "retstart": str(start), "retmax": str(step),
            **({"api_key": key} if key else {}),
        })
        for r in csv.DictReader(io.StringIO(csv_txt)):
            if not r.get("Run"):
                continue
            rows.append({
                "run_accession": r["Run"],
                "experiment_accession": r.get("Experiment", ""),
                "sample_accession": r.get("Sample", ""),
                "study_accession": r.get("BioProject", ""),
                "secondary_study_accession": r.get("SRAStudy", ""),
                "scientific_name": r.get("ScientificName", ""),
                "tax_id": r.get("TaxID", ""),
                "instrument_platform": r.get("Platform", ""),
                "instrument_model": r.get("Model", ""),
                "library_strategy": r.get("LibraryStrategy", ""),
                "library_selection": r.get("LibrarySelection", ""),
                "library_source": r.get("LibrarySource", ""),
                "library_layout": (r.get("LibraryLayout") or "").upper(),
                "read_count": r.get("spots", ""),
                "base_count": r.get("bases", ""),
                "experiment_title": r.get("Experiment_Title", "") or r.get("LibraryName", ""),
                "sample_title": r.get("SampleName", ""),
                "study_title": r.get("Study_Pubmed_id", "") and "" or "",
                "first_public": r.get("ReleaseDate", ""),
                "center_name": r.get("CenterName", ""),
                "source": "entrez",
            })
        LOG.info("  Entrez runinfo %d/%d -> %d run(s)", min(start + step, total), total, len(rows))
        time.sleep(pause)
    return rows


# --- union, filter, write ---
def search(cfg: Config) -> list[dict]:
    src = cfg["query.source"]
    rows: list[dict] = []
    if src in ("ena", "both"):
        rows += search_ena(cfg)
    if src in ("entrez", "both"):
        try:
            rows += search_entrez(cfg)
        except Exception as exc:  # noqa: BLE001
            if src == "entrez":
                raise
            LOG.warning("Entrez search failed (%s); continuing with the ENA hits alone", exc)

    # union on run_accession; the first record (ENA's) wins
    merged: dict[str, dict] = {}
    for r in rows:
        acc = r["run_accession"]
        if acc in merged:
            if merged[acc].get("source") != r.get("source"):
                merged[acc]["source"] = "ena+entrez"
            for k, v in r.items():        # fill blanks from the other archive
                if v and not merged[acc].get(k):
                    merged[acc][k] = v
        else:
            merged[acc] = dict(r)

    # Entrez-only runs lack fastq_ftp / study_title: ask ENA for their record.
    # Runs ENA does not know keep their Entrez fields.
    thin = [a for a, r in merged.items() if not r.get("fastq_ftp")]
    if thin:
        LOG.info("filling ENA metadata for %d Entrez-only run(s)", len(thin))
        for r in metadata.fetch_runs(thin):
            acc = r["run_accession"]
            merged[acc] = {**merged[acc], **{k: v for k, v in r.items() if v}}

    return list(merged.values())


def _filter(cfg: Config, rows: list[dict]) -> list[dict]:
    """The tiered text filter, then the eligibility filters and the run caps."""
    stats = {"strong": 0, "weak": 0, "excluded": 0, "none": 0}
    kept = []
    for r in rows:
        tier = classify(r)
        stats[tier] += 1
        if tier in ("strong", "weak"):
            r["match_tier"] = tier
            kept.append(r)
    LOG.info("text filter: %d strong, %d weak, %d weak-but-excluded, %d no ribo-seq phrase",
             stats["strong"], stats["weak"], stats["excluded"], stats["none"])

    strategies = set(cfg["query.library_strategies"] or [])
    layout = (cfg["query.layout"] or "any").upper()
    min_reads = int(cfg["query.min_read_count"] or 0)
    ex_runs = set(cfg["query.exclude_runs"] or [])
    ex_studies = set(cfg["query.exclude_studies"] or [])

    out = []
    dropped = {"strategy": 0, "layout": 0, "depth": 0, "excluded": 0}
    for r in kept:
        if r["run_accession"] in ex_runs or r.get("study_accession") in ex_studies:
            dropped["excluded"] += 1
            continue
        if strategies and r.get("library_strategy") not in strategies:
            dropped["strategy"] += 1
            continue
        if layout != "ANY" and (r.get("library_layout") or "").upper() != layout:
            dropped["layout"] += 1
            continue
        try:
            n = int(r.get("read_count") or 0)
        except ValueError:
            n = 0
        if min_reads and n < min_reads:
            dropped["depth"] += 1
            continue
        out.append(r)
    LOG.info("eligibility: dropped %d on library_strategy, %d on layout, %d below %s reads, "
             "%d explicitly excluded",
             dropped["strategy"], dropped["layout"], dropped["depth"],
             f"{min_reads:,}", dropped["excluded"])

    # at most N runs per study, deepest first (replicates share a protocol)
    per_study = int(cfg["query.runs_per_study"] or 0)
    if per_study:
        by_study: dict[str, list[dict]] = {}
        for r in out:
            by_study.setdefault(r.get("study_accession") or r["run_accession"], []).append(r)
        trimmed = []
        for runs in by_study.values():
            runs.sort(key=lambda r: (-int(r.get("read_count") or 0), r["run_accession"]))
            trimmed += runs[:per_study]
        LOG.info("runs_per_study=%d: %d -> %d run(s)", per_study, len(out), len(trimmed))
        out = trimmed

    out.sort(key=lambda r: (r.get("study_accession") or "", r["run_accession"]))
    cap = int(cfg["query.max_runs"] or 0)
    if cap and len(out) > cap:
        LOG.info("max_runs=%d: capping %d -> %d run(s)", cap, len(out), cap)
        out = out[:cap]
    return out


COLUMNS = ["run_accession", "study_accession", "experiment_accession", "sample_accession",
           "scientific_name", "library_strategy", "library_selection", "library_layout",
           "instrument_model", "read_count", "base_count", "fastq_bytes",
           "match_tier", "source", "first_public", "center_name",
           "study_title", "experiment_title", "sample_title",
           "library_construction_protocol", "fastq_ftp"]


def run_query(cfg: Config) -> str:
    rows = search(cfg)
    LOG.info("archives returned %d unique run(s)", len(rows))
    kept = _filter(cfg, rows)
    studies = len({r.get("study_accession") for r in kept})
    LOG.info("=> %d candidate ribo-seq run(s) across %d stud(y/ies)", len(kept), studies)

    import os
    path = os.path.join(cfg.dir("meta"), "candidates.tsv")
    write_tsv(path, kept, COLUMNS)
    return path
