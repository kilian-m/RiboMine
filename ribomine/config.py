"""Configuration: defaults, loading, validation, and derived paths.

One JSON file describes a run. `DEFAULTS` is its schema and its documentation:
every key a user may set, with the value used when it is absent. A user's config
is deep-merged onto it, and an unknown key is an error.
"""
from __future__ import annotations

import copy
import json
import os
from dataclasses import asdict, dataclass
from typing import Any

from fqdissect.infer import Thresholds

DEFAULTS: dict[str, Any] = {
    "project": {
        "name": "ribomine",
        "workdir": "work",     # relative paths resolve against the config file's directory
        "threads": 8,          # threads per sample (STAR, bowtie2, cutadapt)
        "jobs": 8,             # samples processed concurrently
        "seed": 20260712,      # read sampling
    },

    "pipeline": {
        "start": "query",         # query | accessions | fastq
        "end": "bam",             # qc | architecture | bam
        "accession_list": None,   # start == "accessions": one run accession per line
        "fastq_dir": None,        # start == "fastq": directory of *.fastq[.gz]
        "resume": True,           # skip any sample whose stage output already exists
        # QC verdicts carried into the later stages (TI-seq is a Ribo-seq variant)
        "keep_verdicts": ["RIBO-SEQ", "TI-SEQ"],
    },

    "reference": {
        "species": "Homo sapiens",
        "taxon_id": 9606,
        "genome_fasta": None,        # required
        "gtf": None,                 # required (Ensembl-style, with CDS frames)
        "star_index": None,          # required; build it with the GTF to get gene counts
        "contaminant_fasta": None,   # null = the bundled human rRNA/tRNA/snRNA/snoRNA/Mt set
        "contaminant_index": None,   # bowtie2 index prefix; null = built on first use
        "annotation_index": None,    # cached GTF index; null = built on first use
    },

    "query": {
        "source": "both",            # ena | entrez | both (neither is a superset)
        # Single words only: ENA's text index is tokenised. These are the broad recall
        # net; `sra/query.py` applies the phrase filters and QC has the final word.
        "terms": ["ribo", "ribosome", "riboseq", "RPF", "translatome",
                  "footprint", "harringtonine", "lactimidomycin", "monosome",
                  "ARTseq", "RiboLace"],
        "taxon_id": None,            # null = reference.taxon_id
        # Ribo-seq has no library_strategy of its own; submitters use these
        "library_strategies": ["RNA-Seq", "OTHER", "ncRNA-Seq", "miRNA-Seq"],
        "instrument_platform": "ILLUMINA",   # "any" to drop the filter
        "layout": "SINGLE",          # SINGLE | PAIRED | any
        "min_read_count": 1_000_000,
        "max_runs": 0,               # 0 = no cap
        "runs_per_study": 0,         # 0 = all; N = the N deepest runs per study
        "exclude_runs": [],
        "exclude_studies": [],
        "ncbi_api_key": None,        # raises the Entrez rate limit from 3/s to 10/s
        "ncbi_email": None,
    },

    # Full runs come from ENA over HTTPS, falling back to the SRA mirrors for runs
    # ENA does not have (docs/DOWNLOAD.md).
    "download": {
        "connections": 8,        # parallel connections per file
        "max_retries": 4,
        "retry_backoff_s": 5,
        "tmpdir": None,          # null = <workdir>/tmp
    },

    # What survives the run. The JSONs, TSVs and plots are always kept.
    "keep": {
        "bam": True,             # trimmed, filtered, mapped reads
        "fastq": False,          # the raw run FASTQ
        "trimmed_fastq": False,
        "clean_fastq": False,    # after contaminant removal: the reads that were mapped
        "sra": False,            # the .sra container (SRA fallback routes only)
        "qc_fastq": False,       # the QC read sample
        "qc_bam": False,         # its local alignments (the evidence for verdict and call)
    },

    "qc": {
        "sample_reads": 200_000,     # reads kept per run for QC/architecture
        "scan_reads": 1_000_000,     # reads streamed to sample from (0 = whole run)
        "max_reads_scored": 200_000,
        # verdict thresholds, tuned on a 187-sample human cohort
        "footprint_len_lo": 25,
        "footprint_len_hi": 36,
        "read_len_peak_frac_min": 0.40,
        "cds_enrich_min": 0.20,      # fraction of genic reads in CDS to look like Ribo-seq
        "cds_strong": 0.55,          # strong CDS enrichment excuses weak periodicity
        # hard gates, whatever the periodicity
        "cds_region_min": 0.50,      # CDS fraction of all reads
        "read_len_min_frac": 0.75,   # fraction of reads in the footprint length window
        "periodic_min": 0.42,        # in-frame fraction (chance = 1/3)
        "periodic_strong": 0.50,
        "tvd_min": 0.10,             # TVD of the frame distribution to uniform
        "tvd_strong": 0.20,
        "min_cds_reads": 150,
        "single_locus_max": 0.50,    # more than this share at one 5' locus = contaminant
        "locus_vs_cds": 1.0,         # one locus holding as many reads as the whole CDS
        "locus_min_to_judge": 0.05,
        "min_unique_frac": 0.15,     # below this the library maps too poorly to use
        # start-codon peak / CDS body above this = TI-seq (elongating Ribo-seq tops
        # out around 30)
        "tiseq_ratio_min": 30,
        "tiseq_min_peak": 200,
    },

    "contaminants": {
        "enabled": True,
        "min_entropy": 1.1,          # low-complexity read filter (bits)
        "max_base_frac": 0.85,
        # position pile-up filter: a 5' position is removed when it is over-represented
        # AND its reads are concentrated on one length (a footprint pile spreads over
        # lengths; an adapter dimer or a miRNA does not) ...
        "pileup_min_count": 30,
        "pileup_min_frac": 0.005,
        "pileup_len_conc": 0.75,
        # ... or when it holds this share of the reads on its own
        "pileup_max_frac": 0.10,
    },

    "architecture": {
        # fqdissect's calling thresholds (documented in fqdissect.infer.Thresholds)
        **asdict(Thresholds()),
        # carry runs whose architecture is `undetermined` into the bam stage with a
        # best-effort trim?
        "process_undetermined": False,
    },

    "process": {
        "min_len": 20,               # drop reads shorter than this after trimming
        # Drop reads that do not show the 3' adapter: their 3' end is the read's end,
        # not the footprint's, and a 3' UMI was never sequenced.
        "discard_untrimmed": True,
        # nt of the 3' scaffold (barcode + adapter) that must be visible to cut on it.
        # Lower it to recover depth when the scaffold runs off most reads, at the cost
        # of specificity.
        "adapter_min_overlap": 7,
        "umi_dedup": False,
        "umi_dedup_tool": "umicollapse",     # umicollapse | umi_tools
        # directional | adjacency | cluster, or unique | percentile (umi_tools only)
        "umi_dedup_method": "directional",
        "umi_dedup_mem_gb": 8,       # JVM heap for umicollapse
        "filter_contaminants": True, # bowtie2 rRNA/tRNA removal before mapping
        "filter_pileups": True,      # pile-up removal after mapping
    },

    # The alignment of the trimmed reads: the output BAM. (The QC alignment's
    # settings are fqdissect's and are not configurable.)
    "mapping": {
        "align_ends_type": "Local",
        "multimap_nmax": 10,         # 1 = unique only
        "mismatch_nmax": 3,
        "mismatch_noverlmax": 0.1,
        "match_nmin_over_lread": 0.9,
        "align_intron_max": 0,       # 0 = STAR default (spliced); 1 = unspliced only
        "extra_args": [],
    },

    "plots": {
        "enabled": True,
        "format": "png",             # png | pdf | svg
        "dpi": 110,
    },
}


def deep_merge(base: dict, over: dict) -> dict:
    """Recursively overlay `over` on a copy of `base`."""
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = v
    return out


class ConfigError(Exception):
    pass


@dataclass
class Config:
    """A validated config plus the derived paths every stage agrees on."""

    data: dict
    path: str = ""

    # -- dotted access ------------------------------------------------------
    def get(self, dotted: str, default=None):
        cur: Any = self.data
        for part in dotted.split("."):
            if not isinstance(cur, dict) or part not in cur:
                return default
            cur = cur[part]
        return cur

    def __getitem__(self, dotted: str):
        sentinel = object()
        v = self.get(dotted, sentinel)
        if v is sentinel:
            raise KeyError(dotted)
        return v

    # -- derived paths ------------------------------------------------------
    @property
    def workdir(self) -> str:
        return self._abs(self.get("project.workdir", "work"))

    def _abs(self, p: str | None) -> str | None:
        """Resolve a path relative to the config file's directory."""
        if p is None:
            return None
        if os.path.isabs(p):
            return p
        base = os.path.dirname(os.path.abspath(self.path)) if self.path else os.getcwd()
        return os.path.normpath(os.path.join(base, p))

    def dir(self, *parts: str) -> str:
        """A workdir subdirectory, created on demand."""
        p = os.path.join(self.workdir, *parts)
        os.makedirs(p, exist_ok=True)
        return p

    @property
    def tmpdir(self) -> str:
        t = self.get("download.tmpdir")
        return self._abs(t) if t else self.dir("tmp")

    @property
    def annotation_index(self) -> str:
        p = self.get("reference.annotation_index")
        if p:
            return self._abs(p)
        stem = os.path.splitext(os.path.basename(self.get("reference.gtf") or "annotation.gtf"))[0]
        return os.path.join(self.dir("refs"), f"{stem}.idx.pkl")

    def ref(self, key: str) -> str | None:
        return self._abs(self.get(f"reference.{key}"))

    # -- validation ---------------------------------------------------------
    def validate(self, *, need_reference: bool = True, need_inputs: bool = True) -> None:
        """Check the config makes sense for what is about to run. `need_inputs` is
        False for sub-commands that do not consume the start point (`ribomine query`)."""
        start = self.get("pipeline.start")
        end = self.get("pipeline.end")
        if start not in START_POINTS:
            raise ConfigError(f"pipeline.start must be one of {START_POINTS}, got {start!r}")
        if end not in END_POINTS:
            raise ConfigError(f"pipeline.end must be one of {END_POINTS}, got {end!r}")
        if need_inputs:
            if start == "accessions":
                p = self.get("pipeline.accession_list")
                if not p:
                    raise ConfigError("pipeline.start='accessions' needs pipeline.accession_list")
                if not os.path.isfile(self._abs(p)):
                    raise ConfigError(f"accession list not found: {self._abs(p)}")
            if start == "fastq":
                p = self.get("pipeline.fastq_dir")
                if not p:
                    raise ConfigError("pipeline.start='fastq' needs pipeline.fastq_dir")
                if not os.path.isdir(self._abs(p)):
                    raise ConfigError(f"fastq dir not found: {self._abs(p)}")
        if not need_reference:
            return
        for key in ("genome_fasta", "gtf", "star_index"):
            p = self.ref(key)
            if not p:
                raise ConfigError(f"reference.{key} is required")
            if not os.path.exists(p):
                raise ConfigError(f"reference.{key} does not exist: {p}")


START_POINTS = ("query", "accessions", "fastq")
END_POINTS = ("qc", "architecture", "bam")
# stage order; a run stops after the stage named by pipeline.end
STAGES = ("query", "qc", "architecture", "bam")


def stages_to_run(start: str, end: str) -> list[str]:
    """The stages a run executes, given its start and end point."""
    first = 0 if start == "query" else 1          # a list/dir skips the query stage
    return list(STAGES[first : STAGES.index(end) + 1])


def load(path: str | None, overrides: dict | None = None) -> Config:
    """Load a config file, merge it onto the defaults, apply CLI overrides."""
    user: dict = {}
    if path:
        if not os.path.isfile(path):
            raise ConfigError(f"config not found: {path}")
        with open(path) as fh:
            try:
                user = json.load(fh)
            except json.JSONDecodeError as exc:
                raise ConfigError(f"{path}: invalid JSON: {exc}") from exc
    unknown = _unknown_keys(user, DEFAULTS)
    if unknown:
        raise ConfigError(
            "unknown config key(s): " + ", ".join(sorted(unknown))
        )
    data = deep_merge(DEFAULTS, user)
    data = deep_merge(data, overrides or {})
    return Config(data=data, path=path or "")


def _unknown_keys(user: dict, schema: dict, prefix: str = "") -> list[str]:
    """Keys in `user` that the schema does not define. Keys starting with `_` are
    comments and are ignored."""
    bad = []
    for k, v in (user or {}).items():
        if k.startswith("_"):
            continue
        if k not in schema:
            bad.append(f"{prefix}{k}")
        elif isinstance(v, dict) and isinstance(schema[k], dict):
            bad += _unknown_keys(v, schema[k], f"{prefix}{k}.")
    return bad


# What `null` means, per key: mostly "work it out for me", not "off".
# `write_example` puts this into the generated config.
NULL_MEANS = {
    "reference.genome_fasta": "REQUIRED: the genome FASTA.",
    "reference.gtf": "REQUIRED: the annotation GTF.",
    "reference.star_index": "REQUIRED: a STAR index of the genome, built with the GTF.",
    "reference.contaminant_fasta":
        "use the human contaminant reference bundled with RiboMine. Contaminant "
        "filtering is on either way; set contaminants.enabled=false to turn it off.",
    "reference.contaminant_index":
        "build the bowtie2 index from contaminant_fasta into <workdir>/refs/ on first use.",
    "reference.annotation_index":
        "build the cached GTF index into <workdir>/refs/<gtf-stem>.idx.pkl on first use.",
    "query.taxon_id": "use reference.taxon_id.",
    "query.ncbi_api_key": "no key: NCBI Entrez is rate-limited to 3 requests/s instead of 10.",
    "query.ncbi_email": "not sent to NCBI (optional, they ask for it on heavy use).",
    "download.tmpdir": "use <workdir>/tmp.",
    "pipeline.accession_list": "not used (only needed when pipeline.start = 'accessions').",
    "pipeline.fastq_dir": "not used (only needed when pipeline.start = 'fastq').",
}


def write_example(path: str) -> None:
    """Write a config with every default, and with what `null` means."""
    doc = {
        "_README": [
            "Every key below is RiboMine's default; delete what you do not change. "
            "Keys starting with '_' are comments. An unknown key is an error.",
        ],
        "_null_means": NULL_MEANS,
        **DEFAULTS,
    }
    with open(path, "w") as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")
