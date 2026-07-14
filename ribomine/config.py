"""Configuration: defaults, loading, validation, and derived paths.

A RiboMine run is fully described by one JSON file. `DEFAULTS` below is that
file's schema *and* its documentation -- every key a user may set appears here
with the value the pipeline uses when the key is absent. A user's config is
deep-merged onto it, so a config that sets only `reference.genome_fasta` is
legal and everything else falls back.

Thresholds live here rather than as module constants because they are the
knobs a user tunes per organism / per cohort. The values are the ones tuned on
a 187-sample human cohort; changing them changes calls, so they are surfaced
rather than buried.
"""
from __future__ import annotations

import copy
import json
import os
from dataclasses import dataclass
from typing import Any

# --- the schema, with every default ----------------------------------------
DEFAULTS: dict[str, Any] = {
    "project": {
        "name": "ribomine",
        # everything the run writes lands under here (relative paths are resolved
        # against the config file's directory, so a config is portable)
        "workdir": "work",
        "threads": 8,          # threads per sample-level job (STAR, bowtie2)
        "jobs": 8,             # samples processed concurrently
        "seed": 20260712,      # reservoir sampling / any shuffling
    },

    "pipeline": {
        # where to start: "query"      -- search the SRA/ENA for ribo-seq runs
        #                 "accessions" -- take a given list of run accessions
        #                 "fastq"      -- take a directory of local FASTQ files
        "start": "query",
        # where to stop:  "qc"           -- after the ribo-seq QC verdict (+ TSV, plots)
        #                 "architecture" -- after the read-architecture call (+ TSV, plots)
        #                 "bam"          -- after full download, trim, filter, mapping
        "end": "bam",
        "accession_list": None,   # start == "accessions": path to a one-per-line list
        "fastq_dir": None,        # start == "fastq": directory of *.fastq[.gz]
        # resume: skip any sample whose stage output already exists
        "resume": True,
        # carry only runs whose QC verdict is in this set into the later stages.
        # TI-SEQ is a ribo-seq variant, not junk, so it is kept by default.
        "keep_verdicts": ["RIBO-SEQ", "TI-SEQ"],
    },

    "reference": {
        "species": "Homo sapiens",
        "taxon_id": 9606,
        "genome_fasta": "/path/to/data/homo_sapiens.90.fasta",
        "gtf": "/path/to/data/homo_sapiens.90.gtf",
        "star_index": "/path/to/data/index_files/STAR-index",
        # bowtie2 index PREFIX of rRNA/tRNA/snRNA/snoRNA/Mt contaminant sequences.
        # null => built from `contaminant_fasta` into <workdir>/refs/ on first use.
        "contaminant_index": None,
        # null => the HUMAN contaminant reference bundled with RiboMine
        # (ribomine/data/human_riboseq_contaminants.rRNA_tRNA_snRNA_snoRNA_Mt.fa).
        # For another organism, point this at that organism's sequences.
        "contaminant_fasta": None,
        # cached GTF index used by the QC stage; null => <workdir>/refs/<gtf>.idx.pkl
        "annotation_index": None,
    },

    "query": {
        # "ena"    -- ENA portal search (one request, and it carries the fastq URLs)
        # "entrez" -- NCBI Entrez over db=sra (catches SRA-only and very fresh studies)
        # "both"   -- union of the two; neither is a superset of the other
        "source": "both",
        # SINGLE WORDS ONLY. ENA's text index is tokenised, so a multi-word wildcard
        # (`study_title="*ribosome profiling*"`) matches nothing while `"*ribosome*"`
        # matches thousands. These tokens are the RECALL net -- deliberately broad,
        # ~190k human runs, most of them Ribo-Zero RNA-seq. Precision is then applied
        # in ribomine/sra/query.py by the strong/weak/exclude phrase regexes, and the
        # final word belongs to the QC stage, which reads the actual periodicity.
        "terms": ["ribo", "ribosome", "riboseq", "RPF", "translatome",
                  "footprint", "harringtonine", "lactimidomycin", "monosome",
                  "ARTseq", "RiboLace"],
        "taxon_id": None,            # null => reference.taxon_id
        # a run must be one of these to be considered at all (ENA library_strategy).
        # Ribo-seq has no strategy of its own; submitters use these.
        "library_strategies": ["RNA-Seq", "OTHER", "ncRNA-Seq", "miRNA-Seq"],
        "instrument_platform": "ILLUMINA",   # "any" to drop the filter
        "layout": "SINGLE",          # SINGLE | PAIRED | any -- ribo-seq is single-end
        "min_read_count": 1_000_000,
        "max_runs": 0,               # 0 = no cap
        "runs_per_study": 0,         # 0 = all; N = the N deepest runs per study
        "exclude_runs": [],
        "exclude_studies": [],
        "ncbi_api_key": None,        # raises the Entrez rate limit from 3/s to 10/s
        "ncbi_email": None,
    },

    "download": {
        # The FULL-dataset download always goes to ENA over HTTPS -- it is several
        # times faster than any other route (docs/DOWNLOAD.md), so there is nothing
        # to choose. Runs ENA has not mirrored fall back to the SRA mirrors
        # automatically; that is availability, not a setting.
        "connections": 8,        # parallel connections per file (aria2c / range-GET)
        "max_retries": 4,
        "retry_backoff_s": 5,
        "tmpdir": None,          # null => <workdir>/tmp
    },

    # --- what survives the run ---------------------------------------------
    # Mining the SRA means hundreds of runs at 1-10 GB each, so the default is to
    # keep ONLY the deliverable BAM. Everything else here is an intermediate that
    # the pipeline can regenerate from the accession, and nothing that is deleted
    # is ever needed to READ the results: the JSONs, the TSVs and the plots -- the
    # full record of every number and every call -- are always kept.
    "keep": {
        "bam": True,             # the deliverable: trimmed, filtered, mapped reads
        "fastq": False,          # the raw run FASTQ, as downloaded
        "trimmed_fastq": False,  # after the architecture-driven trim
        "clean_fastq": False,    # after contaminant removal -- the reads that were mapped
        "sra": False,            # the .sra container (only the SRA fallback routes make one)
        # The QC / architecture stage's own working data: the 200k-read sample, the
        # contaminant-filtered sample, and their two local alignments. They are the
        # EVIDENCE for the verdict and the read-architecture call, so keeping them
        # lets you re-examine a call in a browser -- but the call itself, and every
        # number behind it, is already written to qc.json / profile.json / arch.json
        # and to the two TSVs, which is why they go by default.
        "qc_fastq": False,
        "qc_bam": False,
    },

    # --- stage 2: read sample + QC -----------------------------------------
    "qc": {
        "sample_reads": 200_000,     # reads kept per run for QC/architecture
        "scan_reads": 1_000_000,     # reads streamed to sample from (0 = whole run)
        "max_reads_scored": 200_000,
        # verdict thresholds (tuned on a 187-sample human cohort)
        "footprint_len_lo": 25,
        "footprint_len_hi": 36,
        "read_len_peak_frac_min": 0.40,
        "cds_enrich_min": 0.20,      # fraction of genic reads in CDS to look ribo-seq
        "cds_strong": 0.55,          # strong CDS enrichment excuses weak periodicity
        "periodic_min": 0.42,        # in-frame fraction (chance = 1/3)
        "periodic_strong": 0.50,
        "tvd_min": 0.10,             # TVD of the frame distribution to uniform
        "tvd_strong": 0.20,
        "min_cds_reads": 150,
        "single_locus_max": 0.50,    # >half the reads at one 5' locus = contaminant
        "locus_vs_cds": 1.0,         # one locus holding as many reads as the whole CDS
        "locus_min_to_judge": 0.05,
        "min_unique_frac": 0.15,     # below this the library maps too poorly to use
        "tiseq_ratio_min": 40,       # start-codon peak / CDS body => TI-seq
        "tiseq_min_peak": 200,
    },

    "contaminants": {
        "enabled": True,
        "min_entropy": 1.1,          # low-complexity read filter (bits)
        "max_base_frac": 0.85,
        # position pile-up filter (adapter dimers etc; data-driven, no adapter seqs)
        "pileup_min_count": 30,
        "pileup_min_frac": 0.005,
        "pileup_len_conc": 0.85,
        "pileup_max_frac": 0.10,
    },

    # --- stage 3: read architecture ----------------------------------------
    "architecture": {
        "min_reads": 5000,
        "min_plateau": 0.55,
        "footprint_uniform_max": 0.85,   # a single length spike = fixed-length artefact
        "single_locus_max": 0.50,
        "genomic_frac": 0.55,
        "umi_match_frac": 0.25,
        "ent_random": 1.70,
        "ent_const": 0.90,
        "footprint_pen": 0.50,
        "struct_pen": 0.85,
        "rt_max_len": 2,
        "rt_pen_hi": 0.65,
        "at_jump": 0.16,
        "at_abs": 0.66,
        "gc_template": 0.82,
        "gc_mid": 0.60,
        "linker_max": 10,
        "min_anchor_frac": 0.15,
        # runs whose architecture is `undetermined` still have a QC verdict; carry
        # them into processing with a best-effort trim plan, or drop them?
        "process_undetermined": False,
    },

    # --- stage 4: full download -> trim -> filter -> map --------------------
    "process": {
        "min_len": 20,               # drop reads shorter than this after trimming
        # A read that does not carry the library's adapter never reached the end of
        # its molecule: its 3' end is set by the read length, not by the footprint,
        # and any 3' UMI beyond it was never sequenced. Such reads are not complete
        # footprints (cutadapt's --discard-untrimmed drops them for the same reason),
        # and keeping them yields a SHORT UMI, which umi_tools refuses outright.
        # Set false to keep them anyway -- the missing UMI bases are then written as N.
        "discard_untrimmed": True,
        # nt of the fixed 3' scaffold (barcode + adapter) that must be visible in a read
        # before it is cut there. 7 nt of a known string, anchored at the read end,
        # matches by chance about once in 16,000 reads. Lower it when the library's
        # molecules barely fit the read -- the scaffold then runs off the end of most
        # reads and they cannot be trimmed at all (RiboMine warns when that happens).
        # It is a specificity/depth trade: a bad cut corrupts the footprint boundary.
        "adapter_min_overlap": 7,
        "umi_dedup": False,          # UMI dedup of the BAM -- OFF by default
        # "umicollapse" | "umi_tools". Same algorithms, same `_<UMI>` read-name tag,
        # same answer -- but UMICollapse indexes each position's UMIs in a BK-tree
        # instead of comparing every pair, and ribo-seq is where that matters: a
        # well-translated start codon piles hundreds of distinct UMIs onto one
        # coordinate, which is umi_tools' worst case. umi_tools is kept for its
        # `unique` and `percentile` methods, which UMICollapse does not have.
        "umi_dedup_tool": "umicollapse",
        # directional | adjacency | cluster (both tools), or unique | percentile
        # (umi_tools only)
        "umi_dedup_method": "directional",
        "umi_dedup_mem_gb": 8,       # JVM heap for umicollapse; ignored by umi_tools
        "filter_contaminants": True, # bowtie2 rRNA/tRNA removal before mapping
        "filter_pileups": True,      # data-driven pile-up removal after mapping
        # Every BAM RiboMine leaves on disk is coordinate-sorted and indexed -- the
        # deliverable ones and the QC-stage ones alike -- so any of them can be opened
        # in a genome browser without a further step. This is not optional: an
        # unindexed BAM is a BAM nobody can look at, and UMI deduplication requires a
        # sorted, indexed input anyway.
        # What is KEPT on disk afterwards is the `keep` block's business, not this one.
    },

    "mapping": {
        # end-to-end alignment of the TRIMMED reads: the real output BAM
        "align_ends_type": "EndToEnd",
        "multimap_nmax": 1,          # 1 = unique only; raise to keep multimappers
        "mismatch_nmax": 3,
        "mismatch_noverlmax": 0.1,
        "match_nmin_over_lread": 0.9,
        "align_intron_max": 0,       # 0 = STAR default (spliced); 1 = unspliced only
        "extra_args": [],
        # the QC/architecture stage's own alignment is LOCAL and permissive on
        # purpose (the adapter is an output, not an input) -- not user-tunable here.
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
        stem = os.path.splitext(os.path.basename(self.get("reference.gtf", "annotation.gtf")))[0]
        return os.path.join(self.dir("refs"), f"{stem}.idx.pkl")

    def ref(self, key: str) -> str | None:
        return self._abs(self.get(f"reference.{key}"))

    # -- validation ---------------------------------------------------------
    def validate(self, *, need_reference: bool = True, need_inputs: bool = True) -> None:
        """Check the config makes sense for what is about to run.

        `need_inputs` is False for sub-commands that do not consume the pipeline's
        start point -- `ribomine query` searches the archive and has no use for an
        accession list, so demanding one (which the query is often the thing that
        produces) would be absurd.
        """
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
    """The stages a run executes, given its start and end point.

    `start` selects how the sample list is obtained (query / given list / local
    FASTQs); the analysis stages that follow are always the same, truncated at
    `end`.
    """
    first = 0 if start == "query" else 1          # a list/dir skips the query stage
    last = STAGES.index({"qc": "qc", "architecture": "architecture", "bam": "bam"}[end])
    return list(STAGES[first : last + 1])


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
            + "\n(a typo here would otherwise be silently ignored)"
        )
    data = deep_merge(DEFAULTS, user)
    data = deep_merge(data, overrides or {})
    return Config(data=data, path=path or "")


def _unknown_keys(user: dict, schema: dict, prefix: str = "") -> list[str]:
    """Keys in `user` that the schema does not define (typo guard).

    Keys starting with `_` are comments and are ignored -- JSON has no comment
    syntax, and a config you cannot annotate is a config nobody can read.
    """
    bad = []
    for k, v in (user or {}).items():
        if k.startswith("_"):
            continue
        if k not in schema:
            bad.append(f"{prefix}{k}")
        elif isinstance(v, dict) and isinstance(schema[k], dict):
            bad += _unknown_keys(v, schema[k], f"{prefix}{k}.")
    return bad


# What `null` means, per key. In JSON a null reads as "nothing / switched off", but
# for every key below it means "work it out for me" -- the opposite. That is not a
# distinction a reader can make from the file alone, so `write_example` spells it out
# in the file itself.
NULL_MEANS = {
    "reference.contaminant_fasta":
        "use the HUMAN contaminant reference bundled with RiboMine "
        "(rRNA/tRNA/snRNA/snoRNA/Mt + 45S pre-rRNA + rDNA repeat). Contaminant "
        "filtering is ON either way -- set contaminants.enabled=false to turn it off.",
    "reference.contaminant_index":
        "build the bowtie2 index from contaminant_fasta into <workdir>/refs/ on first "
        "use (takes ~1s). Set this only to reuse an index you already have.",
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
    """Write a config with every default -- and with what `null` means, inline.

    The `_null_means` block is a comment: `load()` ignores any key starting with `_`.
    """
    doc = {
        "_README": [
            "Every key below is RiboMine's own default; delete anything you do not "
            "want to change. Keys starting with '_' are comments and are ignored.",
            "An unknown key is an ERROR, not a silent no-op -- a typo'd threshold "
            "that gets ignored is worse than a crash.",
        ],
        "_null_means": NULL_MEANS,
        **DEFAULTS,
    }
    with open(path, "w") as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")
