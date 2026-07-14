# Internal module contract

Every stage reads the previous stage's JSON out of the per-sample directory and
writes its own. `ribomine.utils.Sample` owns all the filenames; no module builds
a path by string concatenation.

Thresholds come from the config (`ribomine.config.Config`), never from module
constants — `cfg.get("qc.periodic_min")`. Modules take a `cfg` (or a small
thresholds object built from it) rather than reaching for a global.

Nothing here prints to stdout for control flow: modules return dicts, the
pipeline logs. `LOG = logging.getLogger("ribomine.<module>")`.

## ribomine.sra.query

```python
def search(cfg: Config) -> list[dict]
    """Candidate ribo-seq runs. Rows use ENA read_run field names
    (run_accession, study_accession, experiment_accession, sample_accession,
     scientific_name, tax_id, library_strategy, library_selection, library_layout,
     library_construction_protocol, instrument_platform, instrument_model,
     read_count, base_count, fastq_ftp, fastq_bytes, study_title,
     experiment_title, sample_title, first_public, center_name, source)
    `source` is 'ena', 'entrez' or 'ena+entrez'."""

def run_query(cfg: Config) -> str
    """search() + filters + dedup -> writes <workdir>/meta/candidates.tsv, returns path."""
```

## ribomine.sra.metadata

```python
def fetch_runs(accessions: list[str], batch: int = 100) -> list[dict]   # ENA portal, same fields
def fastq_urls(acc: str) -> list[str]                                   # https URLs, R1 first
def read_count(acc: str) -> int | None
```

## ribomine.sra.download

```python
def sample_reads(source: str, out_fastq: str, *, n: int, scan: int, seed: int) -> dict
    """`source` is an SRR/ERR/DRR accession (streamed, nothing stored) or a local
    fastq(.gz) path. Reservoir-samples `n` reads out of the first `scan`
    (scan=0 => whole file). Returns {'n_sampled', 'n_scanned', 'sorted_warning', 'source'}."""

def download_full(acc: str, out_fastq_gz: str, cfg: Config, *, log: str = "") -> dict
    """Full run over ENA+aria2c, falling back to the SRA mirrors only for runs ENA
    has not mirrored (see docs/DOWNLOAD.md). Returns {'route', 'bytes', 'seconds',
    'mb_per_s', 'attempts'}. Idempotent: returns immediately if out_fastq_gz
    already exists and is non-empty."""
```

## ribomine.qc.annotation

```python
def build_index(gtf: str, out_pkl: str) -> str
def ensure_index(cfg: Config) -> str    # builds cfg.annotation_index if absent, returns it
```

## ribomine.qc.contaminants

```python
def filter_fastq(in_fq: str, out_fq: str, cfg: Config, *, label: str = "", threads: int = 8,
                 log: str = "") -> dict
    """bowtie2 --very-sensitive-local against the contaminant index (skipped with a
    warning if reference.contaminant_index is unset) + a low-complexity/homopolymer
    filter. Keeps FULL untrimmed reads. Returns the stats dict written to
    Sample.contam_json: n_input, n_contaminant_rRNA_tRNA_etc, n_low_complexity,
    n_kept, frac_contaminant_structured_rna, frac_low_complexity, frac_kept."""
```

## ribomine.qc.pileups

```python
def filter_bam(in_bam: str, out_bam: str, cfg: Config, *, label: str = "") -> dict
    """Data-driven position pile-up removal (no adapter sequences). Returns the
    dict written to Sample.pileup_json: n_reads_in, n_reads_removed, n_reads_kept,
    frac_removed, n_pileup_positions, n_dominant_positions, top_pileups[]."""
```

## ribomine.qc.profile

```python
def profile_bam(bam: str, fasta: str, *, label: str, max_reads: int = 150_000) -> dict
    """Positional match/composition profile. Returns the dict written to
    Sample.profile_json (keys as in the reference implementation: p5_match,
    p5_comp, clip5_hist, adap_match, adap_comp, anchor_kind, anchor_offset,
    n_anchored, frac_anchored, t3_match, t3_comp, footprint_len_hist,
    read_len_hist, fpend_to_adapter_gap_hist, top5p_locus_frac,
    n_distinct_5p_loci, panel_hits, seed_kmer, n_used, ...)."""
```

## ribomine.qc.verdict

```python
def qc(bam: str, index_path: str, cfg: Config, *, star_log: str = "", label: str = "",
       contam: dict | None = None, pileup: dict | None = None,
       total_reads: int | None = None) -> dict
    """Three-way verdict: 'RIBO-SEQ' | 'TI-SEQ' | 'NOT RIBO-SEQ or LOW QUALITY',
    plus is_riboseq, verdict_reason, reasons[], signals{}, periodicity_inframe_frac,
    periodicity_tvd_uniform, cds_frac_of_genic, top5p_locus_frac, read_len_mode,
    region_frac{}, metagene_start{}, metagene_stop{}, frame_by_len{}, mapping{},
    contaminants{}, usable{projected_usable_reads,...}. Written to Sample.qc_json."""

def periodicity(bam: str, index_path: str, *, max_reads: int = 200_000,
                label: str = "") -> dict
    """3-nt periodicity of a FINISHED BAM -- trimmed, filtered, deduplicated. Decides
    nothing (qc() makes the call, on the untrimmed QC sample); this is the measurement
    OF THE OUTPUT, and the only place a mis-trimmed footprint boundary shows up.
    Reads are taken at a fixed stride, not from the front -- the BAM is coordinate-
    sorted, so its first reads are its first chromosome. Returns {'n_reads_in_bam',
    'n_reads_scored', 'stride', 'read_len_mode', 'periodicity_inframe_frac',
    'periodicity_tvd_uniform', 'n_cds_reads', 'cds_frac_of_genic', 'region_frac'}."""
```

## ribomine.qc.plot / ribomine.arch.plot

```python
def plot_qc(qc: dict, out_path: str, *, dpi: int = 110) -> str
def plot_arch(profile: dict, call: dict, out_path: str, *, dpi: int = 110) -> str
```

## ribomine.arch.infer

```python
@dataclass
class Thresholds:            # built from cfg["architecture"], one field per config key
    ...
    @classmethod
    def from_config(cls, cfg: Config) -> "Thresholds"

def infer(profile: dict, thr: Thresholds | None = None) -> dict
    """The architecture call. status 'ok' | 'undetermined' | 'error'; on ok also a
    `functional` block with trim_5p, footprint_retains_rt_nt, trim_3p_adapter,
    trim_3p_construct, dedup_umi_len, dedup_umi_mask, segments_5p[], segments_3p[].
    Written to Sample.arch_json."""

def architecture_string(call: dict) -> str    # "5'-[UMI,2nt]-[footprint,~30nt]-[TruSeq]-3'"
```

## ribomine.arch.trim

```python
def trim_fastq(in_fq: str, call: dict, out_fq: str, *, min_len: int = 20,
               label: str = "") -> dict
    """Apply the inferred architecture: trim the 5' construct (keeping the enzymatic
    RT base), trim the 3' construct + adapter / poly(A) tail, and move the
    random-templated UMI content into the read name as `_<UMI>` (umi_tools style).
    Handles .gz in and out. Returns the dict written to Sample.trim_json."""
```

## ribomine.process.star

```python
def align_local(fastq: str, outdir: str, cfg: Config, *, threads: int = 8, log: str = "") -> str
    """Permissive LOCAL alignment of untrimmed reads (QC/architecture stage). The
    5' UMI / RT nt / 3' UMI / barcode / adapter land in the soft clips. Returns the BAM."""

def align_final(fastq: str, outdir: str, cfg: Config, *, threads: int = 8, log: str = "") -> str
    """End-to-end alignment of the TRIMMED reads (the deliverable BAM), using the
    cfg['mapping'] parameters. Also writes ReadsPerGene.out.tab (--quantMode GeneCounts)
    when the index has the annotation. Returns the BAM."""

def has_annotation(cfg: Config) -> bool  # was the STAR index built with --sjdbGTFfile?
def sort_index(bam: str, out_bam: str, *, threads: int = 4) -> str
def parse_log(star_log: str) -> dict    # n_input, n_unique, n_multi, frac_* (unique+multi+unmapped = 1)
```

## ribomine.process.counts

```python
def path(star_final_dir: str) -> str    # where STAR left this run's ReadsPerGene.out.tab

def read_counts(tab: str, *, label: str = "") -> tuple[dict[str, int], dict]
    """({gene_id: count}, stats). Counts are STAR's SENSE column (a footprint is a
    piece of the mRNA). stats = {'n_in_genes', 'n_genes_detected', 'n_no_feature',
    'n_ambiguous', 'n_antisense', 'sense_over_antisense', 'frac_in_genes'} -- i.e.
    what the counts do NOT contain. Warns when the library does not read as
    sense-stranded. STAR's N_multimapping row is deliberately NOT reported: under
    multimap_nmax=1 it reads 0 no matter how multimapping the library is, because the
    multimappers were dropped before counting. Log.final.out has that number."""

def matrix(out_tsv: str, star_index: str, columns: list[tuple[str, str]]) -> str
    """Join per-run count tables into the genes x runs matrix. Rows come from the
    index's geneInfo.tab -- every gene, including the all-zero ones, so two matrices
    are comparable. A run with no counts is left OUT as a column (a blank is not a 0)."""
```

Counts are taken DURING the alignment, so they precede the pile-up filter and any UMI
dedup, and they are over the whole gene (all exons, any biotype), not the CDS.

## ribomine.process.dedup

```python
def dedup(bam: str, out_bam: str, cfg: Config, *, log: str = "") -> dict
    """UMI dedup of a coordinate-sorted+indexed BAM whose read names carry `_<UMI>`.
    Backend is process.umi_dedup_tool: "umicollapse" (default) or "umi_tools" -- same
    algorithms, same answer, but UMICollapse indexes each position's UMIs instead of
    comparing every pair, which is what ribo-seq's deep positions need.
    OFF by default (process.umi_dedup). Returns {'n_in','n_out','frac_kept','method','tool'}.
    Raises a clear error if the reads carry no UMI (architecture found none)."""
```

## ribomine.reports

```python
def qc_tsv(cfg: Config, accs: list[str]) -> str            # -> <workdir>/qc/qc_summary.tsv
def arch_tsv(cfg: Config, accs: list[str]) -> str          # -> <workdir>/architecture/architecture.tsv
def process_tsv(cfg: Config, accs: list[str]) -> str       # -> <workdir>/mapping_summary.tsv
```

## ribomine.pipeline

```python
def run(cfg: Config) -> dict
    """Drives stages_to_run(start, end). Per-sample work runs in a process pool of
    cfg['project.jobs']; a sample that fails is logged, recorded in
    <workdir>/failed.tsv, and does not take the run down. `resume` skips a sample
    whose stage output already exists."""
```
