# RiboMine

Mine Ribo-seq datasets out of the SRA: search the archive, decide from the reads
which hits really are ribosome profiling, work out how each library was built,
then download, trim and map the ones worth keeping.

```
search SRA/ENA ─► sample reads ─► QC verdict ─► read architecture ─► download ─► trim ─► map
                  (streamed)      RIBO-SEQ /    UMIs, barcode,                  UMIs to   BAM +
                                  TI-SEQ /      adapter (fqdissect)             read name counts
                                  rejected
```

Neither question can be answered from metadata. Ribo-seq has no
`library_strategy` of its own (it is deposited as `RNA-Seq` or `OTHER`), and
protocol descriptions of adapters and UMIs are often incomplete or wrong. RiboMine
therefore measures both on a sample of the reads: 3-nt periodicity and CDS
enrichment for the verdict, and the positions of UMIs, barcodes and adapter for the
read architecture. When the reads do not support an answer, it says so instead of
guessing.

RiboMine was built to assemble the Ribo-seq cohort for PRICE 2.

## Install

```bash
git clone https://github.com/kilian-m/RiboMine.git && cd RiboMine
conda env create -f environment.yml
conda activate ribomine
pip install -e .            # also installs fqdissect
```

You need a genome FASTA, its GTF (Ensembl-style, with CDS frames), and a STAR index
built with that GTF (`--sjdbGTFfile`; without it there are no gene counts). The
human contaminant reference (rRNA, tRNA, snRNA, snoRNA, Mt, 45S pre-rRNA, rDNA
repeat) is bundled; for another organism set `reference.contaminant_fasta`.

## Quick start

```bash
ribomine init-config config.json     # every setting, with its default
$EDITOR config.json                  # set reference.genome_fasta, .gtf, .star_index
ribomine setup -c config.json        # build the annotation and contaminant indexes

ribomine run -c config.json          # search the archive and process the hits
```

A run can start and stop at different points:

| `--from` | input |
|---|---|
| `query` | search SRA/ENA for candidate runs (`ribomine query` runs only this) |
| `accessions` | a list of run accessions: `--accessions runs.txt` |
| `fastq` | a directory of single-end FASTQs: `--fastq-dir reads/` |

| `--to` | output |
|---|---|
| `qc` | `qc/qc_summary.tsv` and a QC figure per run |
| `architecture` | `architecture/architecture.tsv` and an architecture figure per run |
| `bam` | `bams/*.bam`, `mapping_summary.tsv`, `counts/gene_counts.tsv` |

```bash
ribomine run -c config.json --to qc                          # screen only
ribomine run -c config.json --accessions runs.txt --to architecture
ribomine run -c config.json --fastq-dir reads/ --umi-dedup   # local FASTQs to BAMs
```

`--to qc` never downloads a full run, so screening is cheap; do it first on a large
search. Runs are resumable: a sample whose stage output exists is skipped, and a
failing sample is recorded in `failed.tsv` without stopping the others.

## Output

Everything is written under `project.workdir`.

**`qc/qc_summary.tsv`**: one row per run, with the verdict, its reason, and the
measurements behind it. The figure in `qc/plots/` shows the same signals with the
thresholds drawn in.

```
run_accession  verdict                      read_len_mode  periodicity_inframe  cds_frac_of_genic
SRR24210493    RIBO-SEQ                     30             0.44                 0.60
SRR30357177    NOT RIBO-SEQ or LOW QUALITY  32             0.54                 0.64
```

**`architecture/architecture.tsv`**: what each read is made of, how much to trim,
and how many nt are random and usable for deduplication.

```
run_accession  architecture                                                                        trim_5p  dedup_umi_len
SRR24210493    5'-[UMI,2nt]-[footprint,~30nt]-[UMI,5nt]-[barcode,TAGAC]-[TruSeq,AGATCGGAAGAG…]-3'  2        7
```

**`bams/`**: coordinate-sorted, indexed BAMs of the trimmed, contaminant-filtered
reads, with UMIs in the read name (`@name_UMI`).

**`mapping_summary.tsv`**: per run, what was removed at each step (trimming,
contaminants, pile-ups, duplicates), the mapping rates, and the periodicity of the
finished BAM. That periodicity is measured independently of the QC verdict; if it
is clearly worse than the QC value, suspect the trim.

**`counts/gene_counts.tsv`**: genes × runs matrix of STAR's gene counts
(`--quantMode GeneCounts`, sense strand, whole gene). They are taken during
alignment, so before pile-up removal and deduplication.

By default only the BAMs, tables, figures and per-sample JSONs are kept; FASTQs and
the QC alignments are deleted. The `keep` block of the config changes that.

## How it works

**QC verdict.** 200,000 reads are sampled from the first million (streamed from
ENA, nothing stored), filtered against the contaminant reference with bowtie2, and
aligned untrimmed with STAR in local mode. The verdict rests on the reading frame
of the read 5′ ends within annotated CDS, supported by CDS enrichment, read length,
the start-codon profile and the unique-mapping rate. TI-seq (harringtonine,
lactimidomycin) is recognised by its start-codon peak and kept. Details and
thresholds: [docs/QC.md](docs/QC.md).

**Read architecture.** The same local alignment is passed to
[fqdissect](https://github.com/kilian-m/fqdissect). Non-genomic sequence ends up in
the soft clips, so the position at which read bases start and stop matching the
genome gives the boundaries of the footprint; base composition on either side tells
UMI from barcode from adapter. The adapter is found de novo. A library that cannot
be read is reported as `undetermined` with a reason and is not processed further.
The method is described in fqdissect's
[METHOD.md](https://github.com/kilian-m/fqdissect/blob/main/docs/METHOD.md).

**Processing.** The full run is downloaded from ENA
([docs/DOWNLOAD.md](docs/DOWNLOAD.md)), trimmed with cutadapt according to the
architecture, filtered for contaminants, aligned with STAR, and cleared of
single-position pile-ups (adapter dimers, miRNAs). UMI deduplication is off by
default and refused for libraries without a UMI, because identical footprints are
then real signal (`--umi-dedup`, UMICollapse or umi_tools).

## Configuration

One JSON file describes a run; command-line flags override it.
`ribomine init-config` writes every key with its default, and
[`ribomine/config.py`](ribomine/config.py) documents them. An unknown key is an
error. Keys starting with `_` are comments.

`null` usually means "use the default", not "off": for example
`reference.contaminant_fasta: null` selects the bundled human reference. The
generated config explains each such key in its `_null_means` block.

## Running a cohort on SLURM

`slurm/` contains the job chain used for the full cohort on LRZ CoolMUC-4
(search and split, one RiboMine per node, merge). See
[docs/SLURM.md](docs/SLURM.md).

```bash
slurm/master.sh config/config_lrz.json
```

## Repository layout

```
ribomine/
  cli.py            command line
  config.py         settings and defaults
  pipeline.py       stages, process pool, resume
  architecture.py   hand-off to fqdissect (profile, call, figure, trim)
  reports.py        the TSV tables
  sra/              archive search, metadata, read sampling, download
  qc/               annotation index, contaminant filter, pile-up filter, verdict, figure
  process/          STAR alignment, gene counts, UMI deduplication
  data/             bundled human contaminant reference
slurm/, scripts/    SLURM job chain; cohort split and merge
docs/               QC.md, DOWNLOAD.md, SLURM.md
tests/              pytest suite (no genome needed)
```

## License

MIT, see [LICENSE](LICENSE).
