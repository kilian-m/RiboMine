# QC: is it Ribo-seq, and what is removed

The QC stage works on a sample of each run: 200,000 reads drawn uniformly from the
first million (`qc.sample_reads`, `qc.scan_reads`; `scan_reads: 0` samples the
whole run). Archive FASTQs are in spot order, so this prefix is representative
unless the deposit was sorted, which is detected and logged.

```
sample ─► contaminant filter ─► STAR local alignment ─► pile-up filter ─► verdict
                                                                      └─► profile ─► architecture (fqdissect)
```

The reads are aligned untrimmed, in local mode, so that UMIs, barcodes and adapter
end up in the soft clips. The verdict and the read architecture are both measured
on this alignment.

## Contaminant removal

Before alignment (`qc/contaminants.py`):

| removed | how |
|---|---|
| rRNA, tRNA, snRNA, snoRNA, mitochondrial rRNA/tRNA | `bowtie2 --very-sensitive-local` against the contaminant FASTA |
| low-complexity reads (poly-A/G, simple repeats) | Shannon entropy below `contaminants.min_entropy` (1.1 bits) or one base above `contaminants.max_base_frac` (85 %) |
| reads STAR cannot take (longer than 300 nt, malformed records) | dropped while screening |

The bundled FASTA is human and includes the 45S pre-rRNA and the rDNA repeat,
whose transcribed spacers are absent from mature rRNA sequences. Without this
filter rRNA reads map to hundreds of genomic copies and inflate the multimapping
rate.

After alignment (`qc/pileups.py`): adapter dimers, miRNAs and other fixed
molecules pile up on one genomic 5′ position. No adapter sequences are needed to
find them. A 5′ position is removed when

* it holds at least `pileup_min_frac` (0.5 %) of the reads and `pileup_min_count`
  (30), **and** at least `pileup_len_conc` (75 %) of its reads have one length
  (footprints at a translated codon spread over lengths; a single molecule does
  not), or
* it holds at least `pileup_max_frac` (10 %) of the remaining reads, whatever
  their lengths. This is applied until no such position is left. In libraries
  that pass QC the busiest position holds a median of 0.5 % of the reads.

Both filters run again on the full run in the `bam` stage.

## The verdict

Measured on up to `qc.max_reads_scored` alignments against a cached index of the
GTF (`qc/annotation.py`; it expects Ensembl feature names such as `five_prime_utr`
and the `gene_biotype` attribute):

| signal | definition |
|---|---|
| periodicity | frame of the read 5′ end within annotated CDS, each read length aligned to its own dominant frame. Reported as the in-frame fraction (chance 1/3) and as the total variation distance of the frame distribution to uniform (0 to 0.67) |
| CDS enrichment | fraction of genic reads in CDS |
| region composition | CDS, 5′UTR, 3′UTR, ncRNA, intron, intergenic, mitochondrial |
| read length | mode, and fraction within `footprint_len_lo`–`footprint_len_hi` (25–36 nt) |
| start-codon ratio | 5′-end density at the start-codon peak over the in-frame CDS body downstream |
| top 5′ locus | share of reads at the single busiest 5′ position |
| unique mapping | uniquely mapped fraction, from the STAR log |

A library is **translating** if

* no single 5′ locus holds more than `single_locus_max` (50 %) of the reads, or as
  many reads as the whole CDS (`locus_vs_cds`), and
* at least `cds_enrich_min` (20 %) of genic reads are in CDS, and
* either periodicity is strong (in-frame ≥ `periodic_strong` 0.50 and TVD ≥
  `tvd_strong` 0.20, on at least `min_cds_reads` 150 CDS reads), or the read
  length is footprint-like and periodicity is at least weak (in-frame ≥
  `periodic_min` 0.42, TVD ≥ `tvd_min` 0.10) or CDS enrichment is strong
  (≥ `cds_strong` 0.55).

It is **usable** if, in addition, at least `min_unique_frac` (15 %) of reads map
uniquely, CDS makes up at least `cds_region_min` (50 %) of all reads, and at least
`read_len_min_frac` (75 %) of reads are footprint-length.

| verdict | condition |
|---|---|
| `RIBO-SEQ` | usable |
| `TI-SEQ` | usable, and the start-codon peak holds at least `tiseq_min_peak` (200) reads and is at least `tiseq_ratio_min` (30) times the CDS body. Initiation inhibitors (harringtonine, lactimidomycin) produce this; elongating libraries stay below about 30 |
| `NOT RIBO-SEQ or LOW QUALITY` | everything else; `verdict_reason` says which condition failed |

`pipeline.keep_verdicts` selects which verdicts continue to the architecture and
`bam` stages (default: `RIBO-SEQ` and `TI-SEQ`). All thresholds are in the `qc`
block of the config and were tuned on a cohort of 187 human runs.

Two further observations are reported but do not change the rules above:

* **Antisense deposit**: reads on CDS are mostly on the opposite strand. Footprints
  are sense, so this is either stranded RNA-seq or a reverse-complemented deposit.
* **Mitochondrial-dominated**: at least 30 % of the reads are mitochondrial
  (mitoribosome profiling). The verdict is decided on the nuclear CDS reads; the periodicity of
  the mitochondrial CDS reads is reported separately.

`projected_usable_reads` scales the reads that survived all filters from the sample
to the whole run.

## Periodicity of the finished BAM

`mapping_summary.tsv` repeats the periodicity measurement on the BAM the pipeline
delivers, after trimming, filtering and deduplication. The QC value says whether
the library is Ribo-seq; this one says whether the processed reads still are. A
footprint boundary cut in the wrong place shows up as a drop between the two.

## The read architecture

The profile and the architecture call are made by
[fqdissect](https://github.com/kilian-m/fqdissect) from the pile-up-filtered local
alignment; see its
[METHOD.md](https://github.com/kilian-m/fqdissect/blob/main/docs/METHOD.md). Its
thresholds can be overridden in the `architecture` block of the config.
