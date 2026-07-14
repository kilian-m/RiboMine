# RiboMine

**Mine Ribo-seq datasets out of the SRA.** Search the archive, work out which hits
are *actually* ribosome profiling, work out how each library was built, then
download, trim and map the ones worth keeping.

```
query the SRA ──► sample reads ──► QC & classify ──► read architecture ──► download ──► trim ──► map
                  (streamed,        RIBO-SEQ /        what each read is       (full)      (UMIs →   (BAM)
                   nothing           TI-SEQ /          made of                            header)
                   stored)           reject
```

The point of the middle two stages is that **neither the archive's metadata nor the
paper can be trusted** for either question.

*Is it Ribo-seq?* Ribo-seq has no `library_strategy` of its own — submitters deposit
it as `RNA-Seq` or `OTHER`, the same values a plain RNA-seq run carries. So RiboMine
does not ask the metadata. It asks the reads: 3-nt periodicity within CDS, CDS
enrichment, footprint length, the start-codon metagene. Translating ribosomes step
one codon at a time, and nothing else does.

*How was the library built?* The adapter, the UMIs, the barcode, the non-templated
base the reverse transcriptase adds — published protocol descriptions get these
wrong often enough that you cannot process a cohort on documentation. RiboMine reads
the architecture **out of the reads themselves** (§ *How it works*), and refuses to
answer rather than guess when they cannot support one.

---

## Quick start

```bash
conda env create -f environment.yml
conda activate ribomine
pip install -e .

ribomine init-config config.json     # every default, documented, in one file
$EDITOR config.json                  # point `reference` at your genome/GTF/STAR index
ribomine setup -c config.json        # build the annotation + contaminant indexes (once)

ribomine run -c config.json          # the whole thing
```

The **human contaminant reference is bundled** (3,858 sequences: rRNA / tRNA / snRNA /
snoRNA / Mt, plus the 45S pre-rRNA and the rDNA repeating unit — those two carry the
transcribed spacers, which are excised during rRNA maturation and so appear in no
mature sequence), so the only thing you have to supply is the genome, the GTF and a
STAR index. That matters more than it sounds: a Ribo-seq library that is *not*
contaminant-filtered looks like ~78 % multimapping junk, because every rRNA has
hundreds of genomic copies and so every rRNA read maps to "too many loci". For a
non-human organism, point `reference.contaminant_fasta` at that organism's sequences.

**Every BAM RiboMine writes is coordinate-sorted and indexed** — the deliverables in
`bams/` and the QC-stage alignments alike — so when a verdict looks wrong you can open
the exact alignment it was computed from. There is deliberately no option to turn that
off.

## Start and end points

You rarely want the whole pipeline. Both ends move:

| start (`pipeline.start`) | |
|---|---|
| `query` | search the SRA/ENA for Ribo-seq runs |
| `accessions` | take a list of run accessions (`--accessions runs.txt`) |
| `fastq` | take a directory of local FASTQs (`--fastq-dir reads/`) |

| end (`pipeline.end`) | deliverable |
|---|---|
| `qc` | **`qc/qc_summary.tsv`** + a QC figure per dataset |
| `architecture` | **`architecture/architecture.tsv`** + an architecture figure per dataset |
| `bam` | **`bams/*.bam`** + `mapping_summary.tsv` |

```bash
ribomine run -c config.json --to qc                       # just screen the archive
ribomine run -c config.json --accessions runs.txt --to architecture
ribomine run -c config.json --fastq-dir reads/ --to bam --umi-dedup
```

Everything is resumable: a dataset whose stage output already exists is skipped, so
an interrupted 500-dataset run picks up where it stopped.

## What comes out

**After QC** — one row per dataset, and a figure showing every signal the verdict was
made on, with the thresholds drawn on top, so the call is auditable rather than
asserted:

```
run_accession  verdict                     periodicity_inframe  cds_frac_of_genic  read_len_mode
SRR12285169    RIBO-SEQ                    0.46                 0.78               30
SRR618773      TI-SEQ                      0.73                 0.32               29
SRR1039508     NOT RIBO-SEQ or LOW QUALITY 0.36                 0.30               63
```

TI-seq (harringtonine / LTM / QTI- / GTI-seq) is **kept, not discarded** — it is
Ribo-seq with the ribosomes frozen at start codons, and it is detected as such (the
start-codon peak towers over an empty CDS body).

**After architecture** — what each read is made of, and the two numbers that drive
processing (how much to trim off the 5' end, and how many nt of the read are random
and therefore usable for deduplication):

```
run_accession  architecture                                                    trim_5p  dedup_umi_len
SRR12285169    5'-[UMI,2nt]-[footprint,~30nt]-[UMI,5nt]-[barcode,AGCTA]-[TruSeq]-3'   2        7
```

That row is exactly the documented McGlincy & Ingolia 2017 structure — recovered from
the reads, with no protocol given to the tool.

**After mapping** — sorted, indexed BAMs, plus a summary of what was thrown away and
why (contaminants, pile-ups, reads left too short, duplicates).

## How it works

### Is it Ribo-seq? (`qc`)

The reads are streamed from ENA — never downloaded in full — contaminant-filtered,
and aligned. Then six signatures decide, of which one is the hallmark: the **frame of
the footprint's 5' end within the CDS**. Translating ribosomes advance one codon at a
time, so their footprints pile into one frame; RNA-seq is flat. Everything else (CDS
enrichment, footprint length ~28–32 nt, the start/stop metagene, region composition,
unique-mapping rate) supports it.

The verdict is deliberately **inclusive** — strong periodicity alone is decisive — so
that Ribo-seq variants with odd footprint lengths or unusual chemistry survive. Only
libraries that are contaminant-dominated, non-periodic, or map too poorly to use are
rejected, and each rejection states its reason.

### What is the read made of? (`architecture`)

A Ribo-seq read is a footprint wrapped in scaffolding, any part of which may be
absent:

```
[5' UMI][RT nt][=========== footprint ===========][3' UMI][barcode][adapter]
  random  A/T        genomic, variable-length        random   fixed    fixed
```

**Nothing is trimmed before aligning.** The raw read is aligned with
`STAR --alignEndsType Local`, so everything non-genomic is pushed into the soft
clips instead of preventing the alignment. The adapter is an *output* of this
pipeline, not an input.

The one idea the whole thing rests on: STAR's clip boundary is **not** where the
footprint starts. If the last base of a 5' UMI happens to match the genome base
beside the footprint — probability ¼ — STAR extends the alignment and the clip
shrinks. *Clip lengths leak.* What does not leak is the genomic coordinate the
alignment implies for **read position 0**: extending leftwards decrements `pos` and
`clip5` together, so `pos - clip5` is invariant. Anchor there, and every read position
maps to a fixed genome base. Then ask the only question that matters:

> does read base *p* match the genome base implied for it?

Footprint bases match ~97 % of the time. UMI, adapter and barcode bases are not
genomic, so they match ~25 % by chance. **The transition between those two regimes is
the architecture.** `docs/ARCHITECTURE.md` is the full walk-through — including how a
non-templated RT base (present in only a *fraction* of molecules) is told apart from a
UMI, and why a lone non-genomic 5' base is never a 1-nt UMI.

When the reads cannot support a call — too few alignments, no genomic plateau, a
single-length footprint spike (adapter dimers), a majority of reads at one locus — the
answer is `undetermined` **with a reason**, never a guess.

### Then: download, trim, map (`bam`)

The full run is downloaded by the fastest route that can serve it, trimmed with the
architecture that was just inferred (the RT base is *kept* — it is part of the
molecule; the UMI content goes into the read name), contaminant-filtered, and aligned
end-to-end.

**UMI deduplication is off by default** (`--umi-dedup` / `process.umi_dedup`). It is
only meaningful if the library actually carries a UMI, and on a library that does not,
it silently collapses genuine duplicate footprints — which in Ribo-seq are real
signal, since a heavily translated codon *is* covered many times. RiboMine refuses to
run it on a library whose architecture found no random-templated content.

## Downloading

Route matters more than bandwidth, so RiboMine picks the route for you: **ENA over
HTTPS with 16 parallel connections**, which measured 7–18× faster than the `prefetch`
route most pipelines default to, because a single TCP stream is throttled server-side
at ~12 MB/s and `prefetch` opens exactly one. There is no route setting — nothing
about that ranking is site-specific enough to be worth a knob. The md5 ENA hands back
in the same call as the URL is verified.

ENA does not mirror quite everything (published human Ribo-seq is 99.97 % covered;
dbGaP and the last few weeks of releases are not), so runs it does not have fall back
to AWS Open Data and then the SRA toolkit, automatically. That chain is about
*availability*, not speed.

`docs/DOWNLOAD.md` has the numbers, the coverage measurements, and the landmines —
including that `prefetch --max-size` **exits 0 when it skips an oversized run**, and
that SRA Lite carries fake quality scores.

## Configuration

One JSON file describes a run; `ribomine init-config` writes it with every default.
CLI flags override it. `config.py` *is* the schema — an unknown key is an error, not a
silent no-op, because a typo'd threshold that gets ignored is worse than a crash.

**`null` means "work it out for me", not "off".** JSON has no way to say that, so the
generated config carries a `_null_means` block spelling out each one (keys starting
with `_` are comments and are ignored). The one that catches people:
`reference.contaminant_fasta: null` selects the **bundled human reference** — it does
*not* disable contaminant filtering. To actually disable it, set
`contaminants.enabled: false`. Either way the run log says which reference it used.

RiboMine is human-first but not human-only: the `reference` block (genome FASTA, GTF,
STAR index, contaminant index, taxon) is all that ties it to a species.

## Layout

```
ribomine/
  cli.py            the command line
  config.py         the schema, the defaults, and their documentation
  pipeline.py       stage orchestration, the process pool, resume
  sra/
    query.py        find Ribo-seq runs (and why that takes three steps)
    metadata.py     ENA portal records
    download.py     read sampling (streamed) + the full-run route chain
  qc/
    annotation.py   GTF -> cached interval index
    contaminants.py rRNA/tRNA/snRNA/Mt (bowtie2) + low-complexity removal
    profile.py      the positional match/composition statistics
    verdict.py      is it Ribo-seq?
    pileups.py      data-driven pile-up removal (no adapter sequences)
    plot.py         the QC figure
  arch/
    infer.py        the architecture caller (every threshold lives here)
    trim.py         apply the architecture; UMIs -> read name
    plot.py         the architecture figure
  process/
    star.py         alignment (local for reading architecture, end-to-end for the BAM)
    dedup.py        optional UMI deduplication
  data/             the bundled human contaminant reference (rRNA/tRNA/snRNA/snoRNA/Mt)
  reports.py        the TSVs
docs/
  DOWNLOAD.md       how to get data out of the SRA fast, with measurements
  ARCHITECTURE.md   how the read-architecture call works
  INTERFACES.md     the internal module contract
```

## Credits

The QC and read-architecture methods are ported from the `read_architecture`
pipeline developed for PRICE2, where they were validated against 186 curated human
Ribo-seq studies (91 % exact full-architecture recovery). The ports are verified
bit-for-bit against those reference implementations.
