# How the read-architecture call works

> Ported into RiboMine from the `read_architecture` pipeline built for PRICE2, where the
> method was validated against 186 curated human Ribo-seq studies (91% exact full-architecture
> recovery, 95% on the 5' trim boundary). Module names below are RiboMine's.

This is the deep-dive companion to `README.md` (which is the quick start). It
explains what each stage computes, **which quantity every decision is made on**,
and how to read those quantities off the diagnostic figure the pipeline now
draws for every sample.

---

## 1. The problem

A ribo-seq read is a footprint wrapped in library scaffolding. Reading 5'→3',
any of these may be present or absent:

```
[5' UMI][RT nt][============ footprint ============][3' UMI][barcode][adapter]
   random  A/T          genomic, variable-length        random   fixed   fixed
```

plus poly(A) tails and template-switch (G/C-rich) motifs. We are given only an
accession or a FASTQ and must report the layout — lengths and identities — with
**no prior knowledge of the protocol**, and refuse when the reads cannot support
an answer.

The output is not the *names* of the parts but what a processing pipeline needs:
how many nt to trim from each end, which stretches are random (so must be used
for UMI deduplication), and which adapter to remove.

---

## 2. The one idea everything rests on

We do **not** trim before aligning. We align the raw read with
`STAR --alignEndsType Local`, so everything non-genomic — 5' UMI, RT additions,
3' UMI, barcode, adapter — is pushed into the **soft clips** instead of
preventing the alignment. The architecture is then read back out of the clips.

The catch: STAR's clip boundary is **not** where the footprint starts. If the
last base of a 5' UMI happens to match the genome base next to the footprint
(probability 1/4), STAR extends the alignment by one and the clip shrinks.
**Clip lengths leak.**

What does *not* leak is the genomic coordinate the alignment implies for **read
position 0**. Extending the alignment leftwards by one decrements *both* `pos`
and `clip5`, so `pos - clip5` is invariant. Anchor on that, and every read
position `p` maps to a fixed genome base. Then ask the only question that matters:

> Does read base `p` match the genome base implied for it?

* Footprint bases are genomic → they match ~95–99 % of the time.
* UMI / adapter / barcode bases are not genomic → they match ~25 % by chance.

**The transition between those two regimes is the architecture.** Everything
downstream is a way of locating that transition precisely and naming what sits
on either side.

---

## 3. The four stages

| stage | script | in → out |
|---|---|---|
| 1. sample | `sra/download.py` | accession/FASTQ → 200k reads (reservoir-sampled, streamed off ENA; nothing downloaded in full) — see note below |
| 1b. de-contaminate | `qc/contaminants.py` | FASTQ → filtered FASTQ (rRNA/tRNA/etc + adapter-dimers + low-complexity removed; §14) |
| 2. align | `process/star.py` | filtered FASTQ → BAM (`STAR --alignEndsType Local`, permissive length filters, no trimming) |
| 3. profile | `qc/profile.py` | BAM + genome → positional statistics JSON (§4) |
| 4. QC | `qc/verdict.py` | BAM + annotation → is-it-ribo-seq verdict + quality metrics (§12) |
| 5. call | `arch/infer.py` | profile JSON → the architecture call (§5–7) |
| 6. plot | `arch/plot.py`, `qc/plot.py` | JSONs → annotated figures (§8, §12) |
| 7. trim | `arch/trim.py` | FASTQ + architecture → trimmed FASTQ (UMIs → read name) + trim stats |
| 8. report | `reports.py` | all JSONs (+ re-aligned trimmed reads) → the metrics list (§13) |

Stages 3, 5, 6 read only the profile JSON, never the BAM, so the caller and the
architecture plot iterate in seconds without re-aligning. The QC stage (§12)
reads the BAM against a cached annotation index. `architecture.sh` runs the whole
pipeline for one sample; `run_batch.sh` → `reprofile_all.sh` → `evaluate.py` runs
the architecture benchmark over the cohort.

> **What the sample is, precisely.** Reservoir sampling gives a *uniform* random
> sample of the reads it scans, and the kept reads are shuffled, so ordering
> *within the scanned window* is gone. But only the first `--scan` (default 1M)
> reads are streamed, so a run larger than that is sampled uniformly over its
> **first 1M reads, not the whole file**. ENA/SRA fastq are stored in spot order
> (random w.r.t. content), so that prefix is representative — the exception is a
> **sorted or collapsed** deposit, where the prefix is a biased slice. Those trip
> a warning (reliable for sequence-sorting, best-effort for length-sorting); pass
> `--scan 0` to reservoir-sample the entire run. Because the architecture is a
> per-read-invariant property of the library prep, the *calls* are robust to this
> even when the prefix is skewed — the exposure is limited to the composition /
> plateau / footprint-length statistics.

---

## 4. The quantities (`qc/profile.py`)

Everything is descriptive here — no thresholds, no calling. For each read the
profiler keeps the base string, a per-position **match mask** against the implied
genome base, and the two soft-clip boundaries, then aggregates:

**5'-anchored** (read position 0 = first base; exact, no leakage):
* `p5_match[p]` — fraction of reads whose base `p` matches the genome. This is
  essentially the **CDF of the 5'-clip length**, lifted by a 1/4 chance floor:
  `P(p inside aligned block) + P(p clipped)·¼`.
* `p5_comp[p]` — the A/C/G/T composition at position `p`.
* `clip5_hist` — the 5'-clip-length histogram, i.e. the **survival function**
  that `p5_match` integrates.

**3'-anchored.** Read position is meaningless on the 3' side (footprint length
varies), and the alignment *end* is worse than useless — local alignment stops
at a mismatch, so the base just past it disagrees with the genome *by
construction* (match rate 0.00). The two honest anchors are:
* the **adapter start** — `adap_match` / `adap_comp`, indexed so `A0` is the
  first adapter base and `A-1, A-2, …` walk left toward the footprint. The
  adapter is discovered from a panel (TruSeq, small-RNA RA3, Ingolia linker,
  Nextera, poly-A/G) *or* de-novo from the most common k-mer in the non-genomic
  tail, and — for libraries that read through a ligated linker into the TruSeq
  primer — the one ligated **nearest the insert** (smallest median start) wins.
* the **read's own 3' end** — `t3_match` / `t3_comp`, the fallback when no
  adapter is present (already-trimmed deposits).

**Shape histograms:**
* `footprint_len_hist` — genuine footprints spread across ~26–34 nt; a single
  sharp spike means the "footprints" are a fixed-length artefact (adapter dimer).
* `fpend_to_adapter_gap_hist` — sharp mode = fixed-length 3' construct; broad =
  a variable poly(A) tail sits in between.
* `top5p_locus_frac` — the fraction of reads at the single most common 5'
  genomic coordinate (leakage-invariant). Genuine footprints spread across
  thousands of loci (< ~2 % at any one); a majority at one coordinate is a
  single over-represented species, not footprints.

Reads with indels or splice junctions break the position→genome map and are
skipped (a few %); adapter-dimers and poly(A) reads map to genomic homopolymers
and are dropped by a low-complexity filter on the aligned core.

---

## 5. Reading the 5' side (`call_5prime`)

Two derived quantities do the work, both relative to the sample's own **genomic
plateau** (the median deep-position match rate — a noisy library plateaus at
0.79 where a clean one reaches 0.97, so all thresholds are scaled to it, never
absolute):

* **penetrance**`[p] = (plateau − p5_match[p]) / (plateau − ¼)`, clipped to
  [0,1]. This is the survival function of the per-read 5'-construct length: a
  fixed *L*-nt element sits at penetrance ≈1 for *L* positions then steps to 0;
  an element present in only a fraction *f* of molecules sits at *f*.
* **composition** — entropy and A/T vs G/C bias per position.

The rules, footprint-first:

| element | signature | fate |
|---|---|---|
| **footprint start** | first position that is, and stays, genomic (penetrance < `FOOTPRINT_PEN`) | the boundary; everything left is trimmed |
| **RT untemplated nt** | 1–2 nt against the footprint, at *partial* penetrance **or** A/T-biased | **kept** (it is enzymatic, part of the molecule) |
| **template-switch / linker** | full penetrance, **G/C-biased** (≥ `GC_TEMPLATE`) | trim |
| **UMI** | full penetrance, balanced, high entropy | trim **+ dedup** |
| **barcode** | full penetrance, **fixed base** (low entropy) | trim |

The key discriminations: a **step** in `p5_match` (not a level) marks the UMI
length — a smooth ramp is just a tail of mis-mapped reads and must not be read as
a long UMI. A UMI base is uniform; an **RT base is A/T-biased and partial**,
which is how a non-templated addition is told apart from a UMI base even when the
two abut. **No protocol uses a 1-nt UMI**, so a lone non-genomic 5' base is an
RT addition, not a UMI.

---

## 6. Reading the 3' side (`call_3prime`)

Walk **left from the adapter start** through `adap_match`, nearest-first:

* Stop at the footprint — a base that is genomic (match ≥ `thr_genomic`) **and
  variable** (entropy > const). The "and variable" matters: a *constant* linker
  base can match the genome well above chance simply because one fixed base is
  compared against a biased genomic neighbourhood, so a high match rate alone
  must not end the walk.
* A real UMI/barcode base is *clearly* non-genomic (match ≈ chance). A position
  with only an *intermediate* match is the footprint's own 3' edge, depressed by
  alignment-edge leakage — not a construct base. If **nothing** between footprint
  and adapter is clearly non-genomic, the walk only caught a footprint edge and
  the 3' construct is dropped (no fabricated 1-nt UMI).
* A **homopolymer** run (poly-A) is constant in composition but variable in
  length, so it is a tail, not a fixed barcode, and is peeled off separately.

What remains between footprint and adapter is segmented by entropy into
`random` (UMI), `degenerate`, and `const` (barcode) blocks, in insert-first
order.

---

## 7. The four functional categories

The named fields (`umi5`, `barcode3`, …) are re-cast into what actually drives
processing — **origin and fate** — because the residual "±1 UMI-vs-spacer"
ambiguity lives *inside* one category and therefore changes none of the derived
quantities:

| category | penetrance | composition | genome match | fate |
|---|---|---|---|---|
| footprint | ~0 | genomic | high | **keep** — the data |
| random-templated (UMI / spacer) | ~1 | balanced | chance | trim **+ dedup** |
| fixed-templated (barcode, linker, adapter) | ~1 | low entropy | chance | trim |
| enzymatic (RT nt / template-switch) | partial (RT) or ~1 (TS) | A/T (RT), G/C (TS) | chance | RT nt **kept**, TS trimmed |

`infer()` emits a `functional` block with the 5'→3' segment list and the two
load-bearing numbers: **`trim_5p`** (nt to clip from the 5' end to reach the
footprint) and **`dedup_umi_len`** (total random-templated content, the mask a
deduplicator needs).

---

## 8. The diagnostic figure (`arch/plot.py`)

Every quantity a decision is made on is plotted, with the inferred segmentation
and the exact thresholds drawn on top — so the figure is a **visual audit of the
call**, not a separate summary. The title reads out the inferred architecture
5'→3' — e.g. `5'-[UMI,2nt]-[footprint,~23nt]-[UMI,5nt]-[barcode,GATCA]-[TruSeq]-3'`
— with the trim / adapter / dedup quantities on the line below. Below is
`SRR24652840`, a McGlincy-Ingolia library:

![example architecture figure](results/plots/SRR24652840.arch.png)

There are two legends: the **nucleotide colours** (A/C/G/T) for the composition
panels B and E, and the **segment categories** for the shaded spans in the
match panels A and D. (The two palettes are similar hues but never share a
panel.) An RT base shows as `[RT,<penetrance>%]` — the fraction of molecules
carrying it — because it is added to only some reads.

| panel | quantity | what to look for |
|---|---|---|
| **A** 5' genome-match & penetrance | `p5_match`, penetrance | the **step** up to the plateau — its position is the 5' UMI length; the shaded span is the trimmed construct, the green line the footprint start |
| **B** 5' composition & entropy | `p5_comp` | balanced+high-entropy = UMI; flat/low-entropy = barcode; A/T spike = RT nt; G/C = template-switch |
| **C** 5' clip-length histogram | `clip5_hist` | penetrance is the survival function of this; a fixed UMI is a spike at its length |
| **D** 3' genome-match, adapter-anchored | `adap_match` | the drop from footprint (green) down through the UMI/barcode construct (blue/orange) to the adapter (grey `A0`); dashed lines are the genomic and UMI thresholds |
| **E** 3' composition & entropy | `adap_comp` | the construct's composition; a sharp entropy cliff at the footprint↔adapter boundary |
| **F** footprint-length distribution | `footprint_len_hist` | genuine footprints spread ~26–34 nt; a single spike above the red **artefact gate** (0.85) means it is refused |

Run it standalone on any profile:

```bash
python bin/plot_architecture.py data/profiles/SRR24652840.json -o fig.png
python bin/plot_architecture.py data/profiles/ -d results/plots     # whole cohort
```

`bin/show_profile.py` prints the same numbers as a text worksheet if you prefer
a terminal view.

---

## 9. Refusing to answer

`undetermined` is returned — never a guess — when the profile carries no
readable boundary: too few usable alignments (< 5 000), no genomic plateau
(deep match rate < 0.55, i.e. alignments too noisy to read), a footprint-length
spike above the artefact gate (adapter dimers / fixed contaminants), **more than
half the reads at a single 5' coordinate** (a single over-represented species —
an adapter dimer that aligns to an adapter-like locus, a spike-in, a contaminant
— where genuine footprints spread across thousands of loci; the busiest genuine
library here is 38 %), only a homopolymer 3' anchor (a poly-A library whose real
adapter is invisible — the tail is still named in the reason), or a 3' tail that
cannot be segmented. On the cohort this fires on ~10 % of samples, all genuinely
unreadable. The refused reason is written into the figure title.

---

## 12. Is it ribo-seq at all? (`qc/verdict.py`)

Before calling architecture, the pipeline asks whether the library is even
ribo-seq, using the signatures RiboseQC / riboWaltz use. It reads the genome BAM
against a cached annotation index (`qc/annotation.py` turns the ~1 GB Ensembl
GTF into a 22 MB pickle of per-chromosome interval arrays — CDS with reading
frame, UTRs, ncRNA exons, gene spans, start/stop codons — loaded in a fraction of
a second and turned into NCLS interval indexes).

| signature | what it measures | ribo-seq |
|---|---|---|
| **read-length peak** | footprint length distribution | tight peak ~28–32 nt |
| **3-nt periodicity** | frame of the read 5′ end within CDS (Ensembl phase), each length aligned to its own dominant frame; reported as the in-frame fraction and the **TVD of the frame distribution to uniform** (0 = flat, 0.67 = perfect) | in-frame ≫ ⅓ |
| **CDS enrichment** | fraction of genic reads in CDS (vs UTR/intron) | ≫ RNA-seq |
| **start/stop metagene** | 5′-end density around annotated codons | periodic comb, P-site offset ~12 |
| **region composition** | CDS / 5′UTR / 3′UTR / ncRNA / intron / intergenic / mito | CDS-dominated |
| **single 5′ locus** | `top5p_locus_frac` (§4) | ≪ 0.5 (else contaminant) |
| **unique-mapping rate** | fraction of reads that map uniquely to the genome (STAR log) | ≥ ~20 % (cohort median 42 %) |

**The verdict has three categories:**

* **`RIBO-SEQ`** — footprint-like, CDS-enriched, periodic (elongating ribosomes).
* **`TI-SEQ`** — translation-initiation sequencing (harringtonine / LTM / QTI- /
  GTI-seq). The drugs freeze ribosomes at start codons and elongation runs off, so
  the start-codon metagene peak **towers over a near-empty CDS body**. Measured as
  the *start-codon enrichment ratio* (peak ÷ downstream in-frame body); above
  `TISEQ_RATIO_MIN` (40) it is TI-seq (GTI-seq `SRR618773` = 120×; elongating
  ribo-seq tops out ~30×). 1/187 cohort samples.
* **`NOT RIBO-SEQ or LOW QUALITY`** — everything else, *with the reason stated*:
  a single-locus contaminant, **low unique-mapping** (< `MIN_UNIQUE_FRAC` = 0.15,
  i.e. dominated by non-genomic / multimapping junk — `SRR1507058`: 9 % unique),
  no periodicity, or atypical length/CDS. Low mapping quality folds in here even
  when the mappable minority is periodic, because the library is not usable.

**Projected usable reads** (verdict panel, section F): the reads that survive all
filters and map (`n_reads_scored`) scaled from the 200 k sample to the whole run —
`(n_reads_scored / n_sampled) × total_run_reads`, with the run's total read count
from ENA (single-sample) or `meta/ena_runs.tsv` (cohort). It answers "how many
usable footprints does this dataset actually contain?" (`SRR1507058`: 4.6 M of
80.7 M — most of the run is junk).

**Periodicity is the defining hallmark** — the frame of the footprint 5′ end,
measured against the CDS reading frame, is what separates translating ribosomes
from RNA-seq/small-RNA. The Ensembl `frame` column gives the phase of every CDS
base directly (`phase = (pos − CDSstart − frame) mod 3`), so no per-transcript
walk is needed. The verdict is deliberately **inclusive** so ribo-seq variants
are not discarded: strong periodicity is decisive on its own (it holds for
TI-seq / QTI-seq, whose initiating ribosomes are still in-frame, and for short or
long footprints); failing that, footprint-like + CDS-enriched reads with any
periodicity (or strong CDS enrichment) pass. Only a contaminant-dominated library
(a majority at one locus) is refused outright. The QC is measured on the local
alignment because soft-clipping the 5′ construct / RT base gives a cleaner
footprint 5′ end — and hence sharper periodicity — than end-to-end alignment of
the trimmed reads. `qc/plot.py` renders all six signatures with the
verdict.

## 13. Trimming and the metrics list (`arch/trim.py`, `reports.py`)

`arch/trim.py` applies the inferred architecture to the raw reads: it removes
the `trim_5p` 5′ overhead and the 3′ adapter / barcode / poly(A) tail, **keeps**
the enzymatic RT base inside the footprint, and moves the random-templated UMI
content (5′ + 3′) into the read name (`_<UMI>`, umi_tools style) for
deduplication. Reads left shorter than `--min-len` are dropped.

The decisions here are RiboMine's own — no external trimmer can express them (the
barcode+adapter scaffold anchor, the poly(A) fallback when the scaffold has run off
the read's end, a 5′ UMI split around a barcode). But the *string matching* under
them is cutadapt's C aligner, imported as a library (`cutadapt.align.Aligner`) and
driven by the logic above; and the gzip goes through `xopen`/pigz rather than
Python's `gzip` module. Together those took the stage from 60k to 271k reads/s on
one core, without changing a single cut point (verified: 0 differences over 300k
reads spanning the absent / partial / error-bearing scaffold cases).

`reports.py` then assembles the metrics list, taking **periodicity / region
metrics from the local alignment** (sharpest) and **mapping metrics from the
trimmed reads re-aligned end-to-end** (the permissive local alignment inflates
multimapping):

```
  ribo-seq?                  YES
  architecture               5'-[rt_untemplated:1]-[footprint~29nt]-3'
  periodicity (in-frame)     76%          periodicity TVD→uniform   0.43
  CDS enrichment             73% of genic reads
  contaminants removed       rRNA/tRNA/etc 71%   low-complexity 0%   position pile-up 3%
  reads mapped to genome     — of the filtered reads —
    uniquely mapped 46%   multimapping 27%   unmapped 27%   (sum to 100)
  region composition         CDS 69%  5'UTR 10%  intron 8%  ...
```

---

## 14. Removing contaminants — no hard-coded adapters

Contaminants are removed in two places, and **no adapter sequences are hard-coded
anywhere**:

**Before mapping (`qc/contaminants.py`)** — sequence-based, using references:

| contaminant | filter |
|---|---|
| **rRNA / tRNA / snRNA / snoRNA / Mt** | `bowtie2 --very-sensitive-local` against a contaminant FASTA (a reference of contaminant *sequences*, not adapters; local so the still-present 5′ construct / 3′ adapter are soft-clipped and the footprint core is matched) |
| **homopolymer / low-complexity** (poly-A/G, simple repeats) | Shannon entropy / single-base dominance on the read |

**After mapping (`qc/pileups.py`)** — position-based, fully data-driven. This
replaces the old hard-coded adapter-dimer detection. Adapter / primer dimers and
fixed contaminants all **collapse onto a single genomic 5′ position** as one
identical, fixed-length molecule — e.g. the TruSeq adapter maps to
**chr15:56,885,929** in many libraries (40–72 % of reads, all 22–27 nt, all
identical, not in CDS). A 5′ position is dropped when it is **both**:

* over-represented — ≥ `--min-frac` (0.5 %) of reads and ≥ `--min-count` (30), and
* a single fixed length — its footprints are ≥ `--len-conc` (85 %) one length.

The second condition is the safeguard: a genuinely translated codon can also be
highly covered, but its ribosome footprints always **spread across lengths**
(~26–34 nt), so it is never removed (a real high-coverage peak with 62 % length
spread survives; the chr15 adapter with 97–100 % one length does not). This
catches `SRR6181542` (55 % at one position) that a sequence-based adapter filter
missed, and needs no adapter sequences at all.

…**or** when it is dominant on its own:

* ≥ `--max-frac` (10 %) of the **surviving** reads at one 5′ coordinate, whatever
  its length spread — applied repeatedly to a fixpoint, so a second pile cannot
  hide behind a larger one (`SRR10491340` has three, `SRR18113903` three).

The fixed-length test assumes the contaminant is one *molecule* of one *length*,
so it misses a pile whose **insert** varies: the reads then align at 21/22/23 nt
from the same 5′ base and look length-diverse. `SRR2096968` (an ARTseq library)
put **76 % of its reads on chr15:56,885,929** at length concentration 0.67 —
under the 85 % gate, so the whole pile survived and the sample was rejected as
*“53 % of reads at one locus”*. The same blind spot let a **19 %** pile at that
same adapter locus through in `SRR29327151`, and a tRNA fragment
(`Y:3,367,775`, aligning at both 25 and 29 nt, concentration 0.60) through in
three others.

**10 % is far above anything real.** Across the cohort the top 5′ position of a
library that passes ribo-seq QC holds a **median 0.5 %** of reads and **3.8 % at
p90**, while *every* pile above ~10 % is demonstrably junk — adapter dimers,
poly(T)+P7 primer dimers, adapter-homologous force-alignments (100 % of their
reads carry a mismatch), and tRNA/rRNA fragments that slipped past the bowtie
contaminant filter. No translated codon holds a tenth of a library.

Widening `--len-conc` to a ±1 nt window would catch these piles too, but it keys
on the wrong property — length spread rather than abundance. The dominance rule
fires on **11/187** samples; on the 165 that pass ribo-seq QC it removes a median
of **0.0 %**. See RESULTS.md.

Removing rRNA is also **why multimapping looked huge** (and why the percentages
didn't add up): each rRNA has hundreds of genomic copies, so every rRNA read maps
to "too many loci". It takes `SRR1630831` from **78 % → 27 %** multimapping. The
report shows uniquely-mapped + multimapping + unmapped, which sum to 100 %.

Requires a bowtie2 index built once from the contaminant FASTA
(`bowtie2-build <contaminants.fa> data/contaminants/human_contaminants`), and
runs bowtie2 `--very-sensitive-local`. One residual limitation: a structured-RNA
species that is **not in the contaminant FASTA** (e.g. a tRNA gene absent from the
reference) can slip through as a short-footprint read; those still show up in the
QC region composition as `ncRNA`. The three contaminant fractions are reported in
the QC figure's verdict panel and the metrics list.
