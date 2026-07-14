# Getting data out of the SRA, fast

Downloading is the slowest part of mining the SRA, and the route matters more than
your bandwidth does. RiboMine does not make you choose one, and does not offer a
knob to change it: one route is several times faster than the others everywhere we
measured, so it is simply *the* route. This is which one and why — all numbers
measured on the machine RiboMine was developed on (a German university host, July
2026), not copied from a blog post.

## The short version

**RiboMine downloads from ENA over HTTPS with 16 parallel connections.** That is the
whole policy.

The catch is that ENA does not mirror *every* run, and a run we cannot fetch is a run
we cannot mine — so two more routes sit behind it, tried in order when the one above
cannot serve the accession at all:

| | route | when it runs |
|---|---|---|
| 1 | **`ena_https`** | always, when ENA has the run — ENA's own gzipped FASTQ |
| 2 | `aws_odp` | ENA has no mirror. SRA's Open Data bucket on S3 (anonymous HTTPS), then `fasterq-dump` |
| 3 | `prefetch` | neither of the above. The SRA toolkit: always works, never fastest |

That is a fallback for **availability, not a preference for speed** — 2 and 3 are the
slow routes, and they run only when route 1 has nothing to give. Do not reorder the
chain without re-measuring; the ordering is what the rest of this document is about.

## Measured

Benchmark of one 500 MB run (SRR618773), whole file, cold:

| route | throughput | wall |
|---|---|---|
| **ENA + `aria2c -x16`** | **27.0 MB/s** | 18.5 s |
| `prefetch` + `fasterq-dump` | 8.4 MB/s | 134.7 s |

and on a second run (SRR12109318, 81 MB), isolating each leg:

| leg | throughput |
|---|---|
| ENA, `aria2c -x16` | **41 MB/s** |
| ENA, single-stream `curl` | 14.9 MB/s |
| AWS ODP `.sra`, `aria2c -x16` | 17.9 MB/s |
| AWS ODP `.sra`, single-stream `curl` | 3.8 MB/s |
| `prefetch` (single stream, by design) | 2.3 MB/s |

**End to end, ENA + aria2c is ~7–18× faster than the `prefetch` + `fasterq-dump`
route that most pipelines default to.**

### Why: the bottleneck is per-connection, not bandwidth

Sustained transfer of a multi-GB file:

| connections | ENA | AWS ODP (us-east-1) |
|---|---|---|
| 1 | 11.9 MB/s | ~3 MB/s |
| 16 | **50.1 MB/s** | ~18–26 MB/s |

There is a hard **~12 MB/s ceiling per TCP stream at ENA**, and ~3 MB/s from S3
to Europe (that one is transatlantic round-trip time against TCP's window, not an
S3 limit). Neither is your link. **Parallel connections are the entire game** —
which is why `aria2c -x16` is not a nice-to-have and `prefetch`, which opens one
stream, cannot be rescued.

Don't over-parallelise *across* runs either: 4 runs × 16 connections saturates the
link, and more concurrency just splits the same ceiling.

### Why ENA first

ENA serves **the FASTQ the submitter uploaded**. There is no `.sra` container to
convert, so the bytes on the wire are the bytes you want, already gzipped.
`fasterq-dump` by contrast emits an *uncompressed* FASTQ — about 6× larger — which
you then have to gzip yourself.

The obvious worry is coverage, and it is real but narrow. Measured:

| cohort | runs | no ENA FASTQ |
|---|---|---|
| **human Ribo-seq (published)** | 3,276 | **0.03 %** (1 run) |
| runs released ≥ 6 weeks ago | 10k–64k | 0.8 – 3 % |
| runs released in the last ~4 weeks | 18k–47k | **74 – 79 %** |

**ENA's mirroring lags NCBI by about 4–6 weeks.** For published Ribo-seq — which
is RiboMine's whole use case — coverage is effectively total. The permanent gaps
are dbGaP/controlled-access data (100 % absent from ENA) and withdrawn records;
both fall through to the next route, which is exactly what the fallback chain is
for.

We verified the two sources agree: ENA's FASTQ and the NCBI-derived FASTQ for the
same run are **identical** (same reads, same sequences — only the header text
differs).

## What we deliberately do *not* do

**Aspera.** NCBI's is **dead** — their resolver now answers `Protocol 'fasp' is
retired. Only 'https' is supported`. ENA's is still alive, but its own docs
suggest ~37 MB/s, which is *slower* than the 50 MB/s plain `aria2c -x16` already
gets from the same servers, and it drags in a licensed client. On a short, clean
European path to Hinxton there is nothing for it to win back.

**`prefetch` as the primary route.** It resolves to the same S3 URL `aws_odp`
uses, then fetches it over one throttled stream. It is strictly dominated. It stays
in the chain only because it is the one route that always exists.

**GCP.** `gs://sra-pub-run-odp` does not exist; the Google SRA buckets are
requester-pays and hand back SRA Lite by default. AWS ODP is free, anonymous, and
full-quality.

## Landmines this code guards against

* **`prefetch --max-size` defaults to 20 GB and exits 0 when it skips an oversized
  run.** No files, no error. Any script that trusts the exit code silently produces
  nothing. RiboMine passes `--max-size u` *and* checks the file exists.
* **SRA Lite (`.sralite`) has fake quality scores** — a flat Q30/Q3 per read. It is
  ~60 % smaller and Google hands it back by default. It would quietly poison any
  quality-aware step. The AWS ODP copy we fetch is the full-quality one.
* **`fasterq-dump -e` does not scale**: `-e1` 8.6 s → `-e4` 6.7 s → `-e16` 6.6 s.
  It is I/O-bound writing the uncompressed FASTQ. RiboMine caps it at 6 threads and
  gives the rest to the other samples in the pool.
* **`fasterq-dump` needs ~10× the `.sra` size in scratch.** RiboMine puts it on
  `/dev/shm` when the run fits there, and falls back to the workdir when it does not
  (filling `/dev/shm` takes the machine down).
* **ENA's `sra_ftp` / `sra_bytes` fields are declared but never populated.** ENA
  holds no `.sra` at all. Don't build on them.
* **Downloads get truncated.** ENA returns the md5 in the same call as the URL, so
  RiboMine checks it. A short file otherwise surfaces as a strange read count ten
  steps later, where nobody thinks to blame the transfer.

## The read sample: don't download the run at all

The QC and architecture stages need a *sample* of reads, not the run — the
architecture is a property of every read, and periodicity is measurable on a few
tens of thousands. gzip is a stream format, so **a byte-prefix of the ENA `.gz`
decompresses to a read-prefix**: RiboMine opens the URL, reads the first
`qc.scan_reads` reads, reservoir-samples `qc.sample_reads` of them, and drops the
connection.

Cost per dataset: **a few seconds and tens of MB**, against gigabytes and minutes
for the run. Screening 500 candidates therefore costs about as much as downloading
one of them.

Two things to know about that sample:

* It is a **uniform random sample of the reads it scanned** (the kept records are
  shuffled), but only over the run's first `scan_reads`. SRA/ENA store FASTQ in
  spot order — the order reads came off the sequencer, random with respect to
  content — so the prefix is representative.
* The exception is a **sorted or collapsed deposit**, where the prefix is a biased
  slice. RiboMine detects that (leading k-mers or lengths monotone across the
  stream) and warns; `qc.scan_reads: 0` then reservoir-samples the whole run.

Runs with no ENA mirror fall back to `fastq-dump -X n`, which does stream and stop
early — correct, just slower.

> **Do not use `fasterq-dump --row-limit` for subsampling.** It is applied
> *per thread*, at each thread's slice offset, so you get `row_limit × threads`
> reads scattered across the run rather than a prefix. This is undocumented.
