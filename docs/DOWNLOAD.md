# Downloading from the SRA

## Full runs

RiboMine fetches a full run from ENA over HTTPS with `aria2c`, using
`download.connections` parallel ranged connections (default 8). ENA serves the
submitter's gzipped FASTQ, so nothing has to be converted, and the file is checked
against ENA's md5. Two fallback routes cover the runs ENA does not have; each is
tried only when the one before it cannot serve the accession:

| | route | used when |
|---|---|---|
| 1 | `ena_https` | ENA has a FASTQ for the run |
| 2 | `aws_odp` | no ENA FASTQ: `.sra` from SRA's Open Data bucket on S3 (anonymous HTTPS, `aria2c`), converted with `fasterq-dump` |
| 3 | `prefetch` | not on S3 either: SRA toolkit `prefetch` + `fasterq-dump` |

Every route ends in one gzipped FASTQ (R1 of a paired deposit). ENA with parallel
connections comes first because throughput is limited per connection, not by the
link. Measured from a German university network, July 2026:

| source | 1 connection | 16 connections |
|---|---|---|
| ENA FASTQ | 11.9 MB/s | 50.1 MB/s |
| AWS Open Data `.sra` | ~3 MB/s | ~18–26 MB/s |
| `prefetch` (always one stream) | 2.3 MB/s | – |

End to end, a 500 MB run took 18.5 s from ENA with 16 connections and 134.7 s with
`prefetch` + `fasterq-dump`, which also has to convert the `.sra` and gzip the
output. ENA coverage of published data is nearly complete: of 3,276 published
human Ribo-seq runs, one had no ENA FASTQ. ENA lags NCBI by about 4–6 weeks, and
controlled-access (dbGaP) runs are never mirrored.

## The QC read sample

The QC stage does not download the run. It streams the first `qc.scan_reads`
reads (default 1,000,000) of the ENA FASTQ, reservoir-samples `qc.sample_reads`
of them (default 200,000) and closes the connection: a few seconds and tens of MB
per run. A broken stream is re-opened and sampled again from the start. Runs
without an ENA FASTQ are sampled with `fastq-dump -X n` (the first `n` spots).

ENA stores reads in sequencer order, so the prefix is representative unless the
deposit is sorted or collapsed. RiboMine warns when the stream looks sorted by
sequence or length; `qc.scan_reads: 0` then samples the whole run.

## Configuration

| key | default | meaning |
|---|---|---|
| `download.connections` | 8 | parallel connections per file (`aria2c` allows at most 16) |
| `download.max_retries` | 4 | attempts per route before the next route is tried |
| `download.retry_backoff_s` | 5 | wait before the first retry, in seconds; doubles each time |
| `download.tmpdir` | `null` | download and conversion scratch; `null` = `<workdir>/tmp` |

## Pitfalls

* Without `aria2c`, a single `curl` stream is used, which is several times slower.
* `download.connections` applies per run. When many runs download at once ENA may
  refuse requests; lower the value (see `docs/SLURM.md`).
* `prefetch --max-size` defaults to 20 GB and exits 0 when it skips a larger run.
  RiboMine passes `--max-size u` and checks that the `.sra` exists.
* SRA Lite (`.sralite`) files carry placeholder quality scores. The `aws_odp`
  route fetches the full-quality `.sra`, not SRA Lite.
* `fasterq-dump` needs about 10× the `.sra` size in scratch space. RiboMine uses
  `/dev/shm` when the run fits there, otherwise `download.tmpdir`.
