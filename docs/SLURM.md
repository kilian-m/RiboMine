# RiboMine on LRZ CoolMUC-4

Mining the archive is a multi-day, multi-terabyte job, so it runs as a SLURM chain
rather than as one `ribomine run`. Three jobs:

```
prep            cm4_tiny, 1 node    query the SRA -> candidates.tsv
                                    build the GTF + contaminant indexes (once)
                                    split the cohort: one accession list + one
                                    workdir per node
   |  afterok
   v
run             cm4_std, ARRAY OF 2 x 4 NODES = 896 cores
                                    each of the 8 nodes runs one `ribomine run`
                                    over its own shard
   |  afterany
   v
merge           cm4_tiny, 1 node    symlink the 8 workdirs into one, write the
                                    cohort tables, and say in logs/status.txt
                                    whether anything is left
```

```bash
# from an LRZ login node
slurm/master.sh config/config_lrz.json          # the whole chain
slurm/master.sh config/config_lrz.json run      # another round (see "The loop")
DRY_RUN=1 slurm/master.sh config/config_lrz.json    # print the sbatch lines, submit nothing
```

---

## What you have to do first

**1. Get RiboMine onto LRZ and build the environment.** In `$HOME`
(`/path/to/`), not on DSS — DSS is for data.

```bash
git clone git@github.com:kilian-m/RiboMine.git && cd RiboMine
conda env create -f environment.yml
conda activate ribomine
pip install -e .
ribomine --version          # it should print
```

**2. Edit `config/config_lrz.json`.** Four things, and `prep` checks three of them
for you and dies in minutes rather than letting a 24 h job die on a typo:

| key | what to check |
|---|---|
| `project.workdir` | a DSS path with **several TB** free. The BAMs are the small part; the transient FASTQs are the big one (see *Disk*). |
| `reference.genome_fasta`, `reference.gtf` | exist, and are the pair the STAR index was built from |
| `reference.star_index` | exists — **and was built with `--sjdbGTFfile`**, or there is no gene-count matrix. `prep` says which. |
| `_slurm.mail_user` | currently `` |

The STAR-index requirement is not a nicety: STAR counts reads into genes *while* it
aligns them, and it will not insert junctions on the fly against a shared-memory
genome — which is what makes 28 concurrent samples per node affordable. So the
annotation is either in the index or there is no `counts/gene_counts.tsv`.
Everything else still works.

**3. Check that a cm4 compute node can reach ENA.** This is the one thing I could
not verify from here, and the whole pipeline is built on it. `price2-expansive`
downloads from `rdp.ucc.ie` on cm4_tiny compute nodes, so outbound HTTPS works —
but confirm the hosts RiboMine actually uses:

```bash
salloc -M cm4 -p cm4_inter -N 1 -t 00:10:00
curl -sI https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR12285169\&result=read_run\&fields=fastq_ftp | head -1
curl -sI https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra\&term=riboseq | head -1
exit
```
Both should be `HTTP/... 200`. If they are not, the compute nodes are behind a
proxy and `http_proxy`/`https_proxy` need exporting in `slurm/node_worker.sh` and
`slurm/prep.sh`.

**4. Do the cheap pass first.** Set `pipeline.end` to `"qc"` and submit. QC never
downloads a run — it streams a 200k-read sample — so screening the whole archive
costs hours instead of days, and `merged/qc/qc_summary.tsv` then tells you how many
of the hits are actually Ribo-seq **before** you spend a cohort's worth of bandwidth
on them. Then set it back to `"bam"` and `slurm/master.sh <config> run`: every QC
verdict is reused, nothing is recomputed.

This matters more than it sounds. The default query is deliberately a wide net
(`query.terms` is the recall net; precision comes later), so it can return
*thousands* of candidates, most of them Ribo-Zero RNA-seq. At 1–10 GB each, an
unscreened `end: "bam"` run is tens of terabytes. Either do the QC pass first, or
set `query.max_runs`.

---

## The loop

The run has a **24 h ceiling** and a real cohort does not fit in it. That is
expected, not a failure:

1. `merge` runs `afterany`, so it runs however the array ended — cleanly, out of
   wall time, or with a node dead. It works out what finished from the per-sample
   JSONs, never from exit codes, so those three cases need no telling apart.
2. It writes `<workdir>/logs/status.txt`: **COMPLETE**, or **RESUBMIT NEEDED** with
   the command.
3. You run `slurm/master.sh config/config_lrz.json run` from a login node. Finished
   runs are skipped, so each round is shorter than the last.

**Why you and not a job:** a cm4 compute node is not permitted to `sbatch` (LRZ
rejects it with *"Access/permission denied"*), so no job in the chain can launch the
next one — and a long-running driver on a login node is also not allowed. This is
the same hard-won conclusion `price2-expansive` reached; do not try to automate it
back.

---

## The shape, and why

**896 cores is an array of two 4-node jobs, not one 8-node job.** `cm4_std` caps a
single job at **4 nodes / 448 cores**, but runs **2 jobs at once**. So:

```
_slurm.jobs = 2  x  _slurm.nodes_per_job = 4  x  112 cores  =  896
```

and each node's shard is its global index across the array,
`SLURM_ARRAY_TASK_ID * 4 + SLURM_PROCID` → 0…7.

**One `ribomine run` per node — not two, and not one per core.** Each node loads the
~28 GB STAR index into its own shared memory once, and all 28 of its concurrent
samples attach to that single copy; without that a node would need 28 × 28 GB. A
second RiboMine on the same node would race it: whichever finished first would call
`STAR --genomeLoad Remove` and pull the index out from under the other.

**28 samples × 8 threads = 224 on a 112-core node, deliberately.** That is the node's
*logical* core count (112 physical + hyperthreads), so it is twice as many samples in
flight as there are real cores — and it should be: this stage spends much of its time
blocked on an ENA download, and a sample waiting on the network is not using a core.
What it costs is memory: ~28 GB for the node's one shared STAR genome plus roughly
8 GB per sample in flight, so budget ~250 GB of the node's 480 GB. `project.jobs` is
the number to cut if a node starts swapping. `master.sh` only warns past 2×, where
there are not even hardware threads left to hold the samples.

`node_worker.sh` also removes any **leaked** STAR segment on the way *in*, not only
on the way out. A segment leaked by a wall-time SIGKILL holds 28 GB on that node and
makes the next job's load fail; clearing it at startup is the recovery path that does
not depend on the previous job having exited politely. It is safe because cm4_std
allocates nodes exclusively.

**Per-node workdirs (`shards/work_NN/`), merged afterwards.** `ribomine run` writes
its cohort tables (`qc_summary.tsv`, `failed.tsv`, …) at fixed paths under its
workdir, so eight nodes sharing one workdir would each overwrite the other seven.
They get one each; `merge` then symlinks all eight into `merged/` and calls
RiboMine's *own* report writers over the whole cohort — so the tables come from the
code that would have produced them on one machine, not from a bespoke concatenator
that would drift from the columns. What the nodes *do* share, read-only, is
`refs/` (the GTF index, the bowtie2 contaminant index) and `meta/candidates.tsv`,
symlinked in by `prep`.

**The split is pinned.** Re-running `prep` after the query has grown does not
reshuffle: an accession that already has a shard keeps it, and only new ones are
packed. Moving a finished accession to another node would strand its results in a
workdir nobody reads any more, and it would be downloaded and mapped a second time.
Packing is longest-first by `read_count` — a run is 1–10 GB, and the node that draws
the deep ones sets the wall time.

---

## Disk

`keep.bam` is the only thing on by default, so what *survives* is small. What
*passes through* is not: each run is downloaded (1–10 GB), trimmed, and
contaminant-filtered, and all three live on DSS at once before the intermediates are
deleted. At 28 concurrent samples × 8 nodes — 224 runs in flight — that is roughly
**2–3 TB of transient DSS at peak**, plus the BAMs (~0.5–2 GB each) that stay.

The download's temp directory stays on the workdir's filesystem on purpose — it is
not moved to the node-local NVMe. The ENA route finishes with `os.replace(tmp,
out.fastq.gz)`, which is a rename and **fails across devices**; keeping it on DSS
makes that rename free, whereas `/tmp` → DSS would be a multi-GB copy of every run.
(Node-local `$TMPDIR` *is* used, for tool scratch.)

---

## Things to watch

**ENA throttling.** `download.connections` is *per run*, and 28 runs download at once
on each of 8 nodes — so it is multiplied by 224. It is set to **4** here, not
RiboMine's own default of 8: that is still ~900 sockets against one host at peak,
and 8 would be ~1,800. It costs less than it looks. The per-run connection count
exists to beat ENA's ~12 MB/s per-stream cap, and with 224 runs already in flight the
link saturates long before any single one of them does — the concurrency that matters
here is across runs, not within one.

If downloads start failing (they show up in `failed.tsv` and `mapping_summary.tsv`),
this is still the first knob to turn — downwards.

The cohort is bandwidth-bound long before it is core-bound, which is also the honest
answer to "do I need 896 cores?" — for `end: "bam"`, probably not; for `end: "qc"`,
the cores are the point.

**A run that fails every round.** `status.txt` separates *still to do* from *tried
and failed*, and lists the latter. They do not retry themselves into eternity —
read `merged/failed.tsv` and decide.

**Files**

| | |
|---|---|
| `slurm/master.sh` | the submitter. Login node only. |
| `slurm/prep.sh` | validate + `ribomine setup` + query + split |
| `slurm/run.sh` | the cm4_std array |
| `slurm/node_worker.sh` | one node: its shard, and its STAR shared memory |
| `slurm/merge.sh` | the cohort tables + `status.txt` |
| `scripts/prepare_run.py` | the split (pinned, LPT-packed by read count) |
| `scripts/merge_results.py` | the merged view, the tables, the verdict |
| `config/config_lrz.json` | one file: RiboMine's settings *and* the cluster shape (`_slurm`, which RiboMine ignores) |
