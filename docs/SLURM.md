# RiboMine on LRZ CoolMUC-4

A large cohort is processed as a chain of three SLURM jobs. The scripts target
CoolMUC-4; on another cluster, adapt the `#SBATCH` headers in `slurm/*.sh` and the
`_slurm` block of the config.

```
prep    cm4_tiny, 1 node    check the reference, build the GTF and contaminant
   |                        indexes, split the cohort into one accession list
   | afterok                and one workdir per node
run     cm4_std, array of `jobs` x `nodes_per_job` nodes
   |                        each node runs one `ribomine run` over its shard
   | afterany
merge   cm4_tiny, 1 node    link the per-node workdirs into merged/, write the
                            cohort tables and logs/status.txt
```

```bash
# on a login node, from the repository root
slurm/master.sh config/config_lrz.json            # prep -> run -> merge
slurm/master.sh config/config_lrz.json run        # another round: run -> merge
slurm/master.sh config/config_lrz.json prep       # query + prep only
slurm/master.sh config/config_lrz.json merge      # merge only
DRY_RUN=1 slurm/master.sh config/config_lrz.json  # print the sbatch lines, submit nothing
```

With `pipeline.start: "query"`, `master.sh` first runs `ribomine query` on the login
node (the compute nodes drop the long-idle connection of the ENA search), unless
`<workdir>/meta/candidates.tsv` exists. Delete that file to search again.

## Before the first submit

**1. Install** in `$HOME`. `master.sh` also uses the environment on the login node.

```bash
git clone https://github.com/kilian-m/RiboMine.git && cd RiboMine
conda env create -f environment.yml
conda activate ribomine
pip install -e .
```

**2. Edit the config:** `project.workdir` (e.g. `<your-dss-path>/ribomine/all`), the
three `reference` paths and `_slurm.mail_user` (`""` for no mail). `prep` stops if a
reference path is missing, and warns if the STAR index was built without
`--sjdbGTFfile`, in which case no `counts/gene_counts.tsv` is written.

**3. Check that a compute node reaches ENA and NCBI**, inside an allocation:

```bash
salloc -M inter -p cm4_inter -N 1 -t 00:10:00
curl -s -o /dev/null -w 'ENA  %{http_code}\n' \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR12285169&result=read_run&fields=fastq_ftp"
curl -s -o /dev/null -w 'NCBI %{http_code}\n' \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=riboseq&retmax=1"
timeout 20 aria2c -x4 -s4 --max-download-limit=5M -d /tmp -o smoke.gz \
  https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR122/069/SRR12285169/SRR12285169.fastq.gz; ls -la /tmp/smoke.gz*
```

Both `curl`s must print 200 and `aria2c` must have fetched data. If not, export
`http_proxy` / `https_proxy` at the top of `slurm/node_worker.sh` and `slurm/prep.sh`.

**4. Run the pilot.** `slurm/master.sh config/config_pilot.json` runs the same chain
on the 18 runs in `config/pilot_runs.txt` (1 job x 2 nodes, 4 h); set its paths as
well. Expect COMPLETE in `logs/status.txt` and BAMs in `merged/bams/`.

**5. Screen with QC first.** The default query is broad and returns thousands of
candidates. With `pipeline.end: "qc"` each run is judged from a streamed sample of
200k reads ([QC.md](QC.md)) and nothing is downloaded. Then set `"bam"` and submit
`run`; the verdicts are reused. `query.max_runs` caps the cohort.

## The `_slurm` block

| key | default | meaning |
|---|---|---|
| `jobs` | 2 | array size. `cm4_std` runs 2 jobs at once; further tasks queue |
| `nodes_per_job` | 4 | nodes per array task; `cm4_std` takes 2-4 |
| `cpus_per_node`, `mem_per_node` | 112, `480G` | cores and memory per node for `run` |
| `time` | `24:00:00` | wall time of `run` (the `cm4_std` maximum) |
| `prep_cpus`, `prep_time` / `merge_cpus`, `merge_time` | 28, `12:00:00` / 28, `04:00:00` | the `prep` and `merge` jobs |
| `conda_env` | `ribomine` | environment activated in every job |
| `job_prefix` | `rm-` | job names: `<prefix>prep`, `<prefix>run`, `<prefix>merge` |
| `mail_user` | none | address for begin/end/fail mails |

The block is read only by `slurm/master.sh`; RiboMine ignores keys starting with
`_`. There are `jobs x nodes_per_job` shards, one per node. Each node runs
`project.jobs` samples at once with `project.threads` threads each. The product may
exceed the physical cores, because a sample waiting on a download uses no core;
`master.sh` warns above twice `cpus_per_node`. `download.connections` is per run, so
ENA sees nodes x `project.jobs` x `download.connections` connections; lower it first
if downloads fail.

## Rounds and `logs/status.txt`

`run` is limited to 24 h, so a large cohort needs several rounds. `merge` decides
from the per-sample JSONs what is finished and writes `<workdir>/logs/status.txt`.
Its first line is `COMPLETE` or `RESUBMIT NEEDED`; in the second case run
`slurm/master.sh <config> run` from a login node (compute nodes cannot call `sbatch`).

- Finished samples are skipped; an interrupted sample has no JSON and is redone.
- The split is pinned: an accession keeps its shard across rounds and re-runs of
  `prep`. Do not reduce the number of shards once results exist.
- `still to do` is what the next round will process. `dropped by QC` and `dropped,
  no arch` are finished. `tried and failed` lists runs with an entry in
  `merged/failed.tsv`; they do not cause another round on their own.
- Results are in `<workdir>/merged/` (`qc/qc_summary.tsv`, `mapping_summary.tsv`,
  `counts/gene_counts.tsv`, `bams/`, ...), SLURM logs in `<workdir>/logs/slurm/`.

## Memory and disk

Each node loads the STAR index into shared memory once (about 30 GB) and all of its
samples attach to it; add roughly 8 GB per sample in flight. If a node runs out of
memory, lower `project.jobs` (`config_lrz.json` uses 14). `node_worker.sh` loads the
genome itself and exits with an error if that fails, because RiboMine's fallback
(one index copy per sample) would exceed the node's memory. At startup it removes
the user's stale shared-memory segments over 1 GiB; a load that hangs for 900 s is
purged and retried once. By hand, on the node: `ipcs -m`, then `ipcrm -m <shmid>`.

Per run, measured on the pilot: 2.08 GB downloaded, 298 MB of BAM kept, 4.85 GB on
disk while it is processed. FASTQs are deleted once the BAM exists. BAMs accumulate
to about 1.0 TB for 3,300 runs (unless `keep.bam: false`). Transient space is
4.85 GB x runs in flight (nodes x `project.jobs`): about 0.55 TB for 8 x 14, 1.1 TB
for 8 x 28. A QC-only pass needs almost no disk.

`project.workdir` must be on a project filesystem (a DSS container), not in the home
quota; check with `dssusrinfo all` and `df -h <your-dss-path>`. Leave the download
tmpdir on the workdir's filesystem (the default): a download ends with a rename.
