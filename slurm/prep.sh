#!/bin/bash
# prep: the steps that must run exactly once, before the nodes start.
#
#   1. check the reference paths
#   2. `ribomine setup`: the GTF index and the bowtie2 contaminant index, into
#      <workdir>/refs/, which every node symlinks and only reads
#   3. the genome .fai (pysam would otherwise build it in place, on every node
#      at the same time)
#   4. split the cohort into one accession list and one workdir per node
#
# The archive query is not run here: slurm/master.sh runs it on the login node
# (see the note there). This job only checks that meta/candidates.tsv exists.
#
# Submitted by slurm/master.sh. Required environment:
#   REPO_DIR, CONFIG_FILE, WORK_DIR, CONDA_ENV, N_SHARDS
#
#SBATCH --job-name=rm-prep
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=100G
#SBATCH --time=12:00:00

set -euo pipefail

echo "[prep] $(hostname) start $(date -Is)"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"
cd "${REPO_DIR}"

cfg() {
    python3 - "${CONFIG_FILE}" "$1" <<'PY'
import json, sys
cur = json.load(open(sys.argv[1]))
for part in sys.argv[2].split('.'):
    cur = cur.get(part) if isinstance(cur, dict) else None
    if cur is None:
        break
print("" if cur is None else cur)
PY
}

GENOME=$(cfg reference.genome_fasta)
GTF=$(cfg reference.gtf)
STAR_INDEX=$(cfg reference.star_index)

# --- 1. reference paths ---
fail=0
for p in "${GENOME}" "${GTF}"; do
    if [[ ! -f "${p}" ]]; then echo "[prep] MISSING FILE: ${p}" >&2; fail=1; fi
done
if [[ ! -d "${STAR_INDEX}" ]]; then
    echo "[prep] MISSING STAR INDEX: ${STAR_INDEX}" >&2
    fail=1
fi
(( fail )) && { echo "[prep] fix reference.* in ${CONFIG_FILE}" >&2; exit 2; }

# The gene-count matrix comes from STAR's --quantMode GeneCounts, which needs the
# annotation in the index (--sjdbGTFfile at build time): STAR cannot insert
# junctions at mapping time against a shared-memory genome. Without it no matrix
# is written, so warn now.
GENE_INFO="${STAR_INDEX}/geneInfo.tab"
N_GENES=0
[[ -f "${GENE_INFO}" ]] && N_GENES=$(head -1 "${GENE_INFO}" | tr -dc '0-9')
if [[ "${N_GENES:-0}" -gt 0 ]]; then
    echo "[prep] STAR index carries ${N_GENES} genes -> the gene-count matrix will be written"
else
    echo "[prep] WARNING: the STAR index has no annotation in it (geneInfo.tab is empty" >&2
    echo "[prep]          or absent). STAR cannot count reads into genes, so there will" >&2
    echo "[prep]          be NO counts/gene_counts.tsv. Everything else is unaffected." >&2
    echo "[prep]          To get one, rebuild the index with --sjdbGTFfile ${GTF}." >&2
fi

echo "[prep] free space on the workdir filesystem:"
df -h "$(dirname "${WORK_DIR}")" || true

# --- 2 + 3. shared indexes: built once here, read by every node ---
echo "[prep] building the annotation + contaminant indexes"
ribomine setup -c "${CONFIG_FILE}"

if [[ ! -s "${GENOME}.fai" ]]; then
    echo "[prep] indexing the genome FASTA (once)"
    samtools faidx "${GENOME}"
else
    echo "[prep] genome .fai present"
fi

# --- 4. the split ---
# With pipeline.start=query, require candidates.tsv: prepare_run.py would otherwise
# run the query from this compute node, where it cannot complete.
START=$(cfg pipeline.start)
if [[ "${START}" == "query" && ! -s "${WORK_DIR}/meta/candidates.tsv" ]]; then
    cat >&2 <<EOF
[prep] FATAL: no ${WORK_DIR}/meta/candidates.tsv

slurm/master.sh writes it, on the login node, before submitting this job -- the ENA
portal search cannot complete from a compute node (it holds one connection idle for
minutes and the outbound path drops it). This job will not retry it from here.

  * submit through slurm/master.sh, which does the query for you; or
  * run it yourself on a LOGIN node:  ribomine query -c ${CONFIG_FILE}
  * or drop an existing candidates.tsv into ${WORK_DIR}/meta/
EOF
    exit 2
fi

echo "[prep] splitting into ${N_SHARDS} shard(s)"
python scripts/prepare_run.py -c "${CONFIG_FILE}" --shards "${N_SHARDS}"

echo "[prep] done $(date -Is)"
