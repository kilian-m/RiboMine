#!/bin/bash
# prep -- everything that must happen exactly once, before eight nodes start at once.
#
#   1. validate the reference paths (a typo here otherwise kills a 24 h job)
#   2. `ribomine setup`: the GTF index and the bowtie2 contaminant index, into
#      <workdir>/refs/, which every node then symlinks and only READS
#   3. the genome .fai -- pysam builds it silently and in place otherwise, so
#      eight workers would build it into the same file at the same time and leave
#      a truncated one, after which every architecture call is quietly fabricated
#   4. the split of meta/candidates.tsv into one accession list + one workdir per node
#
# It does NOT run the query. slurm/master.sh does that on the login node before
# submitting this, because the ENA portal search cannot complete from a cm4 compute
# node -- ENA takes minutes to answer it and the outbound path drops the idle
# connection. See the long note in master.sh. This job asserts the search already
# happened rather than quietly trying it again from the one place it cannot work.
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

# ------------------------------------------------------------------ #
# 1. the reference. Fail here, in minutes, not there, in a day.       #
# ------------------------------------------------------------------ #
fail=0
for p in "${GENOME}" "${GTF}"; do
    if [[ ! -f "${p}" ]]; then echo "[prep] MISSING FILE: ${p}" >&2; fail=1; fi
done
if [[ ! -d "${STAR_INDEX}" ]]; then
    echo "[prep] MISSING STAR INDEX: ${STAR_INDEX}" >&2
    fail=1
fi
(( fail )) && { echo "[prep] fix reference.* in ${CONFIG_FILE}" >&2; exit 2; }

# The gene-count matrix comes from STAR's own --quantMode GeneCounts, which only
# works if the annotation is IN the index. It cannot be supplied at mapping time:
# STAR refuses to insert junctions on the fly against a shared-memory genome, and
# the shared genome is what makes 28 concurrent samples per node affordable. So an
# index built without --sjdbGTFfile means no matrix -- say so now, not in 20 hours.
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

# ------------------------------------------------------------------ #
# 2 + 3. the shared indexes: built once, here, read by every node.    #
# ------------------------------------------------------------------ #
echo "[prep] building the annotation + contaminant indexes"
ribomine setup -c "${CONFIG_FILE}"

if [[ ! -s "${GENOME}.fai" ]]; then
    echo "[prep] indexing the genome FASTA (once)"
    samtools faidx "${GENOME}"
else
    echo "[prep] genome .fai present"
fi

# ------------------------------------------------------------------ #
# 4. the split. The query already happened, on the login node.        #
# ------------------------------------------------------------------ #
# Assert it, rather than let prepare_run.py fall through to searching ENA from
# here: that is the one place the search cannot work, it takes ~10 minutes of
# timeouts to discover, and it takes the job down when it does.
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
