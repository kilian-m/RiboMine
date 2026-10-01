#!/bin/bash
# One node of the run: one `ribomine run` over one shard. Started once per node by
# slurm/run.sh.
#
# The node loads the STAR index (~28 GB) into shared memory once and its concurrent
# samples attach to that copy. A segment left behind by a killed job blocks the
# next load, so stale segments are removed at startup and again on exit (trap).
# This is safe because cm4_std allocates nodes exclusively.
#
# Required environment (inherited from slurm/run.sh through srun):
#   REPO_DIR, CONFIG_FILE, WORK_DIR, CONDA_ENV,
#   N_SHARDS, NODES_PER_JOB, RM_JOBS, RM_THREADS
# Optional: STAR_LOAD_TIMEOUT (seconds, default 900)

set -uo pipefail

ARRAY_TASK="${SLURM_ARRAY_TASK_ID:-0}"
NODE_RANK="${SLURM_PROCID:-0}"
SHARD=$(( ARRAY_TASK * NODES_PER_JOB + NODE_RANK ))
HOST="$(hostname -s)"

if (( SHARD >= N_SHARDS )); then
    echo "[node ${HOST}] shard ${SHARD} >= N_SHARDS ${N_SHARDS} -- nothing assigned"
    exit 0
fi

SHARD_ID=$(printf "%02d" "${SHARD}")
ACC_LIST="${WORK_DIR}/shards/shard_${SHARD_ID}.txt"
SHARD_WD="${WORK_DIR}/shards/work_${SHARD_ID}"

echo "[node ${HOST}] shard=${SHARD_ID} array-task=${ARRAY_TASK} rank=${NODE_RANK} start $(date -Is)"

if [[ ! -f "${ACC_LIST}" ]]; then
    echo "[node ${HOST}] no accession list: ${ACC_LIST} -- did prep run?" >&2
    exit 1
fi
N_ACC=$(grep -cve '^\s*$' -e '^#' "${ACC_LIST}" || true)
echo "[node ${HOST}] ${N_ACC} run(s) in ${ACC_LIST}"
if (( N_ACC == 0 )); then
    echo "[node ${HOST}] empty shard -- done"
    exit 0
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"
cd "${REPO_DIR}"

# One BLAS/OpenMP thread per process: the parallelism comes from `--jobs`, and
# STAR/bowtie2 take their threads from `--threads`.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg

# Node-local scratch for the tools. RiboMine's download tmpdir stays on the
# workdir's filesystem: the ENA download ends with a rename, which fails across
# devices.
export TMPDIR="/tmp/${USER}/ribomine.${SLURM_JOB_ID}.${SHARD_ID}"
mkdir -p "${TMPDIR}"

STAR_INDEX=$(python3 - "${CONFIG_FILE}" <<'PY'
import json, sys
print(json.load(open(sys.argv[1])).get("reference", {}).get("star_index", ""))
PY
)

# A STAR killed while loading the genome leaves its shared segment flagged "load in
# progress", and every later STAR waits on it indefinitely. `--genomeLoad Remove`
# does not reliably clear that state, so this user's segments larger than 1 GiB
# (the genome) are removed with ipcrm; smaller segments are left alone.
star_shm_purge() {
    local id
    for id in $(ipcs -m 2>/dev/null \
                | awk -v u="${USER}" '$3 == u && $5 + 0 > 1073741824 {print $2}'); do
        echo "[node ${HOST}] ipcrm stale shared-memory segment ${id} (a killed STAR left it behind)"
        ipcrm -m "${id}" 2>/dev/null || true
    done
}

star_genome_remove() {
    [[ -n "${STAR_INDEX}" ]] || return 0
    STAR --genomeLoad Remove --genomeDir "${STAR_INDEX}" \
         --outFileNamePrefix "${TMPDIR}/_star_remove." >/dev/null 2>&1 || true
    star_shm_purge
}

# Time limit for the load (seconds): on a stuck segment LoadAndExit waits instead
# of failing. A normal load takes 2-5 min.
STAR_LOAD_TIMEOUT="${STAR_LOAD_TIMEOUT:-900}"

star_genome_load() {
    timeout "${STAR_LOAD_TIMEOUT}" \
        STAR --genomeLoad LoadAndExit --genomeDir "${STAR_INDEX}" \
             --outFileNamePrefix "${TMPDIR}/_star_load." >/dev/null 2>&1
}

cleanup() {
    star_genome_remove
    rm -rf "${TMPDIR}"
    echo "[node ${HOST}] shard=${SHARD_ID} exit $(date -Is)"
}
trap cleanup EXIT

# Remove whatever the previous job on this node left behind.
star_genome_remove

# Load the shared genome here and stop if it fails. If RiboMine's own load fails,
# it falls back to NoSharedMemory: one ~30 GB copy of the index per concurrent
# sample, which exceeds the node's memory. With RIBOMINE_STAR_GENOME_LOAD set,
# RiboMine (star.load_genome()) skips its own load and fallback.
if [[ -n "${STAR_INDEX}" ]]; then
    echo "[node ${HOST}] loading the STAR genome into shared memory (~30 GB, once for this node)"
    rc=0
    star_genome_load || rc=$?

    # 124 = `timeout`: the load hung on a segment left mid-load. Purge again and
    # retry once.
    if (( rc == 124 )); then
        echo "[node ${HOST}] the genome load HUNG for ${STAR_LOAD_TIMEOUT}s -- a killed STAR" >&2
        echo "[node ${HOST}] left its segment flagged 'loading'. Purging it and retrying once." >&2
        star_genome_remove
        rc=0
        star_genome_load || rc=$?
    fi

    if (( rc == 0 )); then
        export RIBOMINE_STAR_GENOME_LOAD=LoadAndKeep
        echo "[node ${HOST}] STAR shared genome loaded"
    else
        cat >&2 <<EOF
[node ${HOST}] FATAL: the STAR genome would not load into shared memory (exit ${rc}).

  index : ${STAR_INDEX}
  needs : ~30 GB resident (the SA file alone is ~25 GB)

Not falling back to NoSharedMemory: with --jobs ${RM_JOBS} that would ask this node
for ${RM_JOBS} separate copies of the index, and it would be OOM-killed on the spot.
Refusing outright loses the same shard, and says why.

Usual causes, in order:
  * the job did not ask for enough memory        -> #SBATCH --mem in slurm/run.sh
  * (exit 124 = it hung, twice) a shared segment stuck mid-load. By hand, ON THIS NODE:
        ipcs -m                 # the ~30 GB one owned by ${USER}
        ipcrm -m <shmid>
  * SHMALL/SHMMAX too small for a 30 GB segment (ask LRZ)
EOF
        exit 1
    fi
fi

RM_PID=""
# On SLURM's SIGTERM, pass it on to ribomine; the EXIT trap releases the genome.
# A sample in flight has written no JSON, so the next round redoes it.
trap 'echo "[node '"${HOST}"'] SIGTERM"; [[ -n "${RM_PID}" ]] && kill -TERM "${RM_PID}" 2>/dev/null; ' TERM

ribomine run \
    -c "${CONFIG_FILE}" \
    --workdir "${SHARD_WD}" \
    --accessions "${ACC_LIST}" \
    --jobs "${RM_JOBS}" \
    --threads "${RM_THREADS}" &
RM_PID=$!

rc=0
wait "${RM_PID}" || rc=$?
echo "[node ${HOST}] ribomine exited ${rc}"
exit "${rc}"
