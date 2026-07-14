#!/bin/bash
# One node of the run: one `ribomine run` over one shard of the cohort.
#
# srun'd once per node by slurm/run.sh. Everything here is about the two things
# that are a node's own private business -- its shard, and its shared-memory STAR
# genome -- so that the eight nodes never touch the same object.
#
# THE STAR SHARED-MEMORY GENOME. Each node loads the ~28 GB index into shared
# memory once and lets its 28 concurrent samples attach to that one copy; without
# it the node would need 28 x 28 GB and die. RiboMine drops the segment in a
# `finally` -- but a `finally` does not run when SLURM SIGKILLs the process at the
# wall-time limit, and a leaked segment holds 28 GB on that node and makes the
# next job's `LoadAndExit` fail. So the segment is removed here on the way OUT
# (the trap) and, more importantly, on the way IN: a stale segment from whatever
# ran here before is cleared before we load ours. That is the recovery path that
# does not depend on the previous job having exited politely. It is safe because
# cm4_std allocates nodes exclusively -- there is no other job here whose genome
# we could be removing.
#
# Required environment (inherited from slurm/run.sh through srun):
#   REPO_DIR, CONFIG_FILE, WORK_DIR, CONDA_ENV,
#   N_SHARDS, NODES_PER_JOB, RM_JOBS, RM_THREADS

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

# 28 samples run at once on this node, and numpy/matplotlib in each of them will
# otherwise each open a BLAS thread pool the size of the node. Cap them: the
# parallelism that matters here is `--jobs`, and STAR/bowtie2 get their threads
# from `--threads`.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg

# Node-local NVMe for tool scratch. NOT for RiboMine's own download tmpdir --
# that one must stay on the workdir's filesystem, because the ENA route finishes
# with os.replace(tmp, out.fastq.gz), which is a rename and fails across devices.
# It loses nothing: a rename within the workdir is free, whereas /tmp -> DSS would
# have been a multi-GB copy of every run.
export TMPDIR="/tmp/${USER}/ribomine.${SLURM_JOB_ID}.${SHARD_ID}"
mkdir -p "${TMPDIR}"

STAR_INDEX=$(python3 - "${CONFIG_FILE}" <<'PY'
import json, sys
print(json.load(open(sys.argv[1])).get("reference", {}).get("star_index", ""))
PY
)

# A STAR killed part-way through loading the genome leaves its shared segment flagged
# "load in progress" with nobody loading it, and every later STAR then waits on that
# flag FOREVER -- "Another job is still loading the genome, sleeping for 1 min", over
# and over. `--genomeLoad Remove` does not reliably clear that state, so the segment
# has to be taken out with ipcrm. Which is worth doing carefully, because a hang is
# strictly worse than a crash here: a crash loses one shard and says so, while a hang
# burns the whole 24 h wall clock in silence and produces nothing.
#
# The >1 GiB filter is what keeps this from being a blunt instrument: it removes the
# genome segments and leaves anything small alone. Only this user's segments are even
# visible to ipcrm, and cm4_std allocates nodes exclusively, so there is no other job
# of ours here to take down with it.
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

# `timeout`, not patience: a stuck segment makes LoadAndExit wait forever rather than
# fail, so without a clock the node hangs instead of erroring. 15 min is many times
# what a real load takes (~2-5 min from a warm page cache) and a fraction of the wall.
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

# Clear whatever the last job on this node left behind, before we load ours.
star_genome_remove

# ------------------------------------------------------------------ #
# Load the shared genome HERE, and refuse to run without it.          #
# ------------------------------------------------------------------ #
# RiboMine loads it itself, and if the load fails it warns and falls back to
# NoSharedMemory -- every worker then loads its OWN ~30 GB copy of the index. On a
# workstation that is a memory bill. Here it is suicide: RM_JOBS is 28, so the
# fallback asks the node for 28 x 30 GB = 840 GB, and the node OOM-kills the whole
# process within seconds. The shard is then lost for the round, and the SLURM log
# says only "Killed" -- with the actual cause, one WARNING line, scrolled far above.
#
# So the decision is taken here instead, where it can be fatal on purpose: load the
# segment, and if it will not load, STOP. Exporting the result in
# RIBOMINE_STAR_GENOME_LOAD is what RiboMine reads back (star.genome_load()), so it
# skips its own load and its own fallback entirely.
if [[ -n "${STAR_INDEX}" ]]; then
    echo "[node ${HOST}] loading the STAR genome into shared memory (~30 GB, once for this node)"
    rc=0
    star_genome_load || rc=$?

    # 124 is `timeout`'s: the load did not fail, it HUNG -- which is the signature of a
    # segment left mid-load by a killed STAR. The purge above should have caught it, so
    # this is the belt to that braces (a segment created between the purge and now, or
    # one that Remove revived). Purge again and take exactly one more run at it: a
    # second hang is a real problem and must not be slept through.
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
# SLURM's pre-kill SIGTERM: pass it on and let the EXIT trap release the genome.
# The sample in flight is simply lost -- it has written no JSON, so next round
# redoes it. That is the whole recovery mechanism, and it is the same one for a
# SIGKILL and for a dead node.
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
