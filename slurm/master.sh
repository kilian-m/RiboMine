#!/bin/bash
# RiboMine on LRZ CoolMUC-4 -- the submitter. Run it on a LOGIN node.
#
#   slurm/master.sh config/config_lrz.json            # prep -> run -> merge
#   slurm/master.sh config/config_lrz.json run        # another round (the usual case)
#   slurm/master.sh config/config_lrz.json prep       # just the query + the split
#   slurm/master.sh config/config_lrz.json merge      # just re-merge + re-report
#   DRY_RUN=1 slurm/master.sh config/config_lrz.json  # print the sbatch lines, submit nothing
#
# The chain, and why it is shaped like this:
#
#   prep   cm4_tiny, 1 node   query ENA/NCBI -> candidates.tsv; build the GTF and
#                             contaminant indexes; split the cohort into one
#                             accession list + one workdir per node.
#                             Everything that must happen exactly once.
#
#   run    cm4_std, ARRAY OF 2 x 4 NODES = 896 cores
#                             cm4_std takes at most 4 nodes (448 cores) in ONE job,
#                             but will run 2 jobs at once -- so 896 cores is an array
#                             of two 4-node jobs, not one 8-node job. Each of the 8
#                             nodes runs one `ribomine run` over its own shard.
#
#   merge  cm4_tiny, afterany run
#                             Symlinks the 8 workdirs into one, writes the cohort
#                             tables, and says in logs/status.txt whether anything
#                             is left. `afterany`, so it runs however the array ended
#                             -- clean, out of wall time, or with a node dead.
#
# The run has a 24 h ceiling and a cohort of hundreds of runs does not fit in it.
# So expect to resubmit: `slurm/master.sh <config> run`, once per round, until
# status.txt says COMPLETE. Each round skips what is already finished. It is manual
# because a cm4 compute node is not permitted to sbatch, so no job in the chain can
# launch the next one.

set -euo pipefail

CONFIG_FILE="${1:-}"
STEP="${2:-all}"

if [[ -z "${CONFIG_FILE}" ]]; then
    echo "usage: $0 <config.json> [all|prep|run|merge]" >&2
    exit 2
fi
if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "config not found: ${CONFIG_FILE}" >&2
    exit 2
fi
case "${STEP}" in
    all|prep|run|merge) ;;
    *) echo "unknown step '${STEP}' (want: all, prep, run, merge)" >&2; exit 2 ;;
esac

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
CONFIG_ABS="$(readlink -f "${CONFIG_FILE}")"

# The config is read with python, not jq: jq is not on every LRZ login node, and
# python3 is. (RiboMine ignores every `_`-prefixed key, which is why the cluster
# shape can live in the same file as the pipeline settings.)
cfg() {
    python3 - "${CONFIG_ABS}" "$1" "${2-}" <<'PY'
import json, sys
cur = json.load(open(sys.argv[1]))
for part in sys.argv[2].split('.'):
    cur = cur.get(part) if isinstance(cur, dict) else None
    if cur is None:
        break
print(sys.argv[3] if cur is None else cur)
PY
}

WORK_DIR=$(cfg project.workdir)
CONDA_ENV=$(cfg _slurm.conda_env ribomine)
JOB_PREFIX=$(cfg _slurm.job_prefix "rm-")
MAIL_USER=$(cfg _slurm.mail_user "")

S_JOBS=$(cfg _slurm.jobs 2)
S_NODES=$(cfg _slurm.nodes_per_job 4)
S_CPUS=$(cfg _slurm.cpus_per_node 112)
S_MEM=$(cfg _slurm.mem_per_node 480G)
S_TIME=$(cfg _slurm.time 24:00:00)

PREP_CPUS=$(cfg _slurm.prep_cpus 28)
PREP_TIME=$(cfg _slurm.prep_time 12:00:00)
MERGE_CPUS=$(cfg _slurm.merge_cpus 28)
MERGE_TIME=$(cfg _slurm.merge_time 04:00:00)

RM_JOBS=$(cfg project.jobs 14)
RM_THREADS=$(cfg project.threads 8)

N_SHARDS=$(( S_JOBS * S_NODES ))
TOTAL_CORES=$(( N_SHARDS * S_CPUS ))

if [[ -z "${WORK_DIR}" || "${WORK_DIR}" == "None" ]]; then
    echo "project.workdir is not set in ${CONFIG_FILE}" >&2
    exit 2
fi
if (( S_NODES < 2 || S_NODES > 4 )); then
    echo "cm4_std takes 2-4 nodes per job; _slurm.nodes_per_job is ${S_NODES}" >&2
    exit 2
fi
if (( S_JOBS > 2 )); then
    echo "note: cm4_std runs only 2 jobs at a time -- array tasks beyond the 2nd" >&2
    echo "      will queue rather than add cores." >&2
fi
# A cm4 node has ${S_CPUS} physical cores and twice that many hardware threads, and
# running up to 2x the physical cores' worth of samples is deliberate: a sample
# blocked on an ENA download is not using a core. Past 2x there are not even
# hardware threads left to hold them, and the node only thrashes.
if (( RM_JOBS * RM_THREADS > 2 * S_CPUS )); then
    echo "WARNING: project.jobs x project.threads = $(( RM_JOBS * RM_THREADS )) is more than" >&2
    echo "         2 x ${S_CPUS} = $(( 2 * S_CPUS )) hardware threads on a node. It will thrash." >&2
fi

LOG_DIR="${WORK_DIR}/logs/slurm"
mkdir -p "${LOG_DIR}"

echo "==> repo        : ${REPO_DIR}"
echo "==> config      : ${CONFIG_ABS}"
echo "==> workdir     : ${WORK_DIR}"
echo "==> conda env   : ${CONDA_ENV}"
echo "==> step        : ${STEP}"
echo "==> shape       : ${S_JOBS} job(s) x ${S_NODES} node(s) x ${S_CPUS} cores = ${TOTAL_CORES} cores"
echo "==>               ${N_SHARDS} shard(s); each node runs ${RM_JOBS} dataset(s) x ${RM_THREADS} thread(s)"
echo "==> end point   : $(cfg pipeline.end bam)"

MAIL_ARGS=()
if [[ -n "${MAIL_USER}" && "${MAIL_USER}" != "None" ]]; then
    MAIL_ARGS=(--mail-type=begin,end,fail --mail-user="${MAIL_USER}")
fi

submit() {
    local label="$1"; shift
    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        echo "[DRY] sbatch $*" >&2
        echo "1"
        return 0
    fi
    echo "[submit] ${label}" >&2
    # LRZ is a multi-cluster SLURM: `--parsable` returns "<jobid>;<cluster>", and
    # --dependency does not parse the suffix.
    sbatch --parsable "$@" | cut -d';' -f1
}

COMMON_EXPORT="ALL,REPO_DIR=${REPO_DIR},CONFIG_FILE=${CONFIG_ABS},WORK_DIR=${WORK_DIR}"
COMMON_EXPORT+=",CONDA_ENV=${CONDA_ENV},N_SHARDS=${N_SHARDS}"

PREP_JID=""
RUN_JID=""
MERGE_JID=""

# ---------------------------------------------------------------- prep ----- #
if [[ "${STEP}" == "all" || "${STEP}" == "prep" ]]; then
    PREP_JID=$(submit "prep" \
        --job-name="${JOB_PREFIX}prep" \
        --cpus-per-task="${PREP_CPUS}" \
        --time="${PREP_TIME}" \
        --output="${LOG_DIR}/prep_%j.out" \
        --error="${LOG_DIR}/prep_%j.err" \
        "${MAIL_ARGS[@]}" \
        --export="${COMMON_EXPORT}" \
        "${REPO_DIR}/slurm/prep.sh")
    echo "prep  job id: ${PREP_JID}"
fi

# ---------------------------------------------------------------- run ------ #
if [[ "${STEP}" == "all" || "${STEP}" == "run" ]]; then
    dep=()
    # afterok, not afterany: without the shard lists and the shared indexes that
    # prep writes, every node would fail at once -- and eight nodes each building
    # the same bowtie2 index into the same path is corruption, not a slow start.
    [[ -n "${PREP_JID}" ]] && dep=(--dependency="afterok:${PREP_JID}")

    RUN_JID=$(submit "run" \
        --job-name="${JOB_PREFIX}run" \
        --array="0-$(( S_JOBS - 1 ))" \
        --nodes="${S_NODES}" \
        --cpus-per-task="${S_CPUS}" \
        --mem="${S_MEM}" \
        --time="${S_TIME}" \
        "${dep[@]}" \
        --output="${LOG_DIR}/run_%A_%a.out" \
        --error="${LOG_DIR}/run_%A_%a.err" \
        "${MAIL_ARGS[@]}" \
        --export="${COMMON_EXPORT},NODES_PER_JOB=${S_NODES},CPUS_PER_NODE=${S_CPUS},RM_JOBS=${RM_JOBS},RM_THREADS=${RM_THREADS}" \
        "${REPO_DIR}/slurm/run.sh")
    echo "run   job id: ${RUN_JID}"
fi

# ---------------------------------------------------------------- merge ---- #
if [[ "${STEP}" == "all" || "${STEP}" == "run" || "${STEP}" == "merge" ]]; then
    dep=()
    # afterany: the array is EXPECTED to hit the 24 h wall and be killed. Merge
    # works out what finished from the per-sample JSONs, not from exit codes, so
    # a killed node and a clean one are the same to it.
    [[ -n "${RUN_JID}" ]] && dep=(--dependency="afterany:${RUN_JID}")

    MERGE_JID=$(submit "merge" \
        --job-name="${JOB_PREFIX}merge" \
        --cpus-per-task="${MERGE_CPUS}" \
        --time="${MERGE_TIME}" \
        "${dep[@]}" \
        --output="${LOG_DIR}/merge_%j.out" \
        --error="${LOG_DIR}/merge_%j.err" \
        "${MAIL_ARGS[@]}" \
        --export="${COMMON_EXPORT}" \
        "${REPO_DIR}/slurm/merge.sh")
    echo "merge job id: ${MERGE_JID}"
fi

cat <<EOF

Submitted:
  prep  : ${PREP_JID:-<skipped>}
  run   : ${RUN_JID:-<skipped>}
  merge : ${MERGE_JID:-<skipped>}

Watch it:
  squeue -M cm4 -u \$USER
  tail -f ${WORK_DIR}/logs/slurm/run_*_0.out

When merge has run, read ${WORK_DIR}/logs/status.txt:
  COMPLETE         -> done; the cohort tables are in ${WORK_DIR}/merged/
  RESUBMIT NEEDED  -> from a login node:  $0 ${CONFIG_FILE} run
                      (repeat until COMPLETE; each round only redoes what is unfinished)
EOF
