#!/bin/bash
# run -- the big one. 896 cores: an ARRAY OF 2 JOBS, each 4 nodes x 112 cores.
#
# Why an array and not one 8-node job: cm4_std caps a single job at 4 nodes (448
# cores) but will run 2 of them concurrently. So the only way to hold 896 cores is
# two 4-node jobs, and `--array=0-1 --nodes=4` is that. (--nodes and --array are
# both set by slurm/master.sh from the config, and override the headers below.)
#
# srun starts exactly ONE `ribomine run` per node. Not two, and not one per core:
# each node loads the ~28 GB STAR index into its own shared memory and every one
# of that node's 28 concurrent samples attaches to that single copy. A second
# ribomine on the same node would race it -- whichever finished first would call
# `STAR --genomeLoad Remove` and pull the index out from under the other.
#
# A node's shard is its global index across the whole array:
#     shard = SLURM_ARRAY_TASK_ID * NODES_PER_JOB + SLURM_PROCID
# (SLURM_PROCID is the node ordinal, because there is one task per node.) That
# must agree with the N_SHARDS that prepare_run.py split the cohort into, or a
# node reads someone else's list.
#
# Required environment (set by slurm/master.sh via --export):
#   REPO_DIR, CONFIG_FILE, WORK_DIR, CONDA_ENV,
#   N_SHARDS, NODES_PER_JOB, CPUS_PER_NODE, RM_JOBS, RM_THREADS
#
#SBATCH --job-name=rm-run
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_std
#SBATCH --qos=cm4_std
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --mem=480G
#SBATCH --time=24:00:00
# SIGTERM 10 min before the wall-time SIGKILL. RiboMine has no signal handler --
# it does not need one: every stage writes its JSON atomically and only *after*
# the work it records is complete, so a sample killed mid-flight simply has no
# JSON and is redone next round. The 10 minutes are for node_worker.sh's own trap,
# which releases the shared-memory STAR genome -- a segment leaked by a SIGKILL
# holds 28 GB on that node until someone removes it by hand, and the next job to
# land there cannot load its own.
#SBATCH --signal=B:TERM@600

set -euo pipefail

echo "[run] array-task=${SLURM_ARRAY_TASK_ID:-0} job=${SLURM_JOB_ID}" \
     "nodes=${SLURM_JOB_NUM_NODES} start $(date -Is)"

srun --ntasks="${SLURM_JOB_NUM_NODES}" \
     --ntasks-per-node=1 \
     --cpus-per-task="${CPUS_PER_NODE}" \
     bash "${REPO_DIR}/slurm/node_worker.sh"

echo "[run] array-task=${SLURM_ARRAY_TASK_ID:-0} done $(date -Is)"
