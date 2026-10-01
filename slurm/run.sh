#!/bin/bash
# run: the array job. Each array task is one multi-node job, and srun starts one
# `ribomine run` per node (slurm/node_worker.sh).
#
# --array and --nodes are set by slurm/master.sh from the config and override the
# headers below. cm4_std allows at most 4 nodes per job and 2 running jobs, so
# 8 nodes (896 cores) are an array of two 4-node jobs.
#
# One process per node: the node loads the STAR index into shared memory once and
# all of its concurrent samples attach to that copy. A second `ribomine run` on
# the same node could remove the genome while the first still uses it.
#
# A node's shard is its index across the whole array:
#     shard = SLURM_ARRAY_TASK_ID * NODES_PER_JOB + SLURM_PROCID
# (one task per node, so SLURM_PROCID is the node ordinal). This must match the
# N_SHARDS that prepare_run.py split the cohort into.
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
# SIGTERM 10 min before the wall-time limit, so that node_worker.sh can release the
# shared-memory STAR genome. RiboMine needs no signal handler: each stage writes
# its JSON atomically after the work is done, so an interrupted sample has no JSON
# and is redone in the next round.
#SBATCH --signal=B:TERM@600

set -euo pipefail

echo "[run] array-task=${SLURM_ARRAY_TASK_ID:-0} job=${SLURM_JOB_ID}" \
     "nodes=${SLURM_JOB_NUM_NODES} start $(date -Is)"

srun --ntasks="${SLURM_JOB_NUM_NODES}" \
     --ntasks-per-node=1 \
     --cpus-per-task="${CPUS_PER_NODE}" \
     bash "${REPO_DIR}/slurm/node_worker.sh"

echo "[run] array-task=${SLURM_ARRAY_TASK_ID:-0} done $(date -Is)"
