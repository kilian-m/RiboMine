#!/bin/bash
# merge: write the cohort tables and logs/status.txt.
#
# Submitted with --dependency=afterany on the run array, so it runs however the
# array ended. What finished is read from the per-sample JSONs, not from exit
# codes. status.txt says COMPLETE, or RESUBMIT NEEDED with the command to run from
# a login node (compute nodes cannot call sbatch).
#
# Submitted by slurm/master.sh. Required environment:
#   REPO_DIR, CONFIG_FILE, WORK_DIR, CONDA_ENV, N_SHARDS
#
#SBATCH --job-name=rm-merge
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=100G
#SBATCH --time=04:00:00

set -euo pipefail

echo "[merge] $(hostname) start $(date -Is)"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"
cd "${REPO_DIR}"

export MPLBACKEND=Agg

python scripts/merge_results.py -c "${CONFIG_FILE}" --shards "${N_SHARDS}"

echo "[merge] done $(date -Is)"
echo
echo "===================== ${WORK_DIR}/logs/status.txt ====================="
cat "${WORK_DIR}/logs/status.txt"
