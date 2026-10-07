#!/usr/bin/env bash
#SBATCH --job-name=m6p1a_fastq
#SBATCH --partition=qcpu_18i
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=0-3
#SBATCH --output=benchmark/real_world/results/m6_p1a/slurm_%A_%a.log
#SBATCH --error=benchmark/real_world/results/m6_p1a/slurm_%A_%a.err

# Four independent smoke runs. Not a scientific ProkaDiff job.
# 4 CPU / 16G / 2h: largest object is about 155 MB and 1.6 million spots.
# Scratch stays on the project filesystem, which has multi-petabyte free space.
set -euo pipefail
RUNS=(SRR28349244 SRR28357593 SRR28349243 SRR28357592)
RUN="${RUNS[$SLURM_ARRAY_TASK_ID]}"
REPO_DIR="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff"
cd "$REPO_DIR"
bash "$REPO_DIR/benchmark/real_world/scripts/stage_m6_p1a_fastq.sh" "$RUN"
