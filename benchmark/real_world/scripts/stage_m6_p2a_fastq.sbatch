#!/usr/bin/env bash
#SBATCH --job-name=m6p2a_fastq
#SBATCH --partition=qcpu_18i
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --array=0-22%6
#SBATCH --output=benchmark/real_world/results/prjna1088182/m6_p2a/staging/slurm_%A_%a.log
#SBATCH --error=benchmark/real_world/results/prjna1088182/m6_p2a/staging/slurm_%A_%a.err

# Stage one eligible SRA1824580 run that is not already a verified P1a FASTQ.
# 4 CPU / 16G / 6h: largest object is about 484 MB and 5.0 million spots.
# Does not run ProkaDiff, Bowtie2, or breseq.
set -euo pipefail
REPO_DIR="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff"
cd "$REPO_DIR"
mapfile -t RUNS < <(python3 "$REPO_DIR/benchmark/real_world/scripts/build_m6_p2a_cohort.py" --list-new-runs)
if [[ "${#RUNS[@]}" -ne 23 ]]; then
    echo "expected 23 new paired runs, found ${#RUNS[@]}" >&2
    exit 2
fi
RUN="${RUNS[$SLURM_ARRAY_TASK_ID]}"
echo "stage $RUN task ${SLURM_ARRAY_TASK_ID} job ${SLURM_JOB_ID:-manual}"
bash "$REPO_DIR/benchmark/real_world/scripts/stage_m6_p1a_fastq.sh" "$RUN" p2a
