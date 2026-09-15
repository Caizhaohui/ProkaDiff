#!/usr/bin/env bash
#SBATCH -p qcpu_18i
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -J pd_sub_test
#SBATCH -o benchmark/logs/%x-%j.out
#SBATCH -e benchmark/logs/%x-%j.err

set -euo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"

CONDA_ENV="${PROKDIFF_CONDA_ENV:-/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/.conda/envs/BactGenome}"
if [[ ! -d "${CONDA_ENV}" ]]; then
  CONDA_ENV="/hpcfs/fhome/caizhh/.conda/envs/BactGenome"
fi
export PATH="${CONDA_ENV}/bin:${PATH:-}"

FIXTURE_DIR="${ROOT}/testdata/parity/synthetic_sub"
OUTDIR="${ROOT}/benchmark/results/synthetic_sub_${SLURM_JOB_ID}"
mkdir -p "${OUTDIR}"

echo "=== Building release prokadiff ==="
cargo build --release -p prokadiff

echo "=== Running benchmark_case on synthetic SUB ==="
bash "${ROOT}/scripts/parity/benchmark_case.sh" \
  --dataset "sub_adjacent" \
  --ref "${FIXTURE_DIR}/reference.fa" \
  --r1 "${FIXTURE_DIR}/reads.fq" \
  --outdir "${OUTDIR}" \
  --breseq-gd "${FIXTURE_DIR}/oracle.gd" \
  --threads 4

echo "=== Done synthetic SUB benchmark ==="
