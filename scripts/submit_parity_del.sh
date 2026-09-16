#!/usr/bin/env bash
#SBATCH -p qcpu_18i
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -J pd_del_test
#SBATCH -o benchmark/logs/%x-%j.out
#SBATCH -e benchmark/logs/%x-%j.err

set -euo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"
CARGO_TARGET_DIR="${ROOT}/target/slurm-${SLURM_JOB_ID}"
export CARGO_TARGET_DIR

CONDA_ENV="${PROKDIFF_CONDA_ENV:-/hpcfs/fhome/caizhh/.conda/envs/prokadiff}"
export PATH="${CONDA_ENV}/bin:${PATH:-}"

FIXTURE_DIR="${ROOT}/testdata/parity/synthetic_del"
mkdir -p "${FIXTURE_DIR}" "${ROOT}/benchmark/logs" "${ROOT}/benchmark/results"

python3 "${FIXTURE_DIR}/generate_del.py"

OUTDIR="${ROOT}/benchmark/results/synthetic_del_${SLURM_JOB_ID}"
mkdir -p "${OUTDIR}"

echo "=== Building release prokadiff ==="
cargo build --release -p prokadiff
export PROKADIFF_BIN="${CARGO_TARGET_DIR}/release/prokadiff"

echo "=== Running benchmark_case on synthetic DEL ==="
bash "${ROOT}/scripts/parity/benchmark_case.sh" \
  --dataset "structural_del" \
  --ref "${FIXTURE_DIR}/reference.fa" \
  --r1 "${FIXTURE_DIR}/reads.fq" \
  --outdir "${OUTDIR}" \
  --threads 4

echo "=== Done synthetic DEL benchmark ==="
