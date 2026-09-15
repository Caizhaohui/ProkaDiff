#!/usr/bin/env bash
#SBATCH -p qcpu_18i
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -J pd_ambig
#SBATCH -o benchmark/logs/%x-%j.out
#SBATCH -e benchmark/logs/%x-%j.err

set -euo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"

CONDA_ENV="${PROKDIFF_CONDA_ENV:-/hpcfs/fhome/caizhh/.conda/envs/BactGenome}"
export PATH="${CONDA_ENV}/bin:${PATH:-}"

mkdir -p benchmark/logs benchmark/results/ambiguous_base

FIXTURE_DIR="${ROOT}/testdata/parity/ambiguous_base"
python3 "${FIXTURE_DIR}/generate_ambiguous.py"

OUTDIR="${ROOT}/benchmark/results/ambiguous_base_${SLURM_JOB_ID}"
mkdir -p "${OUTDIR}"

echo "=== Building prokadiff release ==="
cargo build --release -p prokadiff

echo "=== Running benchmark_case for ambiguous_base ==="
bash "${ROOT}/scripts/parity/benchmark_case.sh" \
  --dataset "ambiguous_n" \
  --ref "${FIXTURE_DIR}/reference.fa" \
  --r1 "${FIXTURE_DIR}/reads.fq" \
  --outdir "${OUTDIR}" \
  --threads 4

cp "${OUTDIR}/breseq.gd" "${FIXTURE_DIR}/oracle.gd" || true
echo "=== Done ambiguous_base job ==="
