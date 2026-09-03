#!/usr/bin/env bash
# Read-only Clonal BAM inspection: S>=14 clips near oracle JC windows.
# Defaults to clonal_2372151 outputs. Do not run on a login node. Submit with:
#   mkdir -p benchmark/logs
#   sbatch scripts/submit_inspect_clonal_bam.sh [aligned.bam] [oracle.gd]
#
#SBATCH -p qcpu_18i
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -J pd_clonal_bam
#SBATCH -o benchmark/logs/%x-%j.out
#SBATCH -e benchmark/logs/%x-%j.err

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "error: submit this script with sbatch (partition qcpu_18i); do not run it on a login node." >&2
  exit 1
fi

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"

CONDA_ENV="${PROKDIFF_CONDA_ENV:-/hpcfs/fhome/caizhh/.conda/envs/BactGenome}"
if [[ -f "${CONDA_ENV}/bin/activate" ]]; then
  # shellcheck disable=SC1091
  source "${CONDA_ENV}/bin/activate" "${CONDA_ENV}" 2>/dev/null || true
fi
export PATH="${CONDA_ENV}/bin:${PATH:-}"

BAM="${1:-${ROOT}/benchmark/results/clonal_2372151/prokadiff/aligned.bam}"
GD="${2:-${ROOT}/benchmark/results/clonal_2372151/breseq/output/output.gd}"
LABEL="$(basename "$(dirname "$(dirname "${BAM}")")")/$(basename "$(dirname "${BAM}")")"

if [[ ! -f "${BAM}" ]]; then
  echo "error: BAM not found: ${BAM}" >&2
  exit 1
fi
if [[ ! -f "${GD}" ]]; then
  echo "error: GD not found: ${GD}" >&2
  exit 1
fi
if ! command -v samtools >/dev/null 2>&1; then
  echo "error: samtools not found on PATH" >&2
  exit 1
fi
if ! command -v python3 >/dev/null 2>&1; then
  echo "error: python3 not found on PATH" >&2
  exit 1
fi

echo "host=$(hostname) job=${SLURM_JOB_ID}"
echo "inspect ${BAM}"
echo "oracle ${GD}"
ls -lh "${BAM}"
samtools view -H "${BAM}" | head -n 20
echo "--- counts ---"
samtools view -h "${BAM}" | python3 "${ROOT}/scripts/inspect_clonal_bam_jc.py" \
  --gd "${GD}" --label "${LABEL}"
