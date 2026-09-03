#!/usr/bin/env bash
# Read-only BAM inspection for the synth_is_mob zero-JC diagnosis.
# Do not run on a login node. Submit with:
#   mkdir -p benchmark/logs
#   sbatch scripts/submit_inspect_is_mob_bam.sh [aligned.bam] [ref.fa]
#
#SBATCH -p qcpu_18i
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -t 00:10:00
#SBATCH -J pd_jc_inspect
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

BAM="${1:-${ROOT}/benchmark/results/synth_is_mob_2371412/prokadiff/aligned.bam}"
REF="${2:-${ROOT}/testdata/generated/synth_is_mob/ref.fa}"
LABEL="$(basename "$(dirname "$(dirname "${BAM}")")")/$(basename "$(dirname "${BAM}")")"

if [[ ! -f "${BAM}" ]]; then
  echo "error: BAM not found: ${BAM}" >&2
  exit 1
fi
if [[ ! -f "${REF}" ]]; then
  echo "error: reference not found: ${REF}" >&2
  exit 1
fi
if ! command -v samtools >/dev/null 2>&1; then
  echo "error: samtools not found on PATH" >&2
  exit 1
fi

echo "inspect ${BAM} (ref ${REF})"
samtools view "${BAM}" | python3 "${ROOT}/scripts/inspect_bam_jc.py" --label "${LABEL}" --ref "${REF}"
