#!/usr/bin/env bash
# Read-only diagnosis of the 2387259 candidate-junction second pass.
# Do not run on a login node. Submit with:
#   mkdir -p benchmark/logs
#   sbatch scripts/submit_inspect_jc_second_pass.sh
#
#SBATCH -p qcpu_18i
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -J pd_jc2p
#SBATCH -o benchmark/logs/%x-%j.out
#SBATCH -e benchmark/logs/%x-%j.err

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "error: submit this script with sbatch (partition qcpu_18i); do not run it on a login node." >&2
  exit 1
fi

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"

JOB="${1:-2387259}"
BASE="${ROOT}/benchmark/results/clonal_${JOB}"
GD="${2:-${BASE}/breseq/output/output.gd}"
JFA="${3:-${BASE}/prokadiff/work/jc/junctions.fa}"
JSAM="${4:-${BASE}/prokadiff/work/jc/aligned.sam}"
REF="${5:-${ROOT}/testdata/layer2/clonal/Clonal_Sample/REL606.gbk}"

for f in "${GD}" "${JFA}" "${JSAM}" "${REF}"; do
  if [[ ! -f "${f}" ]]; then
    echo "error: missing ${f}" >&2
    exit 1
  fi
done
if ! command -v python3 >/dev/null 2>&1; then
  echo "error: python3 not found on PATH" >&2
  exit 1
fi

echo "host=$(hostname) job=${SLURM_JOB_ID} diagnosed=${JOB}"
ls -lh "${GD}" "${JFA}" "${JSAM}" "${REF}"
python3 "${ROOT}/scripts/inspect_jc_second_pass.py" \
  --gd "${GD}" --jfa "${JFA}" --ref "${REF}" --jsam "${JSAM}"
