#!/usr/bin/env bash
# Quick Layer-2 Clonal_Sample: ProkaDiff evidence only + compare vs saved breseq oracle.
# Submit with:
#   sbatch scripts/submit_quick_clonal.sh
#
#SBATCH -p qcpu_18i
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -J pd_quick_clonal
#SBATCH -o benchmark/logs/%x-%j.out
#SBATCH -e benchmark/logs/%x-%j.err

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "error: submit this script with sbatch (partition qcpu_18i); do not run it on a login node." >&2
  exit 1
fi

THREADS="${SLURM_CPUS_PER_TASK:-8}"
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"

CONDA_ENV="${PROKDIFF_CONDA_ENV:-/hpcfs/fhome/caizhh/.conda/envs/BactGenome}"
if [[ -f "${CONDA_ENV}/bin/activate" ]]; then
  # shellcheck disable=SC1091
  source "${CONDA_ENV}/bin/activate" "${CONDA_ENV}" 2>/dev/null || true
fi
export PATH="${CONDA_ENV}/bin:${PATH:-}"

DATA="${ROOT}/testdata/layer2/clonal"
SAMPLE_ROOT="${DATA}/Clonal_Sample"
OUTPUT_ROOT="${DATA}/Clonal_Output"

mapfile -t FQS < <(find "${SAMPLE_ROOT}" -type f \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \) | sort)
mapfile -t REFS < <(find "${SAMPLE_ROOT}" -type f \( -name '*.gbk' -o -name '*.gb' -o -name '*.fa' -o -name '*.fna' -o -name '*.fasta' \) | sort)
ORACLE_GD="$(find "${OUTPUT_ROOT}" -type f -name 'output.gd' | sort | head -n 1 || true)"
SAVED_BRESEQ_GD="${ROOT}/benchmark/results/clonal_2408807/breseq/output/output.gd"

REF="${REFS[0]}"
echo "reference: ${REF}"
echo "reads: ${FQS[*]}"
echo "saved breseq gd: ${SAVED_BRESEQ_GD}"
echo "official oracle gd: ${ORACLE_GD}"

PROKDIFF="${ROOT}/target/release/prokadiff"
JOBOUT="${ROOT}/benchmark/results/clonal_${SLURM_JOB_ID}"
mkdir -p "${JOBOUT}/prokadiff"
TIME=(/usr/bin/time -v)

FQ_ARGS=()
for f in "${FQS[@]}"; do
  FQ_ARGS+=(--fastq "${f}")
done

echo "=== ProkaDiff evidence (clonal with Stage 2 split-read candidate discovery) ==="
"${TIME[@]}" -o "${JOBOUT}/prokadiff.time" \
  "${PROKDIFF}" evidence \
    --ref "${REF}" \
    "${FQ_ARGS[@]}" \
    --threads "${THREADS}" \
    --outdir "${JOBOUT}/prokadiff" \
    --keep-bam \
  | tee "${JOBOUT}/prokadiff.stdout"

RUST_GD="${JOBOUT}/prokadiff/output.gd"
if [[ ! -f "${RUST_GD}" ]]; then
  echo "error: ProkaDiff did not write ${RUST_GD}" >&2
  exit 1
fi

(
  cd "${JOBOUT}"
  gdtools SUBTRACT "${RUST_GD}" "${SAVED_BRESEQ_GD}" || true
  mv -f output.gd subtract_rust_minus_breseq.gd 2>/dev/null || true
  gdtools SUBTRACT "${SAVED_BRESEQ_GD}" "${RUST_GD}" || true
  mv -f output.gd subtract_breseq_minus_rust.gd 2>/dev/null || true
)

python3 "${ROOT}/scripts/layer2_gd_compare.py" \
  --jobout "${JOBOUT}" \
  --workload clonal \
  --threads "${THREADS}" \
  --prokadiff-gd "${RUST_GD}" \
  --breseq-gd "${SAVED_BRESEQ_GD}" \
  --oracle-gd "${ORACLE_GD}" \
  --jc-tol-bp 5 || true

echo "done ${JOBOUT}"
