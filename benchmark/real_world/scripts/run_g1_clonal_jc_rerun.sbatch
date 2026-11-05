#!/usr/bin/env bash
# G1 parity-closure Clonal JC re-run (ProkaDiff only vs frozen breseq GD).
# Does NOT overwrite benchmark/results/clonal_2673729 or verified/*.md.
#SBATCH -p qcpu_18i
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH -t 08:00:00
#SBATCH -J g1_clonal_jc
#SBATCH -o benchmark/logs/%x-%j.out
#SBATCH -e benchmark/logs/%x-%j.err

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "error: must run under Slurm" >&2
  exit 1
fi

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
cd "${ROOT}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
OUT="${ROOT}/benchmark/results/parity_closure/g1_clonal_jc/slurm/${SLURM_JOB_ID}"
FROZEN_BRESEQ="${ROOT}/benchmark/results/parity_closure/g1_clonal_jc/frozen_source/breseq_output.gd"
SAMPLE_ROOT="${ROOT}/testdata/layer2/clonal/Clonal_Sample"
mkdir -p "${OUT}/prokadiff" "${ROOT}/benchmark/logs"

CONDA_ENV="${PROKDIFF_CONDA_ENV:-/hpcfs/fhome/caizhh/.conda/envs/prokadiff}"
export PATH="${CONDA_ENV}/bin:${PATH:-}"
CARGO_TARGET_DIR="${ROOT}/target/slurm-g1-${SLURM_JOB_ID}"
export CARGO_TARGET_DIR

{
  echo "job=${SLURM_JOB_ID}"
  echo "node=$(hostname)"
  echo "date=$(date -Is)"
  echo "git=$(git rev-parse HEAD)"
  echo "threads=${THREADS}"
  breseq --version | head -1 || true
  bowtie2 --version | head -1 || true
} | tee "${OUT}/run_header.txt"

test -f "${FROZEN_BRESEQ}"
test -d "${SAMPLE_ROOT}"

mapfile -t FQS < <(find "${SAMPLE_ROOT}" -type f \( -name '*.fastq' -o -name '*.fastq.gz' \) | sort)
mapfile -t REFS < <(find "${SAMPLE_ROOT}" -type f \( -name '*.gbk' -o -name '*.gb' \) | sort)
REF="${REFS[0]}"
FQ_ARGS=()
for f in "${FQS[@]}"; do
  FQ_ARGS+=(--fastq "${f}")
done

cargo build --release -p prokadiff 2>&1 | tee "${OUT}/cargo_build.log"
BIN="${CARGO_TARGET_DIR}/release/prokadiff"

/usr/bin/time -v "${BIN}" evidence \
  --ref "${REF}" \
  "${FQ_ARGS[@]}" \
  --threads "${THREADS}" \
  --outdir "${OUT}/prokadiff" \
  --keep-bam \
  >"${OUT}/prokadiff.stdout" 2>"${OUT}/prokadiff.time"

python3 "${ROOT}/scripts/parity/jc_acceptance_gate.py" \
  --prokadiff-gd "${OUT}/prokadiff/output.gd" \
  --breseq-gd "${FROZEN_BRESEQ}" \
  --outdir "${OUT}/gate_strict" \
  --strict-no-whitelist \
  --normalize-contig-version \
  | tee "${OUT}/gate_strict_summary.txt"
GATE_RC=${PIPESTATUS[0]}

python3 "${ROOT}/scripts/layer2_gd_compare.py" \
  --jobout "${OUT}" \
  --workload clonal_g1 \
  --threads "${THREADS}" \
  --prokadiff-gd "${OUT}/prokadiff/output.gd" \
  --breseq-gd "${FROZEN_BRESEQ}" \
  --csv "${OUT}/layer2.csv" \
  | tee "${OUT}/layer2_summary.txt" || true

echo "formal_jc_gate_exit=${GATE_RC}" | tee -a "${OUT}/run_header.txt"
exit "${GATE_RC}"
