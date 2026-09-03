#!/usr/bin/env bash
# Layer-1 synth_is_mob: ProkaDiff evidence vs breseq on qcpu_18i.
# Do not run on a login node. Submit with:
#   mkdir -p benchmark/logs
#   sbatch scripts/submit_parity_synth_is_mob.sh
#
#SBATCH -p qcpu_18i
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -J pd_is_mob
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

mkdir -p "${ROOT}/benchmark/logs" "${ROOT}/benchmark/results" "${ROOT}/testdata/generated"

missing=()
for bin in wgsim bowtie2 bowtie2-build breseq gdtools python3; do
  if ! command -v "${bin}" >/dev/null 2>&1; then
    missing+=("${bin}")
  fi
done
if (( ${#missing[@]} )); then
  echo "error: missing binaries: ${missing[*]}. Stop layer-1; do not invent parity numbers." >&2
  exit 1
fi

echo "building release prokadiff (incremental; no-op if fresh)..."
cargo build --release -p prokadiff
PROKDIFF="${ROOT}/target/release/prokadiff"

bash "${ROOT}/testdata/generate.sh" synth_is_mob "${ROOT}/testdata/generated/synth_is_mob"
GEN="${ROOT}/testdata/generated/synth_is_mob"
JOBOUT="${ROOT}/benchmark/results/synth_is_mob_${SLURM_JOB_ID}"
mkdir -p "${JOBOUT}/prokadiff" "${JOBOUT}/breseq"

TIME=(/usr/bin/time -v)

echo "=== ProkaDiff evidence ==="
"${TIME[@]}" -o "${JOBOUT}/prokadiff.time" \
  "${PROKDIFF}" evidence \
    --ref "${GEN}/ref.fa" \
    --fastq "${GEN}/edited_R1.fastq" \
    --fastq "${GEN}/edited_R2.fastq" \
    --threads "${THREADS}" \
    --outdir "${JOBOUT}/prokadiff" \
    --keep-bam \
  | tee "${JOBOUT}/prokadiff.stdout"

echo "=== breseq oracle ==="
"${TIME[@]}" -o "${JOBOUT}/breseq.time" \
  breseq -j "${THREADS}" -o "${JOBOUT}/breseq" \
    -r "${GEN}/ref.fa" \
    "${GEN}/edited_R1.fastq" \
    "${GEN}/edited_R2.fastq" \
  | tee "${JOBOUT}/breseq.stdout"

RUST_GD="${JOBOUT}/prokadiff/output.gd"
BRESEQ_GD="${JOBOUT}/breseq/output/output.gd"
if [[ ! -f "${RUST_GD}" ]]; then
  echo "error: ProkaDiff did not write ${RUST_GD}" >&2
  exit 1
fi
if [[ ! -f "${BRESEQ_GD}" ]]; then
  echo "error: breseq did not write ${BRESEQ_GD}" >&2
  exit 1
fi

(
  cd "${JOBOUT}"
  gdtools SUBTRACT "${RUST_GD}" "${BRESEQ_GD}"
  mv -f output.gd subtract_rust_minus_breseq.gd
  gdtools SUBTRACT "${BRESEQ_GD}" "${RUST_GD}"
  mv -f output.gd subtract_breseq_minus_rust.gd
)

python3 "${ROOT}/scripts/layer2_gd_compare.py" \
  --jobout "${JOBOUT}" \
  --workload synth_is_mob \
  --csv "${ROOT}/benchmark/results/synth_is_mob.csv" \
  --threads "${THREADS}" \
  --prokadiff-gd "${RUST_GD}" \
  --breseq-gd "${BRESEQ_GD}" \
  || true

echo "done ${JOBOUT} (do not invent leftover/wall/RSS if compare failed)"
