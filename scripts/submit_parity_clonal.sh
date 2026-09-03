#!/usr/bin/env bash
# Layer-2 Clonal_Sample: ProkaDiff evidence vs official Clonal_Output.gd AND same-node breseq.
# Engine single-sample only (not product dual-sample). Submit with:
#   mkdir -p benchmark/logs testdata/layer2
#   sbatch scripts/submit_parity_clonal.sh
#
#SBATCH -p qcpu_18i
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH -t 24:00:00
#SBATCH -J pd_clonal
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
mkdir -p "${ROOT}/benchmark/logs" "${ROOT}/benchmark/results" "${DATA}"

missing=()
for bin in bowtie2 bowtie2-build breseq gdtools python3; do
  if ! command -v "${bin}" >/dev/null 2>&1; then
    missing+=("${bin}")
  fi
done
if ! command -v curl >/dev/null 2>&1 && ! command -v wget >/dev/null 2>&1; then
  missing+=("curl-or-wget")
fi
if (( ${#missing[@]} )); then
  echo "error: missing binaries: ${missing[*]}. Stop layer-2 clonal; do not invent leftover / wall / RSS." >&2
  exit 1
fi
if [[ ! -x /usr/bin/time ]]; then
  echo "error: /usr/bin/time not found. Stop layer-2; do not invent timing numbers." >&2
  exit 1
fi

fetch_url() {
  local url="$1"
  local dest="$2"
  if command -v curl >/dev/null 2>&1; then
    curl -fL --retry 3 --retry-delay 5 --connect-timeout 30 -o "${dest}.partial" "${url}"
  else
    wget -O "${dest}.partial" "${url}"
  fi
}

download() {
  local url="$1"
  local dest="$2"
  if [[ -s "${dest}" ]]; then
    echo "skip download (exists): ${dest}"
    return 0
  fi
  echo "downloading ${url} -> ${dest}"
  if ! fetch_url "${url}" "${dest}"; then
    rm -f "${dest}.partial"
    echo "error: download failed (no network or HTTP error): ${url}" >&2
    echo "Stop layer-2 clonal; do not invent leftover / wall / RSS." >&2
    exit 1
  fi
  mv "${dest}.partial" "${dest}"
}

WIKI_PUB="https://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing"
download "${WIKI_PUB}/Clonal_Sample.tgz" "${DATA}/Clonal_Sample.tgz"
download "${WIKI_PUB}/Clonal_Output.tgz" "${DATA}/Clonal_Output.tgz"

if [[ ! -d "${DATA}/Clonal_Sample" ]]; then
  tar -xzf "${DATA}/Clonal_Sample.tgz" -C "${DATA}"
fi
if [[ ! -d "${DATA}/Clonal_Output" ]]; then
  tar -xzf "${DATA}/Clonal_Output.tgz" -C "${DATA}"
fi

SAMPLE_ROOT="${DATA}/Clonal_Sample"
if [[ ! -d "${SAMPLE_ROOT}" ]]; then
  SAMPLE_ROOT="${DATA}"
fi
OUTPUT_ROOT="${DATA}/Clonal_Output"
if [[ ! -d "${OUTPUT_ROOT}" ]]; then
  OUTPUT_ROOT="${DATA}"
fi

mapfile -t FQS < <(find "${SAMPLE_ROOT}" -type f \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \) | sort)
mapfile -t REFS < <(find "${SAMPLE_ROOT}" -type f \( -name '*.gbk' -o -name '*.gb' -o -name '*.fa' -o -name '*.fna' -o -name '*.fasta' \) | sort)
ORACLE_GD="$(find "${OUTPUT_ROOT}" -type f -name 'output.gd' | sort | head -n 1 || true)"

if (( ${#FQS[@]} < 1 )); then
  echo "error: no FASTQ under ${SAMPLE_ROOT} after extract. Stop; do not invent results." >&2
  exit 1
fi
if (( ${#REFS[@]} < 1 )); then
  echo "error: no reference under ${SAMPLE_ROOT} after extract. Stop; do not invent results." >&2
  exit 1
fi
if [[ -z "${ORACLE_GD}" || ! -f "${ORACLE_GD}" ]]; then
  echo "error: Clonal_Output output.gd not found. Stop; do not invent results." >&2
  exit 1
fi

REF="${REFS[0]}"
echo "reference: ${REF}"
echo "reads: ${FQS[*]}"
echo "oracle gd: ${ORACLE_GD}"

echo "building release prokadiff (incremental; no-op if fresh)..."
cargo build --release -p prokadiff
PROKDIFF="${ROOT}/target/release/prokadiff"

JOBOUT="${ROOT}/benchmark/results/clonal_${SLURM_JOB_ID}"
mkdir -p "${JOBOUT}/prokadiff" "${JOBOUT}/breseq"
TIME=(/usr/bin/time -v)

FQ_ARGS=()
for f in "${FQS[@]}"; do
  FQ_ARGS+=(--fastq "${f}")
done

echo "=== ProkaDiff evidence (clonal) ==="
"${TIME[@]}" -o "${JOBOUT}/prokadiff.time" \
  "${PROKDIFF}" evidence \
    --ref "${REF}" \
    "${FQ_ARGS[@]}" \
    --threads "${THREADS}" \
    --outdir "${JOBOUT}/prokadiff" \
    --keep-bam \
  | tee "${JOBOUT}/prokadiff.stdout"

echo "=== breseq oracle (same node; env is breseq 0.40.2) ==="
"${TIME[@]}" -o "${JOBOUT}/breseq.time" \
  breseq -j "${THREADS}" -o "${JOBOUT}/breseq" \
    -r "${REF}" \
    "${FQS[@]}" \
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
  gdtools SUBTRACT "${RUST_GD}" "${ORACLE_GD}"
  mv -f output.gd subtract_rust_minus_official.gd
  gdtools SUBTRACT "${ORACLE_GD}" "${RUST_GD}"
  mv -f output.gd subtract_official_minus_rust.gd
)

python3 "${ROOT}/scripts/layer2_gd_compare.py" \
  --jobout "${JOBOUT}" \
  --workload clonal \
  --threads "${THREADS}" \
  --prokadiff-gd "${RUST_GD}" \
  --breseq-gd "${BRESEQ_GD}" \
  --oracle-gd "${ORACLE_GD}" \
  --csv "${ROOT}/benchmark/results/clonal.csv" \
  --jc-tol-bp 5

echo "done ${JOBOUT}"
