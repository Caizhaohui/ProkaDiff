#!/usr/bin/env bash
# Layer-2 Clonal_Sample: ProkaDiff evidence vs official Clonal_Output.gd AND same-node breseq.
# Engine single-sample only (not product dual-sample). Submit with:
#   mkdir -p benchmark/logs testdata/layer2
#   sbatch scripts/submit_parity_clonal.sh
#
#SBATCH -p qcpu_18i
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH -t 24:00:00
#SBATCH -J pd_clonal
#SBATCH -o benchmark/logs/%x-%j.out
#SBATCH -e benchmark/logs/%x-%j.err

set -euo pipefail

time_field() {
  local label="$1"
  local path="$2"

  awk -v label="${label}" '
    {
      line = $0
      sub(/^[[:space:]]+/, "", line)
      if (index(line, label ":") == 1) {
        value = substr(line, length(label) + 2)
        sub(/^[[:space:]]+/, "", value)
        print value
        exit
      }
    }
  ' "${path}"
}

duration_to_seconds() {
  local duration="$1"

  if [[ ! "${duration}" =~ ^([0-9]+:)?[0-9]+:[0-9]+(\.[0-9]+)?$ ]]; then
    echo "error: unsupported /usr/bin/time wall-clock value: ${duration}" >&2
    return 1
  fi

  awk -v duration="${duration}" '
    BEGIN {
      count = split(duration, fields, ":")
      if (count == 2) {
        seconds = fields[1] * 60 + fields[2]
      } else if (count == 3) {
        seconds = fields[1] * 3600 + fields[2] * 60 + fields[3]
      } else {
        exit 1
      }
      printf "%.3f\n", seconds
    }
  '
}

manifest_value() {
  local key="$1"
  local path="$2"

  awk -F= -v key="${key}" '$1 == key { sub(/\r$/, "", $2); print $2; exit }' "${path}"
}

version_matches_pin() {
  local actual="$1"
  local expected="$2"
  local expected_pattern="${expected//./\\.}"

  [[ "${actual}" =~ (^|[^0-9.])${expected_pattern}([^0-9.]|$) ]]
}

write_verified_record() {
  if [[ "$#" -ne 10 ]]; then
    echo "error: write_verified_record requires jobout, record, job id, partition, node, threads, three versions, and a version manifest." >&2
    return 64
  fi

  local jobout="$1"
  local record="$2"
  local job_id="$3"
  local partition="$4"
  local node="$5"
  local threads="$6"
  local prokadiff_version="$7"
  local breseq_version="$8"
  local bowtie2_version="$9"
  local version_manifest="${10}"
  local expected_breseq
  local expected_bowtie2
  local prokadiff_time="${jobout}/prokadiff.time"
  local breseq_time="${jobout}/breseq.time"
  local prokadiff_wall
  local prokadiff_rss
  local breseq_wall
  local breseq_rss
  local prokadiff_wall_s
  local breseq_wall_s

  if [[ "${partition}" != "qcpu_18i" ]]; then
    echo "error: cannot publish verified evidence: expected SLURM_JOB_PARTITION=qcpu_18i, got ${partition:-<unset>}." >&2
    return 1
  fi
  if [[ ! -f "${version_manifest}" ]]; then
    echo "error: cannot publish verified evidence: missing version manifest ${version_manifest}." >&2
    return 1
  fi

  expected_breseq="$(manifest_value breseq "${version_manifest}")"
  expected_bowtie2="$(manifest_value bowtie2 "${version_manifest}")"
  if [[ ! "${expected_breseq}" =~ ^0\.40\.[0-9]+$ ]]; then
    echo "error: cannot publish verified evidence: manifest breseq pin must be 0.40.x, got ${expected_breseq:-<unset>}." >&2
    return 1
  fi
  if [[ ! "${expected_bowtie2}" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "error: cannot publish verified evidence: manifest bowtie2 pin is missing or invalid: ${expected_bowtie2:-<unset>}." >&2
    return 1
  fi
  if ! version_matches_pin "${breseq_version}" "${expected_breseq}"; then
    echo "error: cannot publish verified evidence: breseq version '${breseq_version}' does not match pinned ${expected_breseq}." >&2
    return 1
  fi
  if ! version_matches_pin "${bowtie2_version}" "${expected_bowtie2}"; then
    echo "error: cannot publish verified evidence: bowtie2 version '${bowtie2_version}' does not match pinned ${expected_bowtie2}." >&2
    return 1
  fi

  for path in "${prokadiff_time}" "${breseq_time}"; do
    if [[ ! -f "${path}" ]]; then
      echo "error: cannot publish verified evidence: missing timing file ${path}" >&2
      return 1
    fi
  done

  prokadiff_wall="$(time_field 'Elapsed (wall clock) time (h:mm:ss or m:ss)' "${prokadiff_time}")"
  prokadiff_rss="$(time_field 'Maximum resident set size (kbytes)' "${prokadiff_time}")"
  breseq_wall="$(time_field 'Elapsed (wall clock) time (h:mm:ss or m:ss)' "${breseq_time}")"
  breseq_rss="$(time_field 'Maximum resident set size (kbytes)' "${breseq_time}")"

  if [[ -z "${prokadiff_wall}" || -z "${prokadiff_rss}" || -z "${breseq_wall}" || -z "${breseq_rss}" ]]; then
    echo "error: cannot publish verified evidence: both tools need wall-clock and peak RSS fields." >&2
    return 1
  fi
  if [[ ! "${prokadiff_rss}" =~ ^[0-9]+$ || ! "${breseq_rss}" =~ ^[0-9]+$ ]]; then
    echo "error: cannot publish verified evidence: peak RSS must be an integer number of kbytes." >&2
    return 1
  fi
  if ! prokadiff_wall_s="$(duration_to_seconds "${prokadiff_wall}")"; then
    return 1
  fi
  if ! breseq_wall_s="$(duration_to_seconds "${breseq_wall}")"; then
    return 1
  fi

  mkdir -p "$(dirname "${record}")"
  cat >"${record}" <<EOF
# Verified Clonal parity evidence

- Slurm job: \`${job_id}\`
- Partition: \`${partition}\`
- Node: \`${node}\`
- Threads: \`${threads}\`
- Workload: Methods Mol. Biol. 2014 Clonal_Sample
- Comparator: \`scripts/layer2_gd_compare.py\` exited 0 after the bidirectional comparisons.
- Versions: ProkaDiff \`${prokadiff_version}\`; breseq \`${breseq_version}\`; bowtie2 \`${bowtie2_version}\`.

| Tool | Wall clock (s) | Peak RSS (kB) |
| --- | ---: | ---: |
| ProkaDiff | ${prokadiff_wall_s} | ${prokadiff_rss} |
| breseq | ${breseq_wall_s} | ${breseq_rss} |

This record was generated by \`scripts/submit_parity_clonal.sh\` only after its comparator succeeded and both \`/usr/bin/time -v\` files contained wall-clock and peak-RSS fields. The associated raw job directory and CSV remain local and ignored.
EOF
}

if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  return 0
fi

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "error: submit this script with sbatch (partition qcpu_18i); do not run it on a login node." >&2
  exit 1
fi

THREADS="${SLURM_CPUS_PER_TASK:-8}"
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"
CARGO_TARGET_DIR="${ROOT}/target/slurm-${SLURM_JOB_ID}"
export CARGO_TARGET_DIR

CONDA_ENV="${PROKDIFF_CONDA_ENV:-/hpcfs/fhome/caizhh/.conda/envs/prokadiff}"
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
PROKDIFF="${CARGO_TARGET_DIR}/release/prokadiff"

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

VERIFIED_RECORD="${ROOT}/benchmark/results/verified/clonal_${SLURM_JOB_ID}.md"
write_verified_record \
  "${JOBOUT}" \
  "${VERIFIED_RECORD}" \
  "${SLURM_JOB_ID}" \
  "${SLURM_JOB_PARTITION:-}" \
  "$(hostname -f 2>/dev/null || hostname)" \
  "${THREADS}" \
  "$("${PROKDIFF}" --version)" \
  "$(breseq --version | sed -n '1p')" \
  "$(bowtie2 --version | sed -n '1p')" \
  "${ROOT}/testdata/VERSIONS.txt"

echo "done ${JOBOUT}; verified evidence ${VERIFIED_RECORD}"
