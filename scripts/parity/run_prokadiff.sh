#!/usr/bin/env bash
# scripts/parity/run_prokadiff.sh
# Run prokadiff evidence engine with time/RSS tracking and produce prokadiff.gd

set -euo pipefail

REF=""
R1=""
R2=""
OUTDIR=""
THREADS=8
PROKADIFF_BIN="${PROKADIFF_BIN:-}"
TIME_OUT=""
EXTRA_ARGS=()

usage() {
  cat <<EOF
Usage: $0 --ref <ref.fa|ref.gbk> --r1 <reads_1.fq> [--r2 <reads_2.fq>] --outdir <dir> [--threads <N>] [--prokadiff-bin <path>] [--time-out <time.txt>] [-- <extra-prokadiff-args>]
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref) REF="$2"; shift 2 ;;
    --r1) R1="$2"; shift 2 ;;
    --r2) R2="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --prokadiff-bin) PROKADIFF_BIN="$2"; shift 2 ;;
    --time-out) TIME_OUT="$2"; shift 2 ;;
    --) shift; EXTRA_ARGS=("$@"); break ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1" >&2; usage ;;
  esac
done

[[ -z "${REF}" || -z "${R1}" || -z "${OUTDIR}" ]] && usage

if [[ -z "${PROKADIFF_BIN}" ]]; then
  # Try release build, then debug build, then PATH
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
  if [[ -x "${REPO_ROOT}/target/release/prokadiff" ]]; then
    PROKADIFF_BIN="${REPO_ROOT}/target/release/prokadiff"
  elif [[ -x "${REPO_ROOT}/target/debug/prokadiff" ]]; then
    PROKADIFF_BIN="${REPO_ROOT}/target/debug/prokadiff"
  elif command -v prokadiff &>/dev/null; then
    PROKADIFF_BIN="$(command -v prokadiff)"
  else
    echo "[run_prokadiff] Error: prokadiff binary not found. Please cargo build first." >&2
    exit 1
  fi
fi

# Ensure bowtie2 is in PATH if in BactGenome env
if ! command -v bowtie2 &>/dev/null; then
  if [[ -x "/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin/bowtie2" ]]; then
    export PATH="/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin:${PATH}"
  fi
fi

mkdir -p "${OUTDIR}"
TIME_OUT="${TIME_OUT:-${OUTDIR}/prokadiff.time}"

CMD=("${PROKADIFF_BIN}" "evidence" "--ref" "${REF}" "--threads" "${THREADS}" "--outdir" "${OUTDIR}" "--fastq" "${R1}")
if [[ -n "${R2}" ]]; then
  CMD+=("--fastq" "${R2}")
fi
if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
  CMD+=("${EXTRA_ARGS[@]}")
fi

echo "[run_prokadiff] Running: ${CMD[*]}"
/usr/bin/time -v -o "${TIME_OUT}" "${CMD[@]}"

if [[ -f "${OUTDIR}/output.gd" ]]; then
  cp "${OUTDIR}/output.gd" "${OUTDIR}/prokadiff.gd"
  echo "[run_prokadiff] Wrote ${OUTDIR}/prokadiff.gd"
elif [[ -f "${OUTDIR}/prokadiff.gd" ]]; then
  echo "[run_prokadiff] Verified ${OUTDIR}/prokadiff.gd"
else
  echo "[run_prokadiff] Error: output.gd not produced by prokadiff evidence" >&2
  exit 1
fi
