#!/usr/bin/env bash
# scripts/parity/run_breseq.sh
# Run breseq oracle with time/RSS tracking and produce breseq.gd

set -euo pipefail

REF=""
R1=""
R2=""
OUTDIR=""
THREADS=8
BRESEQ_BIN="${BRESEQ_BIN:-breseq}"
TIME_OUT=""

usage() {
  cat <<EOF
Usage: $0 --ref <ref.fa|ref.gbk> --r1 <reads_1.fq> [--r2 <reads_2.fq>] --outdir <dir> [--threads <N>] [--time-out <time.txt>]
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
    --breseq-bin) BRESEQ_BIN="$2"; shift 2 ;;
    --time-out) TIME_OUT="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1" >&2; usage ;;
  esac
done

[[ -z "${REF}" || -z "${R1}" || -z "${OUTDIR}" ]] && usage

# If breseq is not directly in PATH, check BactGenome conda env
if ! command -v "${BRESEQ_BIN}" &>/dev/null; then
  if [[ -x "/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin/breseq" ]]; then
    BRESEQ_BIN="/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin/breseq"
    export PATH="/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin:${PATH}"
  fi
fi

mkdir -p "${OUTDIR}"
TIME_OUT="${TIME_OUT:-${OUTDIR}/breseq.time}"
BRESEQ_RAW_DIR="${OUTDIR}/breseq_raw"
rm -rf "${BRESEQ_RAW_DIR}"
mkdir -p "${BRESEQ_RAW_DIR}"

CMD=("${BRESEQ_BIN}" "-r" "${REF}" "-j" "${THREADS}" "-o" "${BRESEQ_RAW_DIR}" "${R1}")
if [[ -n "${R2}" ]]; then
  CMD+=("${R2}")
fi

echo "[run_breseq] Running: ${CMD[*]}"
/usr/bin/time -v -o "${TIME_OUT}" "${CMD[@]}"

if [[ -f "${BRESEQ_RAW_DIR}/output/output.gd" ]]; then
  cp "${BRESEQ_RAW_DIR}/output/output.gd" "${OUTDIR}/breseq.gd"
  echo "[run_breseq] Wrote ${OUTDIR}/breseq.gd"
elif [[ -f "${BRESEQ_RAW_DIR}/output.gd" ]]; then
  cp "${BRESEQ_RAW_DIR}/output.gd" "${OUTDIR}/breseq.gd"
  echo "[run_breseq] Wrote ${OUTDIR}/breseq.gd"
else
  echo "[run_breseq] Error: output.gd not produced by breseq" >&2
  exit 1
fi
