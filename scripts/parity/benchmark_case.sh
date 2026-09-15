#!/usr/bin/env bash
# scripts/parity/benchmark_case.sh
# End-to-end execution of breseq oracle + ProkaDiff engine + parity comparison for one case

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

DATASET=""
REF=""
R1=""
R2=""
OUTDIR=""
THREADS=8
SKIP_BRESEQ=false
BRESEQ_GD=""
PROKADIFF_GD=""
JC_TOL_BP=5

usage() {
  cat <<EOF
Usage: $0 --dataset <name> --ref <ref> --r1 <fq1> [--r2 <fq2>] --outdir <dir> [options]
Options:
  --threads <N>         Number of threads (default: 8)
  --breseq-gd <path>    Use precomputed breseq GD instead of running breseq
  --skip-breseq         Skip running breseq (requires --breseq-gd or existing breseq.gd in outdir)
  --jc-tol-bp <N>       Tolerant matching distance for JC in bp (default: 5)
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset) DATASET="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --r1) R1="$2"; shift 2 ;;
    --r2) R2="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --breseq-gd) BRESEQ_GD="$2"; shift 2 ;;
    --skip-breseq) SKIP_BRESEQ=true; shift 1 ;;
    --jc-tol-bp) JC_TOL_BP="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1" >&2; usage ;;
  esac
done

[[ -z "${DATASET}" || -z "${REF}" || -z "${R1}" || -z "${OUTDIR}" ]] && usage

mkdir -p "${OUTDIR}"

# 1. Run or link breseq
if [[ -n "${BRESEQ_GD}" && -f "${BRESEQ_GD}" ]]; then
  echo "[benchmark_case] Using provided breseq GD: ${BRESEQ_GD}"
  cp "${BRESEQ_GD}" "${OUTDIR}/breseq.gd"
elif [[ "${SKIP_BRESEQ}" == "true" && -f "${OUTDIR}/breseq.gd" ]]; then
  echo "[benchmark_case] Using existing breseq GD: ${OUTDIR}/breseq.gd"
else
  echo "[benchmark_case] Running breseq oracle..."
  BRESEQ_ARGS=("--ref" "${REF}" "--r1" "${R1}" "--outdir" "${OUTDIR}" "--threads" "${THREADS}")
  if [[ -n "${R2}" ]]; then
    BRESEQ_ARGS+=("--r2" "${R2}")
  fi
  "${SCRIPT_DIR}/run_breseq.sh" "${BRESEQ_ARGS[@]}"
fi

# 2. Run ProkaDiff
echo "[benchmark_case] Running ProkaDiff engine..."
PROK_ARGS=("--ref" "${REF}" "--r1" "${R1}" "--outdir" "${OUTDIR}" "--threads" "${THREADS}")
if [[ -n "${R2}" ]]; then
  PROK_ARGS+=("--r2" "${R2}")
fi
"${SCRIPT_DIR}/run_prokadiff.sh" "${PROK_ARGS[@]}"

# 3. Compare outputs
echo "[benchmark_case] Running parity comparison..."
python3 "${SCRIPT_DIR}/compare_gd.py" \
  --dataset "${DATASET}" \
  --prokadiff-gd "${OUTDIR}/prokadiff.gd" \
  --breseq-gd "${OUTDIR}/breseq.gd" \
  --outdir "${OUTDIR}" \
  --jc-tol-bp "${JC_TOL_BP}"

echo "[benchmark_case] Completed benchmark for ${DATASET} in ${OUTDIR}"
