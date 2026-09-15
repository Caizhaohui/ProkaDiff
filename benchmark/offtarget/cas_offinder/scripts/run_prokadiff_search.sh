#!/usr/bin/env bash
# run_prokadiff_search.sh
# Runs ProkaDiff offtarget_scan on test reference
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HARNESS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT="$(cd "${HARNESS_DIR}/../../.." && pwd)"

REF="${1:-${ROOT}/testdata/offtarget/cas_offinder/reference.fa}"
OUTFILE="${2:-${HARNESS_DIR}/prokadiff.tsv}"
SPACER="${3:-GAGTCCGAGCAGAAGAAGAA}"
PAM="${4:-NGG}"
MAX_MM="${5:-5}"

echo "[run_prokadiff_search] Building and running offtarget_scan example..."
cargo run --release -p prokadiff-offtarget --example offtarget_scan -- \
  --ref "${REF}" \
  --spacer "${SPACER}" \
  --pam "${PAM}" \
  --editor cas9 \
  --max-mismatches "${MAX_MM}" \
  --output "${OUTFILE}"

echo "[run_prokadiff_search] Output saved to ${OUTFILE}"
