#!/usr/bin/env bash
# run_cas_offinder.sh
# Runner for Cas-OFFinder 2.x oracle against reference.fa
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HARNESS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT="$(cd "${HARNESS_DIR}/../../.." && pwd)"

REF="${1:-${ROOT}/testdata/offtarget/cas_offinder/reference.fa}"
OUTFILE="${2:-${HARNESS_DIR}/cas_offinder.raw.tsv}"
DEVICE="${3:-C}" # 'C' for CPU

if ! command -v cas-offinder &>/dev/null; then
  echo "[run_cas_offinder] ERROR: cas-offinder executable not found in PATH." >&2
  echo "[run_cas_offinder] Status: BLOCKED_EXTERNAL_DEPENDENCY" >&2
  exit 127
fi

INPUT_TXT="${HARNESS_DIR}/input.txt"
cat <<EOF > "${INPUT_TXT}"
${REF}
NNNNNNNNNNNNNNNNNNNNNGG 5
GAGTCCGAGCAGAAGAAGAA 5
EOF

echo "[run_cas_offinder] Running: cas-offinder ${INPUT_TXT} ${DEVICE} ${OUTFILE}"
cas-offinder "${INPUT_TXT}" "${DEVICE}" "${OUTFILE}"
echo "[run_cas_offinder] Output saved to ${OUTFILE}"
