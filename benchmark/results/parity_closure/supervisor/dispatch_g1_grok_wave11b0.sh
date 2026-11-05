#!/usr/bin/env bash
# Cursor CLI grok-4.7-high-fast developer. Failover worker when GLM is rate limited.
# Foreground. The supervisor owns task state and the writer lock.
set -euo pipefail
REPO="${PARITY_REPO:-/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff}"
TASK="${REPO}/benchmark/results/parity_closure/g1_clonal_jc/wave11b0/G1_WAVE11B0_SHADOW.md"
OUT_DIR="${REPO}/benchmark/results/parity_closure/g1_clonal_jc/wave11b0"
MODEL="${PARITY_GROK_MODEL:-grok-4.7-high-fast}"
cd "${REPO}"
mkdir -p "${OUT_DIR}"
exec /hpcfs/fhome/caizhh/.local/bin/cursor-agent -p \
  --output-format json \
  --model "${MODEL}" \
  --trust \
  --auto-review \
  --workspace "${REPO}" \
  -- "$(cat "${TASK}")" \
  >"${OUT_DIR}/grok_wave11b0_stdout.log" \
  2>"${OUT_DIR}/grok_wave11b0_stderr.log"
