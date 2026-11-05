#!/usr/bin/env bash
# Foreground OpenCode worker. The supervisor owns state and the writer lock.
set -euo pipefail
REPO="${PARITY_REPO:-/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff}"
TASK="${REPO}/benchmark/results/parity_closure/g1_clonal_jc/wave11b0/G1_WAVE11B0_SHADOW.md"
LOG_DIR="${REPO}/benchmark/results/parity_closure/g1_clonal_jc/wave11b0"
MODEL="${PARITY_OPENCODE_MODEL:-zhipuai-coding-plan/glm-5.3}"
cd "${REPO}"
mkdir -p "${LOG_DIR}"
exec /hpcfs/fhome/caizhh/.opencode/bin/opencode run --standalone -m "${MODEL}" --auto --format json -- "$(cat "${TASK}")" \
  >"${LOG_DIR}/opencode_wave11b0_stdout.txt" \
  2>"${LOG_DIR}/opencode_wave11b0_stderr.log"
