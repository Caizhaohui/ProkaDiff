#!/usr/bin/env bash
# Wrapper — delegates to parity_closure_supervisor.py
set -euo pipefail
REPO="${PARITY_REPO:-/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff}"
PY="${REPO}/benchmark/results/parity_closure/supervisor/parity_closure_supervisor.py"
exec python3 "${PY}" --repo "${REPO}" "$@"
