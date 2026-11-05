#!/usr/bin/env bash
# Dispatch Grok XHigh Developer for G1 wave-8 (single writer).
set -euo pipefail
REPO="${PARITY_REPO:-/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff}"
TASK="${REPO}/benchmark/results/parity_closure/g1_clonal_jc/engine/G1_GROK_WAVE8.md"
LOG_DIR="${REPO}/benchmark/results/parity_closure/g1_clonal_jc/engine"
STATE="${REPO}/benchmark/results/parity_closure/state/parity_closure_state.json"
OUT="${LOG_DIR}/grok_wave8_stdout.txt"
ERR="${LOG_DIR}/grok_wave8_stderr.log"
MODEL="${PARITY_GROK_MODEL:-grok-4.7-xhigh-fast}"

cd "${REPO}"
PROMPT="$(cat "${TASK}")"
agent --model "${MODEL}" --output-format json --force -p "${PROMPT}" >"${OUT}" 2>"${ERR}"
EC=$?
python3 - <<PY
import json
from pathlib import Path
from datetime import datetime, timezone
def iso():
    return datetime.now(timezone.utc).astimezone().strftime("%Y-%m-%dT%H:%M:%S%z")
p=Path("${STATE}")
s=json.loads(p.read_text())
s["worker_classification"]="WORKER_COMPLETED" if ${EC}==0 else "WORKER_FAILED"
s["developer_model"]="${MODEL}"
s["worker_wave"]="wave-8-grok"
s["status"]="WAITING_FOR_SLURM" if ${EC}==0 else "NEEDS_REPAIR"
s["next_action"]="CONTROLLER_SBATCH_G1_RERUN" if ${EC}==0 else "REDISPATCH_WAVE8"
s["checkpoint"]=iso()
tmp=p.with_suffix(".tmp"); tmp.write_text(json.dumps(s,indent=2)+"\n"); tmp.replace(p)
PY
exit "${EC}"
