#!/usr/bin/env bash
# Wave 11A forensic dispatch. OpenCode may write only the four evidence files.
set -euo pipefail
REPO="${PARITY_REPO:-/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff}"
TASK="${REPO}/benchmark/results/parity_closure/g1_clonal_jc/wave11a/G1_WAVE11A_FORENSIC.md"
LOG_DIR="${REPO}/benchmark/results/parity_closure/g1_clonal_jc/wave11a"
STATE="${REPO}/benchmark/results/parity_closure/state/parity_closure_state.json"
OUT="${LOG_DIR}/opencode_wave11a_stdout.txt"
ERR="${LOG_DIR}/opencode_wave11a_stderr.log"
MODEL="${PARITY_OPENCODE_MODEL:-zhipuai-coding-plan/glm-5.3}"
LOCK="${REPO}/benchmark/results/parity_closure/state/writer.lock"
SNAP="${LOG_DIR}/pre_glm_rust.sha256"

cd "${REPO}"
mkdir -p "${LOG_DIR}"
find crates -name '*.rs' -print | sort | xargs sha256sum > "${SNAP}"

python3 - <<PY
import json
from datetime import datetime, timezone
from pathlib import Path
p=Path("${STATE}")
s=json.loads(p.read_text())
s["worker_classification"]="WORKER_RUNNING"
s["developer_model"]="${MODEL}"
s["worker_wave"]="wave-11a-forensic"
s["status"]="WAVE11A_RUNNING"
s["next_action"]="WAIT_FOR_WAVE11A_FORENSIC"
s["checkpoint"]=datetime.now(timezone.utc).astimezone().strftime("%Y-%m-%dT%H:%M:%S%z")
tmp=p.with_suffix(".tmp")
tmp.write_text(json.dumps(s, indent=2)+"\n")
tmp.replace(p)
PY

printf 'pid=%s\nholder=opencode-wave11a\nstarted=%s\nhost=%s\n' "$$" "$(date -Is)" "$(hostname)" > "${LOCK}"

PROMPT="$(cat "${TASK}")"
set +e
opencode run --standalone -m "${MODEL}" --auto --format json -- "${PROMPT}" >"${OUT}" 2>"${ERR}"
EC=$?
set -e

python3 - <<PY
import hashlib, json, subprocess
from datetime import datetime, timezone
from pathlib import Path
repo = Path("${REPO}")
snap = Path("${SNAP}")
changed = []
current = {}
for line in subprocess.check_output(["bash","-lc","find crates -name '*.rs' -print | sort | xargs sha256sum"], cwd=repo, text=True).splitlines():
    digest, path = line.split(None, 1)
    current[path] = digest
before = {}
for line in snap.read_text().splitlines():
    digest, path = line.split(None, 1)
    before[path] = digest
for path in sorted(set(before) | set(current)):
    if before.get(path) != current.get(path):
        changed.append(path)
p=Path("${STATE}")
s=json.loads(p.read_text())
s["developer_model"]="${MODEL}"
s["worker_wave"]="wave-11a-forensic"
s["checkpoint"]=datetime.now(timezone.utc).astimezone().strftime("%Y-%m-%dT%H:%M:%S%z")
s["wave11a_exit_code"]=int("${EC}")
s["wave11a_production_drift"]=changed
if changed:
    s["worker_classification"]="WORKER_SCOPE_VIOLATION"
    s["status"]="WAVE11A_SCOPE_VIOLATION"
    s["next_action"]="CONTROLLER_REVIEW_SCOPE_VIOLATION"
elif int("${EC}")==0:
    s["worker_classification"]="WORKER_COMPLETED"
    s["status"]="WAVE11A_AWAITING_CODEX"
    s["next_action"]="DISPATCH_CODEX_WAVE11A_REVIEW"
else:
    s["worker_classification"]="WORKER_FAILED"
    s["status"]="WAVE11A_FAILED"
    s["next_action"]="CONTROLLER_DIAGNOSE_WAVE11A"
tmp=p.with_suffix(".tmp")
tmp.write_text(json.dumps(s, indent=2)+"\n")
tmp.replace(p)
lock=Path("${LOCK}")
if lock.exists():
    lock.unlink()
print("WAVE11A_EXIT", ${EC}, "DRIFT", len(changed))
PY
exit "${EC}"
