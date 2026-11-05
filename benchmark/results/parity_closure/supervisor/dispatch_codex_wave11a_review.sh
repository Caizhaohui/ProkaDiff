#!/usr/bin/env bash
# Read-only Codex Sol review of Wave 11A forensic evidence.
set -euo pipefail
REPO="${PARITY_REPO:-/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff}"
REVIEW_DIR="${REPO}/benchmark/results/parity_closure/g1_clonal_jc/review"
PROMPT_FILE="${REPO}/benchmark/results/parity_closure/g1_clonal_jc/wave11a/G1_WAVE11A_CODEX_REVIEW.md"
STATE="${REPO}/benchmark/results/parity_closure/state/parity_closure_state.json"
LOCK="${REPO}/benchmark/results/parity_closure/state/writer.lock"
OUT="${REVIEW_DIR}/codex_sol_wave11a_review.md"
LOG="${REVIEW_DIR}/codex_sol_wave11a_stdout.log"
ERR="${REVIEW_DIR}/codex_sol_wave11a_review.err"

cd "${REPO}"
mkdir -p "${REVIEW_DIR}"
printf 'pid=%s\nholder=codex-wave11a-review\nstarted=%s\nhost=%s\n' "$$" "$(date -Is)" "$(hostname)" > "${LOCK}"
PROMPT="$(cat "${PROMPT_FILE}")"
set +e
codex exec -m gpt-6.1-sol -s read-only --ephemeral \
  -c 'model_reasoning_effort="high"' \
  -o "${OUT}" \
  "${PROMPT}" >"${LOG}" 2>"${ERR}"
EC=$?
set -e
python3 - <<PY
import json
from datetime import datetime, timezone
from pathlib import Path
p=Path("${STATE}")
s=json.loads(p.read_text())
s["checkpoint"]=datetime.now(timezone.utc).astimezone().strftime("%Y-%m-%dT%H:%M:%S%z")
s["wave11a_codex_exit_code"]=int("${EC}")
s["reviewer"]={
    "role": "Independent Reviewer",
    "model": "gpt-6.1-sol",
    "task": "G1_WAVE11A_INDEPENDENT_REVIEW",
    "review_file": "benchmark/results/parity_closure/g1_clonal_jc/review/codex_sol_wave11a_review.md",
    "exit_code": int("${EC}"),
}
if int("${EC}")==0 and Path("${OUT}").exists() and Path("${OUT}").stat().st_size>0:
    s["status"]="WAVE11A_CODEX_RETURNED"
    s["next_action"]="CONTROLLER_INTERPRET_WAVE11A_REVIEW"
    s["worker_classification"]="REVIEW_COMPLETED"
else:
    s["status"]="WAVE11A_CODEX_FAILED"
    s["next_action"]="CONTROLLER_DIAGNOSE_CODEX_WAVE11A"
    s["worker_classification"]="REVIEW_FAILED"
tmp=p.with_suffix(".tmp")
tmp.write_text(json.dumps(s, indent=2)+"\n")
tmp.replace(p)
lock=Path("${LOCK}")
if lock.exists():
    lock.unlink()
print("CODEX_WAVE11A_EXIT", ${EC}, "OUT_BYTES", Path("${OUT}").stat().st_size if Path("${OUT}").exists() else 0)
PY
exit "${EC}"
