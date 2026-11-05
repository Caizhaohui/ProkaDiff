#!/usr/bin/env bash
# PAUSE / RESUME / STOP / STATUS for Parity Closure Supervisor
set -euo pipefail
REPO="${PARITY_REPO:-/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff}"
CTRL="${REPO}/benchmark/results/parity_closure/state/supervisor_control.json"
STATE="${REPO}/benchmark/results/parity_closure/state/parity_closure_state.json"
SUP_LOCK="${REPO}/benchmark/results/parity_closure/state/supervisor.lock"
ts() { date -Is; }

cmd="${1:-status}"
case "${cmd}" in
  pause)
    python3 - <<PY
import json
from pathlib import Path
p=Path("${CTRL}")
d={"pause": True, "stop": False, "pause_reason": "operator_pause", "updated": "$(ts)"}
p.parent.mkdir(parents=True, exist_ok=True)
p.write_text(json.dumps(d, indent=2)+"\n")
print("PAUSED")
PY
    ;;
  resume)
    python3 - <<PY
import json
from pathlib import Path
p=Path("${CTRL}")
d={"pause": False, "stop": False, "pause_reason": "", "updated": "$(ts)"}
p.write_text(json.dumps(d, indent=2)+"\n")
print("RESUMED")
PY
    ;;
  stop)
    python3 - <<PY
import json
from pathlib import Path
p=Path("${CTRL}")
d={"pause": True, "stop": True, "pause_reason": "operator_stop", "updated": "$(ts)"}
p.write_text(json.dumps(d, indent=2)+"\n")
print("STOP_REQUESTED")
PY
    ;;
  status)
    echo "=== control ==="
    cat "${CTRL}" 2>/dev/null || echo "(no control file)"
    echo "=== state (short) ==="
    python3 -c 'import json,sys; s=json.load(open(sys.argv[1])); print(s.get("status"), s.get("next_action"), s.get("supervisor",{}).get("continuation_test_pass"))' "${STATE}" 2>/dev/null || true
    echo "=== supervisor lock ==="
    cat "${SUP_LOCK}" 2>/dev/null || echo "(not running)"
    ;;
  *)
    echo "usage: $0 {pause|resume|stop|status}" >&2
    exit 1
    ;;
esac
