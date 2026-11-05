#!/usr/bin/env bash
# Run Codex Sol independent review for Wave 9 (read-only audit).
set -euo pipefail
REPO="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff"
REVIEW_DIR="${REPO}/benchmark/results/parity_closure/g1_clonal_jc/review"
mkdir -p "${REVIEW_DIR}"
OUT_FILE="${REVIEW_DIR}/codex_sol_wave9_review.md"
ERR_FILE="${REVIEW_DIR}/codex_sol_wave9_review.err"

PROMPT="TASK_ID: G1_WAVE9_INDEPENDENT_REVIEW
Role: Independent Reviewer (Codex gpt-6.1-sol). Read-only code/algorithm audit. Implementer ≠ reviewer.
Review the Wave 9 code changes in ProkaDiff:
1. Examine git diff of crates/prokadiff-evidence/src/engine/jc_cluster.rs and emit.rs.
2. Read FIX_NOTES_OPENCODE_WAVE9.md.
3. Read benchmark/results/parity_closure/g1_clonal_jc/slurm/3065896/gate_strict_summary.txt and gate_strict/jc_acceptance_reconcile.tsv.
Audit objectives:
- Verify clean-room compliance vs breseq: check that NO coordinates (e.g. 1821525, 16972, 3652533) or BR ids are hardcoded in engine logic.
- Verify algorithm correctness of Fix 1 (same-copy repeat_span clustering) and Fix 2 (support tie guard on MOB weak placement).
- Confirm the strict gate outcome: matched=10, over=26, under=8, MOB 5/5 exact. Confirm that Repeat Copy Equivalence 16/18 is NOT G1 PASS.
- Provide clear scientific audit verdict: APPROVE_CODE, REJECT_CODE, or NEEDS_REVISION.
Format your final response in clear markdown and write it to ${OUT_FILE} or output it directly."

echo "Starting Codex Sol review..."
cd "${REPO}"
codex exec -m gpt-6.1-sol -s read-only --ephemeral \
  -c 'model_reasoning_effort="high"' \
  -o "${OUT_FILE}" \
  "${PROMPT}" >"${REVIEW_DIR}/codex_sol_wave9_stdout.log" 2>"${ERR_FILE}"
EC=$?
echo "Codex Sol finished with exit code ${EC}"
exit "${EC}"
