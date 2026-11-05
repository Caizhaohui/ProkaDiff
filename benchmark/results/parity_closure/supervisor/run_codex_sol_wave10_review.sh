#!/usr/bin/env bash
# Run Codex Sol independent review for Wave 10 (read-only audit).
set -euo pipefail
REPO="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff"
REVIEW_DIR="${REPO}/benchmark/results/parity_closure/g1_clonal_jc/review"
mkdir -p "${REVIEW_DIR}"
OUT_FILE="${REVIEW_DIR}/codex_sol_wave10_review.md"
ERR_FILE="${REVIEW_DIR}/codex_sol_wave10_review.err"

PROMPT="TASK_ID: G1_WAVE10_INDEPENDENT_REVIEW
Role: Independent Reviewer (Codex gpt-6.1-sol). Read-only code/algorithm audit. Implementer ≠ reviewer.
Review the Wave 10 code changes in ProkaDiff:
1. Examine git diff of crates/prokadiff-evidence/src/engine/jc_cluster.rs and tests.rs for Wave 10.
2. Read FIX_NOTES_OPENCODE_WAVE10.md.
3. Read benchmark/results/parity_closure/g1_clonal_jc/slurm/3066886/gate_strict_summary.txt and gate_strict/jc_acceptance_reconcile.tsv.
Audit objectives:
- Verify that your previous [P2] finding on interleaved-copy sort order was correctly and completely addressed.
- Verify the new deterministic unit test interleaved_repeat_copy_registers_merge_same_copy_cleanly.
- Audit the FP MOB @ 1269590 findings and mechanism described in Target C. Confirm whether MOB_COVERAGE_5_OF_5 != MOB_PARITY_PASS is correctly upheld.
- Confirm the strict gate outcome: matched=10, over=24, under=8, 5/5 oracle MOB exact, extra FP MOB persists.
- Provide clear scientific audit verdict: APPROVE_CODE, REJECT_CODE, or NEEDS_REVISION, and recommendations for Wave 11.
Format your final response in clear markdown and write it to ${OUT_FILE}."

echo "Starting Codex Sol Wave 10 review..."
cd "${REPO}"
codex exec -m gpt-6.1-sol -s read-only --ephemeral \
  -c 'model_reasoning_effort="high"' \
  -o "${OUT_FILE}" \
  "${PROMPT}" >"${REVIEW_DIR}/codex_sol_wave10_stdout.log" 2>"${ERR_FILE}"
EC=$?
echo "Codex Sol finished with exit code ${EC}"
exit "${EC}"
