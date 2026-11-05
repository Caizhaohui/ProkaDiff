#!/usr/bin/env bash
# Read-only Codex review. Runs only after Wave 11B0 and Wave 11B2 have completed.
set -euo pipefail
REPO="${PARITY_REPO:-/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff}"
REVIEW_DIR="${REPO}/benchmark/results/parity_closure/g1_clonal_jc/review"
OUT="${REVIEW_DIR}/codex_sol_wave11b_review.md"
LOG="${REVIEW_DIR}/codex_sol_wave11b_stdout.log"
ERR="${REVIEW_DIR}/codex_sol_wave11b_review.err"
mkdir -p "${REVIEW_DIR}"
cd "${REPO}"
PROMPT='TASK_ID: WAVE11D_CODEX_REVIEW
Role: Independent Reviewer, Codex gpt-6.1-sol. Read-only.
Read and follow ~/.config/ai-collaboration/POLICY.md before acting. It is not auto-loaded.
Do not edit crates/, docs/, scripts/, or the acceptance gate. Do not run cargo or sbatch.
Review benchmark/results/parity_closure/g1_clonal_jc/wave11b0r/ (the repaired analyzer and its tables) plus wave11b2/WAVE11B2_REGRESSION_REPORT.md against the job 3066886 SAM and FASTA.
Check that every eligible construct was available to the competition, that highest-score ties are preserved, that no coordinate-specific exception is used, and that JC 330, JC 331, the JC 326/321/343 family, and breseq 114 / JC 344 are represented as the SAM supports.
Confirm whether the two defects you raised in your previous Wave 11B review (elif side-order test, len(jc_union) < len(cons) NA heuristic) are actually fixed in wave11b0r, and recheck the numbers you could not verify last time.
Begin your answer with a single line of the form: WAVE11D_VERDICT: PASS or WAVE11D_VERDICT: NEEDS_REVISION.
Do not declare G1 accepted. Job 3072550 is not the acceptance baseline.'
exec /hpcfs/fhome/caizhh/.local/bin/codex exec -m gpt-6.1-sol -s read-only --ephemeral \
  -c 'model_reasoning_effort="high"' \
  -o "${OUT}" \
  "${PROMPT}" >"${LOG}" 2>"${ERR}"
