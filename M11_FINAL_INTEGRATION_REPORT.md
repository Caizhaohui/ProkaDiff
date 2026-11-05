# Milestone M11: Final Integration & Acceptance Report
## 亚克隆低频变异与混群审计模式最终验收报告 (Subclonal / Pool-seq Variant Audit Mode Final Report)

- **Milestone ID**: `M11`
- **Controller**: OpenCode CLI / GLM-5.3-Flash (acting Controller per user-authorized handoff 2026-11-05; former Controller Antigravity / Gemini 3.8 Flash Medium weekly quota exhausted mid-milestone)
- **Scientific Architect**: Claude Opus 5.5 (`claude-opus-5-5`, M11-T01, $0.205 verified spend)
- **Developer P1**: OpenCode CLI / GLM-5.3 (separate task identities per dispatch ledger; combined GLM quota tracked)
- **Independent Reviewer**: OpenAI Codex CLI / `gpt-6.1-sol` (Rounds 1–6, all read-only)
- **Governing Policy**: Multi-AI Collaboration V3.2 (`~/.config/ai-collaboration/POLICY.md`)
- **Execution Window**: 2026-11-05 08:56 – 2026-11-05 (controller handoff recovery audit → acceptance)
- **Slurm**: `qfree` (vnode1; jobs 3130686, 3130717); fallback interactive allocation `3095377` (`bnode2`, `qcpu_18i`) used only for the takeover recovery run with reasons logged
- **Final Verdict**: **`[PASS] Milestone M11 Round 6 Accepted`** — Codex Sol, no open P1/P2 findings (`.prokadiff/codex_m11_r6_review.txt`)
- **Baseline Commit**: `b30fbadfb666d11cfcccac605be742e3ed223a77` (user-authorized, local only)

---

## 1. Executive Summary

Milestone M11 delivers the production-ready **opt-in Subclonal / Pool-seq Variant Audit Mode** (`--subclonal`) for ProkaDiff:

1. **Low-Frequency Statistics Engine (`M11-T02`)**: strand-specific read counting, continuity-corrected Wilson score CIs (Newcombe 1998 Method 4), two-tailed Fisher exact strand-bias test (R `fisher.test` criterion, log-space with Stirling tail), upper-tail binomial error test with local ±50 bp background estimation and error floors (SNP 0.005 / INDEL 0.010), Benjamini–Hochberg FDR, underflow-protected binomial recurrence, u64-widened checked count arithmetic.
2. **Domain Model & Classifier (`M11-T03`)**: `SubclonalMetrics`/`SubclonalVariant` with dual-mode (consensus vs opt-in) parsing, parent-background differential (ΔAF with Newcombe mixed-score CI), and the 7-class deterministic priority classification (`EDITED_SPECIFIC`, `EDITED_ENRICHED`, `PRE_EXISTING_IN_STARTER`, `DEPLETED_IN_EDITED`, `AMBIGUOUS_STARTER_TRACE`, `INDETERMINATE_STARTER_DEPTH`, `NOT_SUBCLONAL`) with strict reportable filtering.
3. **Non-Prescriptive Reporting (`M11-T04`)**: 26-column `subclonal_audit.tsv` (documented in `docs/schema.md`) and `report.md` subclonal section with the mandatory verbatim uncertainty caveat, including zero-result audits; dynamic product publication (`product_output_names`) with stale-optional-product removal and full rollback on failure.
4. **Read-Level Allele Identity (`Rounds 2–4 hardening`)**: per-read pileup identity (`BaseObs.read_id`, `InsertionObs`, `DeletionObs` recorded post-normalization, pre-trim, one event per read per deletion run); per-read REF/ALT/OTHER classification with set-difference semantics; deletion spans structurally immune to homopolymer/end trimming; ≥3 bp deletions excluded from the authorized ≤2 bp scope.
5. **Production Titration Suite (`M11-T05`, `Rounds 3–5 hardening`)**: real-CLI-pipeline gradient (FASTQ → bowtie2 → pileup → extraction → evaluation) at 50%, 25%, 10%, 5% SNP AF plus 1 bp DEL (25%), 2 bp DEL (20%), INS `GT` (15%), and 0% negative control; exact depth/ALT-count assertions, ±0.02 AF tolerance, builder-path regression for the trimming counterexample.
6. **CLI & Validation (`Rounds 1–5 hardening`)**: `--subclonal`, `--min-subclonal-af` (0.01–0.50, default 0.05), `--min-subclonal-alt-reads` (default 3), `--max-strand-bias-pvalue` (default 0.001); `is_finite()` bounds guards with NaN regression tests; Rule 6 starter mandate preserved.

**Default consensus behavior is untouched**: subclonal processing is strictly opt-in-gated; all pre-existing tests pass unchanged; default products contain no subclonal content (provenance-hashed bundle, `subclonal_audit.tsv` verified absent).

## 2. Quality Gates Status Matrix (as amended, contract §6 + §7)

| Quality Gate | Requirement | Evidence | Status |
|---|---|---|---|
| **Gate A: Non-Regression** (revised basis §7) | Opt-in isolation; pre-existing tests unchanged & green; committed pre-M11 oracle parity; provenance-hashed default outputs; explicit residual limitation | User authorization record `.prokadiff/m11_gate_a_authorization.md`; ERRATUM in `.prokadiff/m11_r4_repair_validation.log`; `.prokadiff/m11_gate_a_default_outputs.txt` (10 products hashed, job 3130686); 256 oracle tests green; 494/0/3 workspace suite; mtime proof | **PASS** |
| **Gate B: Opt-in Precision** | Correct low-frequency calling within authorized scope (SNP/SUB/INS/DEL ≤2 bp), negative-control suppression | Per-read counting fixes verified R2–R5; production titration with exact counts; `test_production_read_titration_gradient_with_negative_control`; mixed/disjoint deletion tests; masked-anchor insertion tests | **PASS** |
| **Gate C: Scientific Language** | `audit_release_gates.py` strict pass; mandatory caveat; no causal/sensitivity overclaims | 3 gates PASS / 0 findings / exit 0 (Controller run + Codex independent rerun R2–R6); `EDITED_SPECIFIC` wording uncertainty-aware; zero-result caveat emission | **PASS** |
| **Gate D: HPC Execution** | Builds/tests/synthetic alignment on Slurm `qfree` | qfree vnode1 jobs 3130686 & 3130717 (14 GB, `RUST_TEST_THREADS=4`); SLURM_* env-conflict root cause + OOM diagnosis recorded; fallback runs disclosed | **PASS** |
| **Gate E: Independent Review** | Codex Sol (`gpt-6.1-sol`) independent sign-off | **Round 6: `[PASS] Milestone M11 Round 6 Accepted`** — no open findings (`.prokadiff/codex_m11_r6_review.txt`) | **PASS (ACCEPTED)** |

## 3. Multi-AI Review Trail (Rounds 1–6)

| Round | Verdict | Findings | Resolution |
|---|---|---|---|
| **R1** (09:43) | REJECTED | 9 blocking findings | Fixed pre-R2 (production integration, dynamic publication, authorized variant-type boundary, error floors, binomial underflow, checked adds, reportable filtering, CLI bounds, zero-result caveat) |
| **R2** (10:19) | REJECTED | 4×P1 + 2×P2 | Fixed 10:20–10:27 (insertion REF/ALT mutual exclusion; 2 bp DEL single event; starter detection decoupled from subclonal ceiling; production titration; NaN guards; stale-product cleanup). Session interrupted by Gemini quota exhaustion at the full-suite step |
| **R3** (≈11:20) | REJECTED | 2×P1 + 3×P2 | Repair round 1 (attempt 1/3): read-level pileup identity; per-read classification; production INS/DEL loci + exact counts; sha256 snapshot manifest; qfree partition reconciliation (DISP-M11-012) |
| **R4** | REJECTED | 1×P1 + 2×P2 | Repair round 2 (attempt 2/3): structural `DeletionObs` spans immune to trimming + builder-path counterexample regression; Gate A default-output hash bundle + 256 committed oracles; schema wording fix (DISP-M11-014) |
| **R5** | REJECTED | 1×P2 (Gate A basis only) | ERRATUM appended; pre-M11 comparison feasibility exhaustively investigated (impossible); **user authorized Gate A basis revision** (contract §7); baseline commit `b30fbad` (DISP-M11-016) |
| **R6** | **[PASS] Accepted** | none | Formal acceptance against amended contract; frozen-criteria preservation verified |

## 4. Deterministic Validation Evidence (final snapshot)

| Item | Result | Log |
|---|---|---|
| `cargo test --workspace` | **494 passed / 0 failed / 3 ignored**, exit 0 — qfree vnode1, job 3130717 | `.prokadiff/m11_r4_repair_validation.log` |
| `cargo clippy --workspace --all-targets -- -D warnings` | 0 warnings, exit 0 (qfree) | same |
| Python oracle suites (58 parity + 3 layer-2 + 195 release-gate) | **256 passed**, exit 0 | same |
| `audit_release_gates.py` | 3 gates PASS / 0 findings / exit 0 (Controller + Codex independent reruns) | same + `.prokadiff/codex_m11_r6_review.txt` |
| Default-mode output bundle | 10 products sha256-hashed with full provenance; `subclonal_audit.tsv` absent; schema V2 confirmed | `.prokadiff/m11_gate_a_default_outputs.txt` |
| Snapshot manifest | 162 tracked-file entries bound to HEAD `b30fbad`, manifest sha256 `6aae6c16…` | `.prokadiff/m11_snapshot_manifest.txt` |

## 5. Recorded Limitations (honesty contract)

1. **Gate A revised basis**: literal byte-for-byte equivalence to the M10-final binary is NOT established (no pre-M11 state existed; see contract §7 and authorization record). From M12 onward, content-addressed baselines are committed before new development, making literal byte-diff comparisons standard practice.
2. **Synthetic validation only**: titration tests use pseudorandom **non-repeating** synthetic references; this is NOT evidence of performance in repetitive genomic regions.
3. **AF 0.05 is a configurable threshold, not an experimental sensitivity floor**; detection limits depend on depth, strand balance, mapping quality, and local error rate. No "detects down to 5%" claims are made anywhere.
4. **Authorized low-frequency scope**: SNP, SUB, INS/DEL ≤2 bp only. Low-frequency MOB/JC/AMP calling is not implemented and remains out of scope.
5. **qfree nodes are memory-bounded** (16 CPU / 15.5 GB): full workspace validation requires `--cpus-per-task=4 --mem=14000` with `RUST_TEST_THREADS=4`; unbounded runs OOM (recorded).
6. **Three pre-existing ignored tests** remain ignored (untouched, pre-M11).

## 6. Preserved Invariants (verified by Round 6)

- Default single-clone consensus behavior — unchanged (opt-in gating verified).
- DV1 EventId semantics — untouched.
- M10 batch audit compatibility — `batch_e2e` 10/10 green.
- Mutation classification contracts — 138 classify tests green.
- G1 Clonal JC `FROZEN_HOLD` (`MATCHED=10, OVER=24, UNDER=8, PASS=false`) — untouched.
- Validation Registry `PENDING` states (PRJNA1088182, PRJNA884016) — untouched.
- Rule 6 starter mandate; Rule 9 noodles BAM I/O; Rule 14 read-only real data — upheld (no new dependencies added).

## 7. Milestone Status

**`AWAITING_USER_APPROVAL`** — technical work complete and independently accepted. The user (`FINAL_AUTHORIZATION_AUTHORITY`) retains push, tag, and release authority. Next roadmap milestone requires a fresh single-round user authorization per `STAGE_ROADMAP.md`.
