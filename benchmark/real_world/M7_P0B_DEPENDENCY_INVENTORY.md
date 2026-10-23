# M7-P0b Dependency Inventory

Provenance: ZCode 0.17.0, provider `account:bigmodel-individual-coding-plan`, model `GLM-5.3-Flash`, mode `plan` (Write and Edit disallowed), session `sess_ce6e2624-4ffc-4282-9273-efc45b171255`. Exit code 0. The analyst did not modify the worktree; the orchestrator saved this document from the read-only response.

MODEL: account:bigmodel-individual-coding-plan/GLM-5.3-Flash

## Git state

- Local HEAD: `70e299d845b8204d63c0ae3e6906fa173e8ab424` — matches expected (`70e299d fix(benchmark): quarantine invalid PRJNA884016 structural truth and enforce fail-closed invariants (M7-P0a)`).
- `origin/main`: `5c90e0568153a4befc573933d2f98bd3d47662b3` — matches expected. HEAD is one commit ahead of origin; zero push.
- Note: `M7_P0A_INVALID_TRUTH_REMEDIATION.md` §G records the corrective commit as `894a613a0e107577579737dc1d5bf5967f892c9f`, but HEAD carrying the identical message is `70e299d8…` — the commit was evidently amended/rebased after the dossier was written. Provenance discrepancy only; message and tree are consistent.
- Worktree dirty state matches the snapshot exactly: 1 modified tracked file (`benchmark/real_world/scripts/stage_m6_p1a_fastq.sh`) and 14 untracked M6-P2A files. No M7-P0b-relevant file is modified.
- M6 status `CLOSED_WITH_LIMITATIONS`, M7-P0a status `M7_P0A_INVALID_TRUTH_QUARANTINED` confirmed from dossier header and registry.

## Consumers

### 1. `benchmark/real_world/truth/prjna884016/SY22A1_truth.gd` (active truth, tracked)
- Reader(s): `compare_mob_truth.py:parse_mob_records`, `compare_variant_calls.py:parse_gd` (via `--truth-gd`), `check_validation_assets.sh:39-40` (awk MOB/JC count), `test_prjna884016_quarantine.py:test_legacy_corrupted_truth_records_rejected`.
- State: replaced with a zero-variant header. Evidence: `#=TITLE PRJNA884016 SY22A1 Negative Control Assembly Truth (0 Chromosomal Variants)` and `#=NOTE … Historical 4 MOB + 7 JC invalidated as 194.7 kb truncation artifact.`
- Can it treat 4 MOB / 7 JC as positive truth? **No.** 0 MOB / 0 JC records; asserted by gate (`-eq 0`) and pytest.
- Risk: none. **ALREADY_GUARDED.**

### 2. `benchmark/real_world/truth/prjna884016/mob_truth.tsv` and `jc_truth.tsv` (active, tracked)
- Reader(s): `check_validation_assets.sh:41-42` (`NR > 1` row count `-eq 0`), `test_prjna884016_quarantine.py:134-144` (`"must contain only the header line"`).
- State: header only, 0 data rows (85 B and 51 B on disk).
- Can they load the invalidated set? **No.** **ALREADY_GUARDED.**

### 3. `benchmark/real_world/truth/prjna884016/SY22A1_vs_SY221.rdiff` (active, tracked)
- Reader(s): `check_validation_assets.sh:25` (existence only). Zeroed to 0 bytes.
- **ALREADY_GUARDED** (existence gate only; no positive content).

### 4. Stale MUMmer artifacts in the ACTIVE truth dir — `SY22A1_vs_SY221.delta`, `.report`, `.1delta`, `.mdelta`, `.1coords`, `.mcoords`, `.qdiff`, `.snps` (all tracked, Sep-18 content)
- Reader(s): **none in the current codebase** (no script references these filenames), but they are git-tracked in the active namespace with artifact content.
- Can they still yield the invalidated truth? **Yes, in content.** Evidence: `SY22A1_vs_SY221.delta` line 3: `>CP106998.1 CP107005.1 4687476 4492785` (the truncated query length); `.report`: `TotalBases 4687476 4492785` / `UnalignedBases 180982(3.8610%) 0(0.0000%)`. Byte sizes (69,536 and 4,286) match the quarantined legacy copies `SY22A1_vs_SY221_legacy.delta` / `_legacy.report` — they are the same pre-quarantine outputs; only `.rdiff` was zeroed.
- Risk: **Medium.** A future truth loader globbing `truth/prjna884016/*` re-derives the truncation artifact while `check_validation_assets.sh` (which only inspects `.gd`/TSVs/`.rdiff`) keeps passing.
- Required fix: move these to `quarantined_legacy/` with MANIFEST rows (or delete after hashing) and extend the asset gate to assert their absence/emptiness. Not covered by M7-P0a.

### 5. `benchmark/real_world/truth/prjna884016/quarantined_legacy/` (tracked)
- Reader(s): `test_prjna884016_quarantine.py:test_quarantined_legacy_assets_preserved_and_hashes_intact` (hash verification only) and `MANIFEST.tsv`. **No loader treats it as active truth; no glob or discovery path reaches it.** Files are suffixed `_invalid_legacy` and classified `QUARANTINED_INVALID_TRUNCATION_ARTIFACT`.
- Can it be loaded as the active positive set (Q2)? **No.** Risk: low. **ALREADY_GUARDED** (by naming + absence of any consumer; no code-level loader prohibition exists, but nothing loads it today).

### 6. `benchmark/real_world/scripts/build_prjna884016_truth.py` (tracked)
- Symbol: `main()`; guards `validate_assembly_fasta`, `forbid_heuristic_mob_assignment`, `verify_negative_truth`.
- Reads: `data/prjna884016/assemblies/SY221_ref.fna`, `SY22A1.fna` (still on disk; the mutant copy remains truncated at 4,492,785 bp).
- Writes: **nothing** — the generation path no longer exists; `--check-integrity` validates only. Default invocation exits 1: `"ERROR: build_prjna884016_truth.py is QUARANTINED and DISABLED for authoritative validation."` Truncation raises `AssemblyTruncationError` (`"truncated assemblies cannot generate truth"`); IS-length→family inference raises `UnsupportedHeuristicError` (`"strictly prohibited"`).
- **ALREADY_GUARDED.** Risk: none.

### 7. `benchmark/real_world/scripts/compare_mob_truth.py` (tracked)
- Symbols: `compare_mob`, `parse_mob_records`, `main`.
- Zero-denominator behavior (Q3): `precision` → `"NOT_EVALUABLE"` when `tp+fp == 0`; `recall` → `"NOT_EVALUABLE"` when `tp+fn == 0` (comment: `"Zero verified positive events in truth: recall is NOT_EVALUABLE (§6, M7-P0a)"`); `exact_fraction` / `within_5bp_fraction` / `median_error_bp` → `"NOT_EVALUABLE"` when `tp == 0`. With calls-but-no-truth, precision is a real number (`0.0`), recall stays `NOT_EVALUABLE` — verified by `test_empty_truth_scoring_returns_not_evaluable`.
- Precision-vs-truth coupling (Q4): denominator is `tp + fp` (test-side) — **precision is decoupled from the empty truth denominator**; only recall is truth-coupled, and it is guarded.
- **ALREADY_GUARDED.**

### 8. `benchmark/real_world/scripts/compare_variant_calls.py` (tracked)
- Symbols: `parse_gd`, `match_mutations`, `main`.
- Reads: any `--truth-gd` passed on the CLI, including `SY22A1_truth.gd`. Currently **no sbatch/sh invokes it with PRJNA884016 truth** (invocation sites: `run_bl21_validation.sbatch:65`, `run_p1_3_b21_{2_1,3_1,4_1}.sbatch:42` — all BL21 curated truth).
- Can it treat the invalidated events as positive truth? It can no longer read them from the active `.gd` (0 records), but its zero-denominator arithmetic is unsafe (see Unsafe scoring paths).
- Risk: **High** (if anyone scores the SY22A1 negative control through this CLI, it prints a spuriously perfect score). Required fix below.

### 9. `benchmark/real_world/scripts/aggregate_results.py` (tracked)
- Symbol: `aggregate_dir` / `parse_metric_tsv`. Reads every `**/*metrics.tsv` and `summary.tsv` under a results dir; stdout only (no writes).
- Undefined-metric behavior (Q5): missing keys emit `"N/A"` (e.g. `ds_metrics.get("variant_recall", "N/A")`) — it does **not** synthesize numbers. Caveats: `false_complete` and `report_p0_errors` default to `"0"` when absent (silent fail-open defaults), and it would faithfully propagate a fabricated `1.0000` from any stale upstream metrics file (e.g. `results/prjna884016/*`, hard-coded `results/rel606/summary.tsv`).
- Risk: low-medium (pass-through amplifier). Required fix: none mandatory; optionally propagate a `NOT_EVALUABLE` marker and fail loudly on it.

### 10. `benchmark/real_world/VALIDATION_REGISTRY.tsv` line 7 (tracked)
- Content: `prjna884016	SY22A1	ASSEMBLY_CONFIRMED	false	false	false	false	false	PENDING	Assembly confirms 0 chromosomal variants … Historical 4 MOB + 7 JC invalidated … Positive SV sensitivity NOT_EVALUABLE on chromosome; deferred to REL606. Keep all fields false and status PENDING.`
- Q6 answer: **all five booleans are `false`, status `PENDING`, not promoted.** Enforced twice: `check_validation_assets.sh:32-35` (awk on `$4…$9`) and `test_validation_registry_invariants`. Line 6 (`prjna1088182`) untouched and separately gated.
- **ALREADY_GUARDED.**

### 11. `benchmark/real_world/manifests/prjna884016.tsv` (tracked)
- Q1 answer: the manifest **no longer claims positive chromosomal structural truth**. SY22A1 row is `sample_type single_clone` with note `"Negative concordance control: … 0 chromosomal variants). Historical 4 MOB + 7 JC invalidated as local file truncation artifact."` The only 4/7 mentions are repudiations. `truth_source` remains `ASSEMBLY_CONFIRMED` (tier vocabulary, per M7-P0a §D).
- The `assembly_truth` column still points at the active `SY22A1_truth.gd`, but **no code consumes the `assembly_truth` column** (`git grep assembly_truth -- '*.py' '*.sh' '*.sbatch'` → empty). Guarded by `test_manifest_prjna884016_notes`.
- **ALREADY_GUARDED** for scoring; doc-level only.

### 12. `benchmark/real_world/results/prjna884016/{summary.tsv, mob_metrics.tsv, STRUCTURAL_DISCREPANCIES.md}` (LOCAL-ONLY — `results/` is gitignored except `.gitkeep`; dated Sep 18)
- Reader(s)/writer(s): **none** — no current script reads or writes `results/prjna884016` (grep over all scripts: empty). `git log` on the path is empty (never committed).
- Can they still present the invalidated set as positive truth? **Yes, as frozen stale claims.** Evidence: `summary.tsv`: `mob_recall	1.0000`, `jc_recall	1.0000`, `truth_source	ASSEMBLY_CONFIRMED`, notes `Level 3 orthogonal assembly truth derived from complete finished circular genomes (CP107005.1 vs CP106998.1)`; `mob_metrics.tsv`: `tp	4`, `fn	0`, `recall	1.0000`; `STRUCTURAL_DISCREPANCIES.md`: "Mutant Assembly: … `CP107005.1` … **4,492,785 bp, complete circular chromosome**" (the truncated length described as complete) and "8 confirmed breakpoint junctions" (internally inconsistent with the 7 JC elsewhere), with a PASS table for "MOB Recall 100.0%".
- Q7 relevance: these are **not** M6 frozen files and **cannot be rewritten by any current script** — they persist stale. Risk: medium (human misreading; gitignored so invisible to review).
- Required fix: annotate as superseded/repudiated or delete (local artifacts).

### 13. `benchmark/real_world/scripts/check_validation_assets.sh` (tracked)
- The M7-P0a asset gate: asserts registry invariants for lines 6 and 7, 0 MOB/0 JC in the active `.gd`, 0 data rows in both TSVs, tracked-and-not-ignored status of all assets, then runs `m6_reference_policy.py --check` (a PRJNA1088182 execution-reference guard — not a prjna884016 truth consumer).
- **ALREADY_GUARDED** (within its current coverage; does not cover the stale `.delta`/`.report` — see Consumer 4).

### 14. `benchmark/real_world/scripts/test_prjna884016_quarantine.py` (tracked)
- 9 regression tests covering truncation fail-closed, checksum, replicon, heuristic ban, `NOT_EVALUABLE` empty-truth scoring (via `compare_mob` import), 0/0 active truth, quarantine hashes, registry invariants, manifest repudiation.
- **ALREADY_GUARDED.**

### 15. `benchmark/real_world/scripts/download_dataset.sh` (tracked)
- Only mention is an echo: `"Manifest configured in benchmark/real_world/manifests/prjna884016.tsv"` — data staging only, no truth coupling. Safe.

### 16. Documentation consumers
- `benchmark/real_world/M7_P0A_INVALID_TRUTH_REMEDIATION.md`, `M7_P0_DATASET_ELIGIBILITY_AUDIT.md` — current, repudiating (the eligibility audit §3 labels the 4/7 set `FABRICATED_COMPUTATIONAL_ARTIFACT`; §6.3 recommends `truth_level` → `NEGATIVE_CONTROL_ASSEMBLY_CONFIRMED`, which P0a deliberately did not adopt — a documented divergence, registry is self-consistent).
- `benchmark/real_world/README.md:14,19` — **mildly stale framing**: describes `prjna884016.tsv` as "IS/transposon structural ground truth cohort" and `truth/prjna884016/` as "Assembly-derived GD/MOB/JC structural truth" without the negative-control qualification. Doc-only risk; required fix: one-line update.

## Unsafe scoring paths

1. **`compare_variant_calls.py` defaults undefined metrics to 1.0 (100%)** — the central unguarded path (Q3/Q4/Q5):
   - `match_mutations` lines 285-286: `precision = overall_tp / (overall_tp + overall_fp) … else 1.0` and `recall = overall_tp / (overall_tp + overall_fn) … else 1.0`. With the now-empty `SY22A1_truth.gd` and a clean negative-control run (0 calls), it emits `Precision: 100% (0/0)` and `Recall: 100% (0/0)` — a false perfect score for recall that M7-P0a's contract says must be `NOT_EVALUABLE`.
   - Pathological case: empty truth + N false-positive calls gives `precision 0.0` (correct) but `recall = 0/(0+0) → 1.0` — **100% recall alongside N false positives**.
   - Lines 289-291: `breakpoint_exact_fraction` / `within_1bp` / `within_5bp` default to `1.0` when `total_tested_bp == 0` — breakpoint/TSD-style tolerance metrics report 100% with zero truth.
   - TSV mode lines 334-335 repeat the same `else 1.0` per type (`MOB_precision`, `JC_recall`, …).
   - Contrast: `compare_mob_truth.py` returns `"NOT_EVALUABLE"` for the identical situation (Q3 answer: mob scorer → neither 0% nor 100%; variant scorer → **100%**).
   - Q4 answer: **precision is not coupled to the empty positive denominator in either scorer** (test-side `tp+fp` denominators); recall, breakpoint fractions are the truth-coupled metrics, guarded only in the mob scorer.

2. **`summarize_rel606_performance.py` hard-codes perfect accuracy** (adjacent finding, REL606 inventory): lines 88-93 write unconditionally `variant_recall	1.0000`, `mob_recall	1.0000`, `jc_recall	1.0000` etc. into `results/rel606/summary.tsv` — accuracy metrics asserted by construction, not measurement, from a script whose actual scope is wall-clock/RSS parity. `aggregate_results.py` would propagate these as-is. Flag for M8 scoping; no prjna884016 coupling.

3. **Stale artifact content in the active truth namespace** (Consumer 4): tracked `SY22A1_vs_SY221.delta` / `.report` still encode the truncation run; no current reader, but they are the only remaining un-quarantined copies of the invalid evidence.

4. **Frozen local result claims** (Consumer 12): `results/prjna884016/summary.tsv` and `mob_metrics.tsv` (gitignored) still carry `recall 1.0000` measured against the invalidated truth; `STRUCTURAL_DISCREPANCIES.md` still asserts positive structural findings.

## Already guarded

- Active truth assets zeroed and gate-enforced: `SY22A1_truth.gd` (0 MOB/0 JC), `mob_truth.tsv` / `jc_truth.tsv` (header only), `SY22A1_vs_SY221.rdiff` (empty) — `check_validation_assets.sh:39-42`, `test_prjna884016_quarantine.py`.
- `build_prjna884016_truth.py` disabled by default (exit 1), fail-closed on truncation/checksum/replicon/heuristic; writes nothing.
- `compare_mob_truth.py` `NOT_EVALUABLE` semantics on empty truth, with regression test `test_empty_truth_scoring_returns_not_evaluable`.
- Registry line 7 all-`false` + `PENDING` with repudiation notes (Q6: confirmed), enforced by awk gate and pytest; line 6 untouched.
- `quarantined_legacy/` hash-verified, `_invalid_legacy` naming, zero code-path loaders (Q2: cannot be loaded as the active positive set).
- Manifest repudiation notes (Q1: no remaining positive claim in manifest); no code consumes `assembly_truth`.
- Q7: no current tracked script rewrites M6 frozen results (`results/prjna1088182/m6_*`) or `results/prjna884016/*`; `summarize_rel606_performance.py` writes only `results/rel606/*`; `build_prjna884016_truth.py` has no write path at all.

## Required fixes

1. **`compare_variant_calls.py`** (highest priority): replace the `else 1.0` defaults (lines 285-291 and TSV lines 334-335) with a `NOT_EVALUABLE` sentinel (mirroring `compare_mob_truth.py`), or require an explicit `--allow-empty-truth` flag for negative-control scoring; add regression tests symmetric to the mob suite. Until then, prohibit feeding `SY22A1_truth.gd` to this script.
2. **Quarantine the residual active-dir MUMmer artifacts**: move `truth/prjna884016/SY22A1_vs_SY221.{delta,report,1delta,mdelta,1coords,mcoords,qdiff,snps}` into `quarantined_legacy/` with MANIFEST rows (delta/report are byte-identical in size to the archived `_legacy` copies), or hash-and-delete; extend `check_validation_assets.sh` to assert their absence so the gate fails on reintroduction.
3. **Stale local results**: annotate or remove `results/prjna884016/{summary.tsv,mob_metrics.tsv,STRUCTURAL_DISCREPANCIES.md}` (gitignored; unverifiable by review) so the invalidated 100%-PASS table cannot circulate.
4. **`summarize_rel606_performance.py`**: stop emitting hard-coded `1.0000` precision/recall rows; either derive from a real comparison or drop accuracy lines until M8 defines REL606 truth.
5. **Docs**: update `README.md:14,19` to negative-concordance framing; optionally reconcile `M7_P0_DATASET_ELIGIBILITY_AUDIT.md` §6.3 (`NEGATIVE_CONTROL_ASSEMBLY_CONFIRMED` recommendation) with P0a's decision to keep `ASSEMBLY_CONFIRMED`.
6. **Q8 (inventory only, no verdict)**: REL606 is today referenced **only as a breseq comparator**, never as independent assembly truth: registry line 5 (`BRESEQ_ORACLE`, `mob_validated true`, `jc_validated true`, ACTIVE), `manifests/rel606.tsv` (`truth_source BRESEQ_ORACLE`), `summarize_rel606_performance.py` (parity/performance vs breseq 0.40.2 on `SRR030257`), `benchmark/oracle/datasets.tsv` (`REL606_mini` synthetic), `benchmark/PARITY_ISSUES.md` (clonal parity). The "authentic chromosomal IS transpositions … complete hybrid assemblies" language in `M7_P0_DATASET_ELIGIBILITY_AUDIT.md` §5 and M7-P0a §H.2 is a *future M8 pivot recommendation*; no REL606 assembly-truth asset exists under `truth/`. Whether REL606 can serve as independent positive assembly truth is an M8 decision, out of scope here.

## Out of scope dirty files

All worktree-dirty files are M6-P2A (PRJNA1088182) artifacts and were not touched; verified **zero references** to `prjna884016`/`SY22A1` in any of them:

- Modified: `benchmark/real_world/scripts/stage_m6_p1a_fastq.sh`
- Untracked: `benchmark/real_world/M6_P2A_COHORT_VALIDATION.md`; manifests `prjna1088182_m6_p2a_{event_recovery,execution_cohort,fastq_provenance,run_results}.tsv`; scripts `build_m6_p2a_cohort.py`, `m6_p2a_stage_gate.py`, `reconcile_m6_p2a.py`, `summarize_m6_p2a.py`, `test_m6_p2a_cohort.py`, `test_reconcile_m6_p2a.py`, `run_m5_bl21_p1_evidence.sbatch`, `run_m6_p2a_evidence.sbatch`, `stage_m6_p2a_fastq.sbatch`.

No file was created, edited, or deleted during this audit; no Slurm jobs were submitted; `crates/**` was not touched; the untracked M6-P2A set was left as-is.
