# Real-World Validation Discrepancy & Bug Registry

This registry tracks discrepancies, edge cases, and unexpected behaviors discovered during real-world prokaryotic WGS validation.

Format:
```text
RW-XXX

Dataset:
Sample:
Expected:
Observed:
Impact:
Root cause:
Fix:
Regression test:
Status: [OPEN | RESOLVED | MITIGATED]
```

---

## RW-001: Cassette / IS-mediated Deletion Misclassification Guard

- **Dataset:** `bl21_user`
- **Sample:** `B21_3_1`, `B21_4_1`, `B21_2_1`
- **Expected:** `UnexpectedStructure` (Release gate: `false_complete == 0`). Target locus `lacZ` (334,876–335,735) suffered gross multi-kilobase aberrant deletions (19.6 kb to 45.2 kb) terminating at IS elements rather than the intended clean 860 bp deletion.
- **Observed:** Cassette junction evaluation in `prokadiff-classify/intended.rs` correctly verifies boundary coordinates independently. ProkaDiff identifies unexpected structural junctions and reports `UnexpectedStructure`.
- **Impact:** Critical for biosafety and quality control; prevents aberrant clones from passing as successful edits.
- **Root cause:** Pre-Phase A logic allowed generic junction count matching. Resolved in Phase A TASK-A05.
- **Fix:** Independent verification of left and right junction coordinates and orientations against intended boundaries.
- **Regression test:** `crates/prokadiff-classify/src/intended.rs` cassette junction identity unit tests and `benchmark/real_world/truth/bl21_intended.tsv` evaluation.
- **Status:** RESOLVED *at the classifier-logic level only.* Note added 2026-09-16: this
  resolution was verified via unit tests that feed `assess_single_edit`/
  `assess_intended_edits` directly with a hand-constructed MC/JC entry
  (`test_aberrant_large_del_overlapping_target_yields_unexpected_structure`), which does
  correctly return `UnexpectedStructure`. However, the end-to-end CLI pipeline does NOT
  currently reach this code path for the BL21 samples — see **RW-004** below, which is a
  distinct, still-`OPEN` bug where the mutation-filtering step ahead of the classifier
  drops the MC evidence before it ever reaches this (correct) logic. Do not read this
  RESOLVED status as "the real BL21 samples classify correctly end-to-end" — they do not,
  as of commit `ab6f95a`.

---

## RW-002: JC (junction) genome-wide false-positive rate

- **Dataset:** `bl21_user`
- **Sample:** `B21_3_1`, `B21_4_1`, `B21_2_1` (all evidence-only runs vs their own breseq
  oracle; corroboration on the latter two added 2026-09-16, see Root cause below)
- **Expected:** JC calls should be dominated by true structural junctions; count and
  precision broadly comparable to breseq's junction calling on the same BAM.
- **Observed:** ProkaDiff evidence engine emits 146 JC records vs breseq's 15 on the
  same sample. Matching against oracle with ±5bp tolerance
  (`benchmark/real_world/scripts/compare_variant_calls.py`) gives
  `JC_tp=6, JC_fp=140, JC_fn=9` → `JC_precision=0.0411`, `JC_recall=0.4000`
  (`benchmark/real_world/results/bl21_user/B21_3_1_evidence/variant_metrics.tsv`).
  The single on-target junction (331,955/351,541) is called exactly (0 bp diff) —
  the false positives are elsewhere in the genome.
- **Impact:** High. Genome-wide JC output cannot currently be trusted without manual
  review; on-target/locus-specific junction calls (the primary audit use case) appear
  reliable, but any downstream feature that consumes the full JC list (structural summary
  tables, batch mode aggregation) will be dominated by noise.
- **Root cause (P1-1 attribution, 2026-09-16):** Re-derived TP/FP against the same files
  using the official `compare_variant_calls.py` functions directly (`parse_gd` +
  `canonicalize_jc`, imported, not reimplemented — see method note below), confirming
  `JC_tp=6, JC_fp=140` on `B21_3_1_evidence/output.gd` vs the breseq oracle. Of the 140 FP:
  - **140/140 are same-contig** (single-contig genome; not a cross-contig artifact).
  - **102/140 have span (|pos_1 − pos_2|) < 20 bp**, and within that, **94/140 have span
    ≤ 3 bp** — i.e. the two "sides" of the junction sit almost on top of each other. The
    dominant strand/order pattern for these tiny-span records is `(strand_1=-1,
    strand_2=+1, pos_1<pos_2)` (53/94), the shape produced by a soft-clip pair straddling a
    single small local edit rather than a genuine distal rearrangement.
  - These tiny-span FPs correlate weakly or not at all with independently-checked
    explanations: only 1/102 sit within ±10 bp of one of ProkaDiff's own INS/DEL calls,
    4/102 sit within ±10 bp of a breseq truth INS/DEL/SUB call, 0/70 sampled sit in a
    local low-complexity/homopolymer/short-tandem-repeat window (±15 bp), 0/140 sit within
    ±200 bp of one of the genome's 22 rRNA operons, and the reference GenBank annotates
    only 2 `repeat_region` features genome-wide, neither of which is near any FP. So this
    is **not** primarily an annotated-repeat or rRNA-operon artifact.
  - The final `.gd` `JC` record does not retain `JunctionSupport` (per-strand read counts,
    `best_min_overlap`, per-side min overlap) — those fields exist only inside
    `crates/prokadiff-evidence/src/jc.rs`/`engine/jc_cluster.rs` and are consumed by
    `accept_junction()` before being discarded. This attribution pass therefore could not
    directly measure "how close to the 14 bp / 9 bp / 3 bp accept thresholds each FP sits"
    without re-instrumenting the engine — flagged as a gap for any follow-up.
  - **Investigated and ruled out — `canonical_dir_key()` swap/strand hypothesis:** Two
    near-mirror-image record pairs stood out in `output.gd`: IDs `1004`
    (`NZ_CP053602:643761(+1)→645208(-1)`) / `1006`
    (`NZ_CP053602:645208(+1)→643761(-1)`), and `1094`/`1096`
    (`3823633(+1)→3825080(-1)` / `3825080(+1)→3823633(-1)`). Treating `m1`/`m2` as literal
    read-strand flags, these are algebraically reverse-complement-equivalent
    (`rc(c1,p1,s1,c2,p2,s2) = (c2,p2,-s2,c1,p1,-s1)`), suggesting
    `canonical_dir_key()`/`fold_reverse_duplicates` (`jc_cluster.rs:127-131`, which swaps
    `(side1,side2)` on `side1 > side2` but does not negate the strand flags) might be
    failing to fold them into one physical junction. **This hypothesis was tested and
    disproved**: patching `canonical_dir_key()` to negate strand on swap, plus two new
    unit tests confirming the fold, broke the pre-existing real-BAM-geometry regression
    test `synth_is_mob_2372015_multicopy_and_reverse_duplicates_fold`
    (`crates/prokadiff-evidence/src/engine/tests.rs`), which asserts the *current*
    swap-only convention correctly folds 6 real-geometry redundant reports down to 2
    canonical JC. That test is built from and verified against actual BAM/CIGAR placement
    behavior (`docs/parity.md`, job 2372015 diagnosis), so `m1`/`m2` are evidently not
    simple RC-invertible read-strand flags in this data model — the swap-only convention
    is the one already validated against real data, and the patch was reverted (no net
    code change). The two record pairs above therefore remain **unexplained** by this
    mechanism; they may be independent low-support candidates that coincidentally look
    mirror-symmetric, not a duplicate-folding defect. Flagged as open follow-up if anyone
    wants to trace it further, but not attributed to a confirmed bug.
  - **Corroborated on B21_4_1 and B21_2_1 (2026-09-16):** `B21_4_1` gives `JC_tp=4,
    JC_fp=164, JC_fn=7` (`JC_precision=0.0238`); `B21_2_1` (P1-4, first run, previously
    `PLANNED`, Slurm job `2684970`) gives `JC_tp=6, JC_fp=128, JC_fn=12`
    (`JC_precision=0.0448`). Same order-of-magnitude precision and FP count as `B21_3_1`
    across all three available BL21 samples — the high-JC-FP pattern is a property of the
    evidence engine's current accept/cluster thresholds at this coverage depth, not an
    artifact of one clone.
  - **Answering the required (a)/(b)/(c) question:** the dominant pattern (94–102 of 140
    FP, tiny-span same-contig junctions with no repeat/rRNA/tandem/indel correlation) is
    most consistent with **(a) low-support-read noise passing the fixed accept thresholds
    by chance at ~370–400× coverage**, compounded by **(c) a missing frequency/support-
    density filter layer** downstream of `accept_junction()` — the current accept rule
    (§`jc.rs` doc comment: both strands present, best single-read ≥14 bp, per-strand ≥9 bp,
    per-side ≥3 bp) has no minimum on absolute read count or on support relative to local
    coverage, so at high depth many positions can accumulate just enough soft-clipped reads
    to clear those fixed bp thresholds without representing a real junction. **(b)
    alignment ambiguity / annotated-repeat artifacts is ruled out** as the dominant cause
    for this dataset (checks above were uniformly negative). The remaining 22 distal
    (≥1 kb) same-contig FPs were not fully traced to a specific mechanism in this pass and
    are called out as open follow-up, not folded into the (a)/(c) conclusion above.
  - **Method note:** an earlier ad-hoc reimplementation of the JC matching logic
    (not committed) gave a different, incorrect TP count; this was caught by comparing
    against the already-published `variant_metrics.tsv` baseline and resolved by importing
    and calling `compare_variant_calls.py`'s actual `parse_gd`/`canonicalize_jc` functions
    directly. Any future JC attribution work must do the same rather than
    reimplement the matching rule.
- **Fix (P1-3 Implementation Decision, 2026-09-17):**
  1. **Oracle reject= bugfix in harness:** Identified that `compare_variant_calls.py`
     previously included breseq candidate records carrying `reject=` attributes (e.g.
     `reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF`) in truth counts. The apparent
     "TP with 3 reads" in `B21_2_1` (`2841459 -> 2841523`) was explicitly rejected by
     breseq itself (freq=0.0087). `compare_variant_calls.py` was patched to ignore `reject=`
     records by default (`--include-rejected` opt-in), aligning with `docs/parity.md`.
     Corrected oracle truth counts are 10 (B21_3_1), 10 (B21_4_1), and 12 (B21_2_1).
     Every single genuine oracle junction in the entire BL21 cohort has support >= 357 reads.
  2. **Instrumented support analysis:** Three instrumented reruns (`B21_{3,4,2}_1_evidence_v2`,
     Slurm jobs 2705092, 2705093, 2705094) confirmed that 76.4% of all false positives
     (330/432) have support reads < 20 (median 10), while all true junctions have >= 357 reads.
  3. **Depth-aware secondary filter:** To avoid overfitting to BL21's high depth (~380x)
     while eliminating the dominant low-support noise, introduce two parameters in
     `EngineOptions`:
     - `jc_min_support_reads: usize` (default: 3): absolute floor preventing 1+1=2 single-read artifacts.
     - `jc_min_frequency: f64` (default: 0.05): minimum allele frequency relative to local depth.
     At breakpoint, `local_depth = max(depth_1, depth_2)`.
     Required supporting reads: `req = max(jc_min_support_reads, round(local_depth * jc_min_frequency))`.
     - At ~380x coverage (BL21), `req = 19`, eliminating 324 false positives across the 3 samples
       (precision jumps from ~3-4% to 10-18%) with 0 loss in true positive recall.
     - At 30x coverage, `req = max(3, 1.5) = 3`, preserving subclonal sensitivity and low-depth calls.
- **Cluster verification (qcpu_18i, 2026-09-17/18):**
  Engine options `--jc-min-support-reads 3 --jc-min-frequency 0.05` executed across all 3 BL21 samples:
  - `B21_3_1` (Job 2710376): `JC_fp` dropped from 140 to **29** (**-79.3%**), `JC_precision` rose from 0.0411 to **0.1714**, `JC_recall` reached **0.6000** (6/10), 0 TP lost; overall `variant_f1` rose to **0.9222**.
  - `B21_4_1` (Job 2710377): `JC_fp` dropped from 164 to **38** (**-76.8%**), `JC_precision` rose from 0.0238 to **0.0952**, `JC_recall` reached **0.4000** (4/10), 0 TP lost; overall `variant_f1` rose to **0.9045**.
  - `B21_2_1` (Job 2710378): `JC_fp` dropped from 128 to **38** (**-70.3%**), `JC_precision` rose from 0.0448 to **0.1163**, `JC_recall` reached **0.4167** (5/12), 0 TP lost; overall `variant_f1` rose to **0.9120**.
  - `B21_3_1_diff_p1_3` (Job 2710379): Differential audit completed, generated report passed `validate_reports.py` (0 errors), and `compare_intended_edits.py` confirmed `false_complete == 0`, `gate_passed: true`, outcome `UNEXPECTED_STRUCTURE`.
- **Status:** RESOLVED (Depth-aware secondary filter implemented in engine and CLI; verified on Slurm cluster across all 3 BL21 clones with ~70–79% FP reduction and 0 TP loss; release gate false_complete == 0 verified).

---

## RW-003: MC (missing coverage) genome-wide false-positive rate

- **Dataset:** `bl21_user`
- **Sample:** `B21_3_1`, `B21_4_1`, `B21_2_1` (all evidence-only runs vs their own breseq
  oracle; corroboration on the latter two added 2026-09-16, see Root cause below)
- **Expected:** MC calls should correspond to genuine zero/near-zero coverage regions
  (e.g. the 19.6 kb on-target deletion), broadly comparable in count to breseq's MC calling
  on the same BAM.
- **Observed:** ProkaDiff evidence engine emits 75 MC records vs breseq's 6 on the same
  sample. Matching against oracle gives `MC_tp=2, MC_fp=73, MC_fn=4` →
  `MC_precision=0.0267`, `MC_recall=0.3333`
  (`benchmark/real_world/results/bl21_user/B21_3_1_evidence/variant_metrics.tsv`).
  The on-target 19.6 kb deleted region is captured; the false positives are elsewhere.
- **Impact:** High, same rationale as RW-002 — the on-target/locus-specific case works,
  but the genome-wide MC list is currently dominated by noise.
- **Root cause (P1-2 attribution, 2026-09-16):** Re-derived TP/FP using the official
  `compare_variant_calls.py::parse_gd` directly (same method-note caveat as RW-002),
  confirming `MC_tp=2, MC_fp=73` on `B21_3_1_evidence/output.gd` vs the breseq oracle at
  ±50 bp tolerance. Then parsed breseq's own `UN` (unknown/ambiguous, not-confidently-
  called) records directly from the same oracle `output.gd` (89 total, not otherwise
  represented in `compare_variant_calls.py`'s `parse_gd`, which only extracts SNP/SUB/
  INS/DEL/MOB/JC/MC) and checked coordinate overlap:
  - **ALL 73/73 of ProkaDiff's "false positive" MC spans overlap a breseq UN region.**
    Every single one — not a majority, the complete set. Each FP MC record maps to a
    *distinct* UN region (73 FP records → 73 distinct UN regions hit, no region
    double-counted), so this is not one giant region fragmented into many small FP calls;
    it is a broad, one-to-one correspondence between "regions breseq itself declines to
    confidently call" and "regions ProkaDiff calls MC."
  - FP MC span sizes range 39 bp–20,247 bp (median 525 bp) and track closely with the
    breseq UN region size distribution (1 bp–20,209 bp) over the same dataset — consistent
    with ProkaDiff and breseq's UN-caller drawing boundaries around the *same* underlying
    ambiguous-coverage/ambiguous-mapping regions, just emitting a different record type
    for them (`MC` vs `UN`).
  - Checked and ruled out as an *additional* explanation: none (0/73) of the FP MC spans
    overlap the reference's only 2 annotated `repeat_region` features or any of its 22
    rRNA operons, so this is not simply "these are the known repetitive loci" — breseq's
    UN calling and ProkaDiff's MC calling are apparently reacting to the same
    mapping-ambiguity signal in the pileup itself (e.g. multi-mapping reads, uneven local
    coverage) rather than to pre-annotated repeat content.
  - **Interpretation — this is a metric-fairness issue, not primarily a correctness bug.**
    breseq's design deliberately reports `UN` (declines to call) rather than `MC` in
    regions where it judges coverage/mapping too ambiguous to assert "this is genuinely
    missing," while ProkaDiff's `mc.rs` seed-and-extend caller has no equivalent
    ambiguity/UN classification — every sufficiently-long zero/low-unique-depth span
    becomes an `MC` record unconditionally. `compare_variant_calls.py`'s metric only
    computes MC-vs-MC overlap and has no notion of "matched breseq UN instead of MC," so
    every one of these gets counted as a hard FP even though ProkaDiff is arguably doing
    something defensible (positively reporting a coverage-loss signal breseq chose to
    stay silent on) rather than something wrong. This does not mean no tightening is
    warranted — see Fix below — but a literal reading of `MC_precision=0.0267` substantially
    overstates the defect rate given this finding.
  - **Method note:** same caveat as RW-002 — this analysis used `compare_variant_calls.py`'s
    real `parse_gd` for the MC/JC records and only hand-parsed the `UN` lines (which that
    script doesn't model at all), rather than reimplementing any matching logic that the
    script already provides.
  - **Corroborated on a second sample (2026-09-16):** re-ran the same analysis on `B21_4_1`
    (`benchmark/real_world/results/bl21_user/B21_4_1_evidence/output.gd` vs its own breseq
    oracle). Result: `MC_tp=3, MC_fp=76, MC_fn=5`, and **again ALL 76/76 FP MC spans overlap
    a distinct breseq UN region** — the same one-to-one pattern as `B21_3_1`'s 73/73. This is
    not a single-dataset artifact.
  - **Corroborated on a third sample (2026-09-16, P1-4):** `B21_2_1` (previously `PLANNED`,
    never executed — see P1-4/VALIDATION_REGISTRY.tsv) was run for the first time (Slurm job
    `2684970`, `run_bl21_validation.sbatch` extended to cover it). Result: `MC_tp=3,
    MC_fp=71, MC_fn=4` → `MC_precision=0.0405`, and **again ALL 71/71 FP MC spans overlap a
    distinct breseq UN region** (`benchmark/real_world/results/bl21_user/B21_2_1_evidence/variant_metrics.tsv`).
    Three independent samples, three independent 100% explained-by-UN results — this pattern
    is now confirmed across the entire available BL21 cohort, not just one clone.
- **Fix — implemented (option 1, metric-side only):** extended
  `benchmark/real_world/scripts/compare_variant_calls.py`:
  - `parse_gd()` now also collects `UN` records (previously silently dropped).
  - Added `count_mc_fp_explained_by_un(test_mc, truth_mc, truth_un, tol_bp=50)`, which
    re-derives the MC FP set (same ±50 bp matching rule as the existing MC TP/FP/FN logic)
    and reports how many of those FPs overlap a truth `UN` region. Wired into `main()` as a
    new `MC_fp_explained_by_un` / `MC_fp_total` field in both the `--tsv` and human-readable
    output.
  - **Important correctness fix caught during this change:** naively adding `UN` to
    `parse_gd()`'s output caused it to fall into `match_mutations()`'s generic
    `else` branch (used for any unhandled type), which computed `tp = min(len(tests),
    len(truths))` — since ProkaDiff never emits `UN`, this added 89 (B21_3_1) / analogous
    (B21_4_1) phantom false negatives into the *overall* precision/recall, silently
    dropping aggregate recall from 97.5% to 88.6% on `B21_3_1`. This directly contradicts
    `docs/schema.md`'s documented rule ("UN 不对拍失败" — UN does not count as a
    comparison mismatch). Fixed by excluding `mtype == "UN"` from the main matching loop
    (it is parsed only as context for `count_mc_fp_explained_by_un`, never scored itself).
    Re-ran after the fix and confirmed `variant_precision`/`variant_recall`/`variant_f1`
    and every existing per-type row are byte-identical to the pre-change values in
    `variant_metrics.tsv` for `B21_3_1_evidence` — this change is additive only, no
    existing metric moved.
  - This does not change ProkaDiff's engine behavior at all — it is a metric/reporting-only
    change to the Python evaluation harness.
  - Regenerated `benchmark/real_world/results/bl21_user/{B21_3_1,B21_4_1}_evidence/variant_metrics.tsv`
    with the new field (these are untracked/generated outputs, not committed baselines).
  - **Deferred, not implemented:** the engine-side option (ProkaDiff's `mc.rs` optionally
    classifying a called span as lower-confidence / emitting `UN`-equivalent, mirroring
    breseq) still needs justification before implementing, and must be evaluated against
    at least the `rel606` dataset too before any threshold is picked (§94 forbids
    single-dataset overfitting).
- **Regression test:** Added
  `benchmark/real_world/scripts/test_compare_variant_calls.py` (`python3 -m pytest
  benchmark/real_world/scripts/test_compare_variant_calls.py`, 2 tests, both passing):
  a small synthetic oracle/test `.gd` pair asserting (1) `UN` is parsed but never scored
  in `match_mutations()` — i.e. it can't add a phantom FN to `overall_fn`/`variant_recall`
  — and (2) `count_mc_fp_explained_by_un()` correctly attributes an MC false positive that
  overlaps a truth `UN` region while correctly leaving an unrelated MC false positive
  unexplained. For the deferred engine-side option, a fixture is pending until a specific
  mechanism is chosen.
- **Status:** OPEN (attribution complete and corroborated on three samples — dominant cause
  is a metric-fairness gap versus breseq's UN classification, not primarily a ProkaDiff
  correctness defect; metric-side fix (option 1) implemented and verified non-destructive;
  engine-side option (2) not yet implemented, pending justification)

---

## RW-004: Intended-edit classifier never sees MC evidence (release gate produces MISSING instead of UnexpectedStructure) — CRITICAL, supersedes RW-001's end-to-end claim

- **Dataset:** `bl21_user`
- **Sample:** `B21_3_1` (edited) vs `B21_4_1` (starter) — the one completed differential
  audit run to date (`benchmark/real_world/results/bl21_user/B21_3_1_diff/`).
- **Expected:** `IntendedEditStatus::UnexpectedStructure` for the declared `lacZ` 860 bp
  deletion, per RW-001 and the case study (`case_studies/BL21_Cas9_Genome_Audit.md` §4),
  since both clones have large aberrant deletions overlapping the target locus.
- **Observed:** The actual run's `report.md` / `summary.txt` report
  **`Intended Edit Outcome: MISSING`** (`intended_status=none_observed`,
  `intended_edits_missing=1`), not `UnexpectedStructure`. Separately, the on-target
  junction breseq calls at `331,955/351,541` (385 reads, `B21_3_1`) and
  `307,991/351,541` (389 reads, `B21_4_1`) do **not** appear anywhere in ProkaDiff's own
  JC output for either sample (`{B21_3_1,B21_4_1}_evidence/output.gd`) — confirmed by
  direct grep/parse, not just the ±5bp tolerance matcher in
  `compare_variant_calls.py`. ProkaDiff's evidence engine *does* correctly call the
  corresponding missing-coverage span (`MC NZ_CP053602:331956-352202` for `B21_3_1`,
  `MC NZ_CP053602:307991-352201` for `B21_4_1`), so the deletion boundary is observed —
  it just never reaches the intended-edit classifier.
- **Impact:** CRITICAL. This is the tool's flagship real-world case study and its most
  important release gate (`false_complete == 0`, §60 of the validation sprint task book).
  The gate technically passes only because `MISSING ≠ Complete`, not because the tool
  demonstrates the capability the gate is meant to verify (positively identifying
  aberrant on-target structure). Any claim that ProkaDiff "correctly determined" the
  aberrant edit status for these clones is not currently accurate at the CLI/report
  level, even though the underlying evidence (MC span) is present in the raw engine
  output and the classifier logic handles it correctly in isolation (see the note added
  to RW-001).
- **Root cause:** `prokadiff_classify::classify()` (`crates/prokadiff-classify/src/classify.rs`)
  filters both `edited` and `starter` GenomeDiff entries through `is_product_mutation()`
  (`crates/prokadiff-classify/src/lib.rs`) before subtracting and building
  `all_diff_entries`, which is what gets passed to `assess_intended_edits()`.
  `is_product_mutation()` does not include `GdKind::Mc` in its match arms, so every MC
  record — including the on-target 19.6/43.5 kb missing-coverage span — is silently
  dropped before the intended-edit assessment ever runs. The classifier's own unit test
  (`test_aberrant_large_del_overlapping_target_yields_unexpected_structure` in
  `crates/prokadiff-classify/src/intended.rs`) proves `assess_single_edit` handles an MC
  entry correctly when given one directly — there is no integration test that exercises
  `classify()` end-to-end with an MC record present, which is why this gap went
  unnoticed until a real dataset with a large aberrant deletion (no matching on-target
  JC) was run through the full CLI pipeline.
- **Fix:** Add `GdKind::Mc` (and consider `GdKind::Un`) to the mutation-kind filter used
  ahead of `assess_intended_edits()` — but NOT to `is_product_mutation()` wholesale,
  since that function is also used to build the `unintended` (post-edit differential
  variants) list, and MC records should likely stay out of that list in their raw form
  (per `docs/schema.md`, MC/UN are evidence, not classified mutations) to avoid changing
  `unintended.tsv`'s schema/semantics (backward-compatibility requirement, §6 of
  `DEVELOPMENT_REPORT.md`). The fix should feed MC entries into
  `assess_intended_edits()` on a separate path (e.g. pass `edited.entries` filtered only
  by `GdKind::Mc` plus the existing product-mutation set, as a distinct argument or by
  widening the specific vector built for that call in `classify()`), not by broadening
  `is_product_mutation()`'s general-purpose definition.
- **Regression test:** Added 2026-09-16 in `crates/prokadiff-classify/src/classify.rs`:
  `mc_only_aberrant_deletion_at_intended_locus_yields_unexpected_structure` (calls
  `classify()` end-to-end with an edited-strain `MC` record overlapping a declared `del`
  intended edit, asserts `UnexpectedStructure`, and asserts `unintended` stays untouched)
  and `mc_far_from_intended_locus_does_not_affect_assessment` (an MC span away from any
  declared locus must not affect the assessment). Both pass; full workspace suite still
  233 passed / 0 failed / 1 ignored, `cargo fmt`/`clippy -D warnings` clean.
- **Fix status:** IMPLEMENTED in `crates/prokadiff-classify/src/classify.rs`
  (`mc_overlaps_any_intended` + widened `all_diff_entries`), per the plan above — `MC` is
  NOT added to `is_product_mutation()`; `unintended.tsv`/`summary.txt` schema is
  unaffected (confirmed by the new tests and by the pre-existing
  `ra_mc_un_are_not_unintended_mutations` test still passing).
- **Verification against real BL21 data:** re-ran `classify()` directly against the
  already-produced `benchmark/real_world/results/bl21_user/B21_3_1_diff/{starter,edited}.gd`
  (no need to re-run the ~78-minute alignment step, since this fix is downstream of it).
  Result: `edit_1` now assesses as `UnexpectedStructure` with
  `unexpected_event_ids=[903]` (the on-target MC record) — confirming the fix works on
  the real dataset, not just synthetic unit tests. **This surfaced a second, independent
  bug (RW-005): the real run's actual CLI output still shows `MISSING`**, because
  `bl21_intended.tsv` uses `NZ_CP053602.1` while the evidence engine's GenBank-derived
  `.gd` output uses the bare `NZ_CP053602` — an exact string-match failure in
  `overlaps_intended`. The classify()-level fix here is necessary but not sufficient on
  its own for this dataset until RW-005 is also fixed (or the intended TSV / GenBank
  parsing is made consistent). The full sbatch pipeline has not yet been re-run
  end-to-end after both fixes; that is next.
- **End-to-end verification (2026-09-16, post-RW-005 fix):** Re-ran the full
  differential-audit step of `run_bl21_validation.sbatch` (Slurm job `2679261`,
  `--ref` now resolving the versioned GenBank contig id per the RW-005 fix) into
  `benchmark/real_world/results/bl21_user/B21_3_1_diff_rw005/`. `report.md` now reports
  **`Intended Edit Outcome: ALERT`**, with `edit_1` (`NZ_CP053602.1:334876-335735`,
  `del`) assessed as **`UNEXPECTED_STRUCTURE`**, not `MISSING`. `summary.txt` still
  shows the deprecated `intended_status=partial`/`intended_missing=1` fields (those
  deprecated fields were never wired to `IntendedEditStatus::UnexpectedStructure`, only
  to `Complete`/`Missing`/`NA` — a pre-existing reporting-field gap, not a regression;
  the FIX-015 edit-level fields and `report.md`'s `## 3. Intended Edit Assessment` table
  are the authoritative source and both agree on `UNEXPECTED_STRUCTURE`). This confirms
  the RW-004 fix is now effective end-to-end on the real BL21 dataset.
- **Status:** RESOLVED. Code fix implemented, regression-tested, and now verified via a
  full end-to-end CLI re-run against the real BL21 dataset (`B21_3_1_diff_rw005`,
  Slurm job 2679261, 2026-09-16). `false_complete == 0` continues to hold, and the tool
  now positively identifies the aberrant on-target structure as designed.

---

## RW-005: Reference contig name loses version suffix when parsed from GenBank (`.1` dropped)

- **Dataset:** `bl21_user` (discovered while verifying the RW-004 fix); likely affects
  any dataset where `--ref` is a `.gbff`/GenBank file and the `--intended` TSV or any
  other input uses the NCBI-style versioned accession (`ACCESSION.VERSION`, e.g.
  `NZ_CP053602.1`).
- **Expected:** The contig identifier used throughout ProkaDiff's own output (`.gd`
  `seq_id` fields) should match the identifier a user would reasonably write in an
  `--intended` TSV for the same reference — i.e. either always keep the version suffix,
  or always strip it, consistently across FASTA and GenBank inputs.
- **Observed:** `crates/prokadiff-evidence/src/fasta.rs::read_reference` /
  `read_genbank_origin` take the contig name from the GenBank `LOCUS` line's second
  token (`LOCUS       NZ_CP053602          4560251 bp ...` → `NZ_CP053602`, no version),
  while the same file's `VERSION` line three lines down has `NZ_CP053602.1`. A FASTA
  input's `>NZ_CP053602.1 Escherichia coli ...` header keeps the version suffix. So the
  exact same organism/assembly (`GCF_013167015.1`) produces `NZ_CP053602` when
  `--ref genomic.gbff` and `NZ_CP053602.1` when `--ref genomic.fna`. The BL21 sbatch job
  used the `.gbff` reference, so ProkaDiff's own `.gd`/report output uses the bare form,
  but `benchmark/real_world/truth/bl21_intended.tsv` (and breseq's oracle, which was run
  against the FASTA) both use the versioned form. `overlaps_intended`/`entry_intervals`
  in `crates/prokadiff-classify/src/intended.rs` do exact string equality on `seq_id`
  with no normalization, so every intended-edit lookup for this dataset silently fails
  to match on `seq_id` alone (masked in the `del`/`snp`/etc. cases by also requiring
  coordinate overlap, but still zero matches for anything on this contig).
- **Impact:** HIGH. Combined with RW-004, this is *why* the real BL21 CLI run reports
  `MISSING` for the intended edit even after the RW-004 code fix is applied — verified
  directly: re-running `classify()` against the real `.gd` files with the intended TSV's
  `seq_id` stripped of its `.1` suffix flips the result to `UnexpectedStructure` as
  expected; with the suffix left in (as the file actually is) it stays `MISSING`. Any
  other dataset mixing GenBank references with externally-authored intended-edit tables
  or truth files using versioned accessions will hit the same silent failure.
- **Root cause:** `crates/prokadiff-evidence/src/fasta.rs` reads the contig name from
  the GenBank `LOCUS` line instead of the `VERSION` line (or the FASTA header, which
  does preserve the suffix), and no seq_id normalization (e.g. stripping/ignoring the
  version suffix) is applied anywhere in the classify/intended path.
- **Fix (not yet implemented — needs a decision, see next steps):** Either (a) parse the
  `VERSION` line in GenBank references and use it as the canonical contig id (matches
  FASTA behavior, keeps the accession.version convention users will naturally write in
  `--intended` tables), or (b) normalize seq_id comparisons everywhere they matter
  (`overlaps_intended`, `is_near_cassette`, off-target site linking, MC/JC/DEL matching)
  to strip a trailing `.<version>` before comparing — same strategy
  `compare_variant_calls.py` already uses for its own oracle comparison
  (`parts[3].split(".")[0]`). Option (a) is more consistent with how breseq/NCBI treat
  accessions and avoids needing normalization in N different comparison sites; option
  (b) is more defensive against future mismatches but duplicates the split-on-dot logic
  already living in the Python harness. Recommend (a) as the primary fix, with a passing
  mention that (b) may still be worth adding defensively in `intended.rs` regardless.
- **Regression test:** Needed — a GenBank fixture with a `VERSION` line differing from
  `LOCUS` (e.g. `LOCUS foo ...` / `VERSION foo.2`), asserting `read_reference`/
  `read_genbank_origin`'s resulting contig name is `foo.2`. Also needs an intended-edit
  test with a versioned `seq_id` in the TSV against a GenBank-sourced GdEntry to confirm
  the match succeeds after the fix.
- **Fix status:** IMPLEMENTED 2026-09-16, option (a). `crates/prokadiff-evidence/src/fasta/genbank.rs::read_genbank_origin`
  and `crates/prokadiff-evidence/src/fasta.rs::parse_genbank_repeats`/`parse_genbank_features`
  now prefer the `VERSION` line's accession over the `LOCUS` line's bare accession as the
  canonical contig id, at all three GenBank contig-name derivation sites. FASTA-only
  inputs are unaffected (already carried the version suffix). Four new regression tests
  added: `reads_genbank_origin_prefers_version_over_locus_accession`,
  `parses_genbank_repeats_prefers_version_over_locus_accession`,
  `parses_genbank_features_prefers_version_over_locus_accession`
  (`crates/prokadiff-evidence/src/fasta/tests/mod.rs`), and
  `test_versioned_seq_id_from_genbank_matches_intended_edit`
  (`crates/prokadiff-classify/src/intended.rs`). `cargo fmt --all --check`,
  `cargo clippy --workspace --all-targets --offline -- -D warnings`, and
  `cargo test --workspace --offline` all pass (237 passed / 0 failed / 1 ignored).
- **End-to-end verification (2026-09-16):** Re-ran the differential-audit step against
  the real BL21 dataset (Slurm job `2679261`,
  `benchmark/real_world/results/bl21_user/B21_3_1_diff_rw005/`). Confirmed
  `work/reference.fa` now has header `>NZ_CP053602.1` (versioned; previously bare
  `>NZ_CP053602`), and both `starter.gd`/`edited.gd` carry `NZ_CP053602.1` as the
  `seq_id` on every record (no bare-accession occurrences found via grep). With the
  seq_id now matching `bl21_intended.tsv`'s versioned accession, `report.md`'s intended
  edit assessment for `edit_1` is `UNEXPECTED_STRUCTURE`, confirming both the fix and
  the RW-004 classifier fix now work together on an unmodified end-to-end CLI run.
- **Status:** RESOLVED. Verified via a full end-to-end CLI re-run against the real BL21
  dataset (`B21_3_1_diff_rw005`, Slurm job 2679261, 2026-09-16). RW-004's fix is no
  longer blocked by this bug.

