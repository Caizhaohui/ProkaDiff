# M7-P0b Scientific Review

**MODEL:** Claude Opus 5.5 at high effort

---

## Decision A: Does PRJNA884016 have any valid positive chromosomal structural truth?

**Verdict:** NO

**Evidence:**
- NCBI assemblies CP107005.1 (SY22A1) and CP106998.1 (SY221) are 100.000% bit-identical across all 4,687,476 bp (M7_P0A §B checksum table, all strains share SHA256 `524449dac020a114fbb372aa453bc6d6a3fb962dd2489407ea2107cc115e7e1a`)
- Active truth file `SY22A1_truth.gd` contains zero MOB and zero JC records (lines 1-4: header only, explicit note "Historical 4 MOB + 7 JC invalidated as 194.7 kb truncation artifact")
- `mob_truth.tsv` and `jc_truth.tsv` contain headers only, 0 data rows (M7_P0B §2 consumer inventory)
- All 4 MOB + 7 JC events in legacy truth resided strictly within coordinates 4,492,786–4,687,476, the truncated tail region (M7_P0A §A.3)

**Limitation:** Genuine IS transpositions occurred exclusively on the extrachromosomal cas-expression plasmid, which was never deposited in NCBI (M7_P0A §A.2; M7_P0 §2 replicon audit shows 1 replicon per assembly, chromosome only).

---


Provenance: Claude Code 2.1.216, model `claude-opus-5-5`, effort `high`, permission mode `dontAsk`, tools restricted to Read/Grep/Glob/Bash, session `80f44b04-2177-4101-850e-58e5168d0a40`, exit 0. An earlier `--permission-mode plan` invocation produced no output for about 40 minutes and was stopped; it is not this review. The reviewer did not modify the worktree; the orchestrator saved this document.



## Decision B: Do chromosome assemblies support descriptive zero-difference concordance?

**Verdict:** YES

**Evidence:**
- Complete circular finished assemblies for both strains (CP107005.1 and CP106998.1, 4,687,476 bp each)
- Bit-identical sequences confirmed by normalized SHA256 match (M7_P0A §B.1)
- Assembly status `Complete Genome, Chromosome` for both accessions (M7_P0 §2)
- VALIDATION_REGISTRY.tsv line 7 notes: "Assembly confirms 0 chromosomal variants between CP107005.1 and CP106998.1 (100.000% bit-identical)"
- This concordance aligns with Sheng et al. 2023 primary finding: "100% sequence identity at 100% chromosome coverage"

**Limitation:** Zero-difference concordance is descriptive assembly fact, not a measure of tool specificity until a differential caller is executed against matched WGS pairs.

---

## Decision C: Which metrics are mathematically defined when the positive truth set is empty?

| Metric | Denominator | Defined? | Evidence |
|:---|:---|:---:|:---|
| **Recall** | `tp + fn` (truth-side) | **NOT_EVALUABLE** | 0 truth → `fn = 0` → denominator 0; compare_mob_truth.py:76-80 returns `"NOT_EVALUABLE"` |
| **Sensitivity** | `tp + fn` | **NOT_EVALUABLE** | Synonym for recall; same denominator |
| **Breakpoint accuracy** (exact/5bp fractions) | `tp` | **NOT_EVALUABLE** | compare_mob_truth.py:82-89 returns `"NOT_EVALUABLE"` when `tp == 0` |
| **TSD accuracy** | `tp` (matched MOB events) | **NOT_EVALUABLE** | Cannot measure TSD concordance with 0 matched events |
| **Precision** | `tp + fp` (test-side) | **INDEPENDENTLY DEFINED** | Denominator is the number of tool calls; defined as long as tool was executed (can be 0.0 if fp > 0, or NOT_EVALUABLE if 0 calls made) |
| **Specificity** | `tn + fp` | **INDEPENDENTLY DEFINED** | Requires independent classification of true negatives via assembly agreement or negative-site audit; not automatically derived from empty truth |
| **False-positive rate** | `fp / (fp + tn)` | **INDEPENDENTLY DEFINED** | Requires independent enumeration of true-negative site space |

**Limitation:** Precision is mathematically defined but biologically uninterpretable when there are zero true events to contextualize the false-positive burden. Specificity and FP rate require independent truth classification of the zero-call expectation at negative loci, not merely the absence of positive truth.

---

## Decision D: When may an unmatched ProkaDiff call be classified as a biological false positive?

**Verdict:** When assembly-confirmed zero-difference truth independently supports the expectation of zero variants at the called locus.

**Evidence:**
- PRJNA884016 chromosomes are bit-identical complete assemblies → any structural call on this pair contradicts assembly truth → biological false positive
- This is distinct from "unmatched by breseq oracle" (algorithmic disagreement) or "unmatched by partial truth" (incomplete gold standard)
- Classification requires independent assembly evidence or validated negative control status, not merely the absence of a matching positive record

**Limitation:** Cannot classify unmatched calls as biological false positives when:
- Truth is incomplete (e.g., PRJNA1088182 publication-reported truth covering only intended edits, not all biological variation)
- Assembly truth is unavailable or low-quality
- Comparison is oracle-vs-oracle (breseq vs ProkaDiff) without independent assembly ground truth

---

## Decision E: Is plasmid IS transposition evidence appropriately separated from chromosomal truth?

**Verdict:** YES

**Evidence:**
- Active truth files (`SY22A1_truth.gd`, `mob_truth.tsv`, `jc_truth.tsv`) contain 0 events (chromosomal truth only)
- `manifests/prjna884016.tsv` line 2 notes: "Positive IS transpositions occurred into plasmid-borne cas genes, which were not deposited in NCBI assemblies"
- VALIDATION_REGISTRY.tsv line 7 notes: "Positive SV sensitivity NOT_EVALUABLE on chromosome; deferred to REL606"
- M7_P0 §3 eligibility matrix classifies plasmid MOB as `NOT_EVALUABLE` due to "No plasmid assembly"
- No plasmid sequence deposited under PRJNA884016 (M7_P0A §A.2; M7_P0 §2 all 12 assemblies: 1 replicon = chromosome only)

**Limitation:** Genuine biological IS transpositions documented by Sheng et al. 2023 cannot be evaluated against WGS data because the plasmid replicon was never assembled or deposited.

---

## Decision F: Does VALIDATION_REGISTRY.tsv row 7 preserve defensible truth-level semantics?

**Verdict:** YES, with documented architectural tension

**Evidence:**
- Line 7: `prjna884016 SY22A1 ASSEMBLY_CONFIRMED false false false false false PENDING`
- `truth_level = ASSEMBLY_CONFIRMED` reflects the methodological tier (complete assembly evidence), not polarity (M7_P0A §D.1)
- All five validation booleans are `false` (variant_validated, intended_validated, mob_validated, jc_validated, report_reviewed)
- Status is `PENDING` (not promoted to ACTIVE)
- Notes explicitly repudiate the truncation artifact and state "Positive SV sensitivity NOT_EVALUABLE on chromosome"
- Enforced by `check_validation_assets.sh:32-35` and `test_validation_registry_invariants` (M7_P0B §10)

**Limitation:**
- M7_P0 §6.3 recommended `truth_level → NEGATIVE_CONTROL_ASSEMBLY_CONFIRMED`, but M7-P0a retained `ASSEMBLY_CONFIRMED` under Claude Opus 5.5 guidance to avoid enum proliferation
- The string `ASSEMBLY_CONFIRMED` without context could imply positive structural truth; defensibility depends on strict co-reading with the all-false validation booleans and PENDING status

---

## Decision G: Does REL606 provide independent assembly truth, or only breseq comparator evidence?

**Verdict:** ONLY breseq comparator evidence at present

**Evidence:**
- VALIDATION_REGISTRY.tsv line 5: `rel606 SRR030257 BRESEQ_ORACLE true false true true true ACTIVE`
- `truth_level = BRESEQ_ORACLE` (not `ASSEMBLY_CONFIRMED`)
- `manifests/rel606.tsv`: `truth_source BRESEQ_ORACLE` (M7_P0B §16, Q8)
- No assembly-truth asset exists under `truth/rel606/` (M7_P0B: "no REL606 assembly-truth asset exists under `truth/`")
- `summarize_rel606_performance.py` compares ProkaDiff vs breseq 0.40.2 output parity and runtime, not against independent assembly truth
- M7_P0 §5 and M7_P0A §H.2 recommend REL606 LTEE evolved isolates as a **future M8 pivot** for positive structural calibration

**Limitation:** Current validation uses breseq as the oracle; breseq-vs-breseq is output parity, not accuracy. Assembly-confirmed positive structural truth for REL606 chromosomal IS insertions and deletions is deferred to Milestone M8.

---

## Contract Decision: ACCEPT with clarifications

### Exact scoring rules for implementation:

When the positive chromosomal structural truth count is zero (`len(truth_mobs) == 0` or `len(truth_jcs) == 0`):

1. **Positive recall:** Return `"NOT_EVALUABLE"` (string sentinel, not numeric 0.0 or 1.0)
   - Denominator: `tp + fn` where `fn = len(truth_events)` → 0
   - Rationale: 0/0 is mathematically undefined; cannot measure sensitivity without positive ground truth

2. **Positive sensitivity:** Return `"NOT_EVALUABLE"` (synonym for recall)

3. **Breakpoint coordinate accuracy (exact match fraction, within-5bp fraction):** Return `"NOT_EVALUABLE"` when `tp == 0`
   - Denominator: `tp` (number of matched events)
   - Rationale: Cannot measure coordinate error distribution with zero matches

4. **TSD accuracy on positive events (length concordance, sequence identity):** Return `"NOT_EVALUABLE"` when `tp == 0`
   - Denominator: `tp` (number of matched MOB events)
   - Rationale: TSD properties are only defined for matched transposition events

5. **Precision:** Follow its own denominator `tp + fp` (test-side calls)
   - If `tp + fp == 0` (tool made no calls), return `"NOT_EVALUABLE"`
   - If `tp + fp > 0`, compute normally: `precision = tp / (tp + fp)` (can be 0.0, defined)
   - Rationale: Precision is test-side metric; mathematically defined as long as tool executed and produced calls

6. **Biological specificity:** Requires independent classification of true negatives via assembly-confirmed zero-difference truth or validated negative control sites
   - Do NOT infer specificity solely from assembly agreement without executing the differential caller
   - Do NOT report specificity = 100% based purely on empty truth; must verify tool produced 0 calls when 0 were expected

7. **Median coordinate error, error quantiles:** Return `"NOT_EVALUABLE"` when the error list is empty (`tp == 0`)

8. **F1-score, F-beta:** Return `"NOT_EVALUABLE"` if either precision or recall is `"NOT_EVALUABLE"`

### Implementation notes:
- Use the string sentinel `"NOT_EVALUABLE"` (all caps, underscore-separated) for consistency with `compare_mob_truth.py:74,80,87-89`
- Do NOT substitute 0% (falsely pessimistic) or 100% (falsely optimistic)
- TSV output: print `NOT_EVALUABLE` verbatim as the value cell
- JSON output: use the string `"NOT_EVALUABLE"`, not null or a missing key
- Aggregation scripts (`aggregate_results.py`) must propagate `NOT_EVALUABLE` as-is and never default it to a number

---

## Inventory claim verification

### Claim 1: `compare_variant_calls.py` currently substitutes 1.0 when a denominator is zero

**Verdict:** CONFIRMED

**Evidence:**
- Line 285: `precision = overall_tp / (overall_tp + overall_fp) if (overall_tp + overall_fp) > 0 else 1.0`
- Line 286: `recall = overall_tp / (overall_tp + overall_fn) if (overall_tp + overall_fn) > 0 else 1.0`
- Lines 289-291: Breakpoint fractions default to `1.0` when `total_tested_bp == 0`
- Lines 334-335 (TSV mode): Per-type `MOB_precision`, `JC_recall` etc. all default to `1.0` via `else 1.0`

**Pathological case confirmed:** Empty truth + N false-positive calls → `recall = 0/(0+0) → 1.0` (100% recall with N false positives)

### Claim 2: `compare_mob_truth.py` already returns NOT_EVALUABLE

**Verdict:** CONFIRMED

**Evidence:**
- Line 74: `precision = "NOT_EVALUABLE"` when `(tp + fp) == 0`
- Line 80: `recall = "NOT_EVALUABLE"` when `(tp + fn) == 0` with comment `"Zero verified positive events in truth: recall is NOT_EVALUABLE (§6, M7-P0a)"`
- Lines 87-89: `exact_fraction`, `within_5bp_fraction`, `median_error` all return `"NOT_EVALUABLE"` when `tp == 0`
- Regression test `test_empty_truth_scoring_returns_not_evaluable` validates this behavior (M7_P0B §14)

---

## Required implementation actions (out of scope for this review, documented for M7-P0c)

1. **`compare_variant_calls.py`** (highest priority): Replace `else 1.0` defaults (lines 285-286, 289-291, 334-335) with `"NOT_EVALUABLE"` sentinel matching `compare_mob_truth.py` semantics; add regression tests mirroring `test_empty_truth_scoring_returns_not_evaluable`

2. **Residual MUMmer artifacts**: Quarantine or delete `truth/prjna884016/SY22A1_vs_SY221.{delta,report,1delta,mdelta,1coords,mcoords,qdiff,snps}` (tracked files with pre-quarantine truncation content); extend `check_validation_assets.sh` to assert their absence (M7_P0B §4 "Stale MUMmer artifacts")

3. **Stale local results**: Annotate or remove `results/prjna884016/{summary.tsv,mob_metrics.tsv,STRUCTURAL_DISCREPANCIES.md}` (gitignored; carry invalidated 100% recall claims) (M7_P0B §12)

---

**SCIENTIFIC_DESIGN_PASS**
