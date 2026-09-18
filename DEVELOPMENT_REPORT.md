# ProkaDiff Development & Repositioning Report

## Post-edit Genome Audit for Engineered Prokaryotes (Phases A – G)

**Date:** 2026-09-16  
**Repository:** [ProkaDiff](https://github.com/Caizhaohui/ProkaDiff)  
**Start Commit:** `ac30888`  
**End Commit:** *(Current sprint on main)*  

---

## 1. Executive Summary of Repositioning

ProkaDiff has successfully transitioned from a narrow "CRISPR off-target predictor" into a WGS-first **Post-edit Genome Audit & Verification Framework for Engineered Prokaryotes**.

The core operational principle is now:
```text
Observed genome changes  ──►  Edit verification  ──►  Genome-wide annotation  ──►  Mechanistic association  ──►  Audit Reporting
```
Rather than predicting off-target sites and hunting for matching mutations, ProkaDiff first observes all differential genomic changes between the starter strain and edited clone (`Post-edit Differential Variants`), rigorously verifies declared on-target editing outcomes, annotates mobile genetic elements and structural variants, and generates both machine-readable audit tables and an executive human-readable Markdown report (`report.md`).

---

## 2. Implemented Phases & Tasks

### Phase A: Product & Terminology Realignment
- Updated `README.md`, `README.zh.md`, and CLI `--help` texts to reflect the Post-edit Genome Audit mission.
- Created `docs/scientific_model.md` defining strict scientific evidence levels (Level 0: Observation $\to$ Level 1: Contextual Association $\to$ Level 2: Mechanistic Hypothesis $\to$ Level 3: Experimentally Supported Causality).
- Locked scientific terminology:
  - Banned unjustified claims of "confirmed CRISPR off-target" from pairwise WGS alone; replaced with `candidate guide-dependent off-target`.
  - Replaced "DSB-induced mutations" with `collateral genome change` / `distal small variant`.
  - Legacy three-class taxonomy (`structural`, `near_homolog`, `scattered_snv`) marked as backward-compatible `legacy_class`.

### Phase B: Intended Edit Verification Overhaul
- Overhauled `crates/prokadiff-classify/src/intended.rs`:
  - Fixed critical bug where `e.fields.first().parse()` erroneously treated `seq_id` as mutation ID; replaced with real `e.id`.
  - Implemented multi-state `IntendedEditStatus`: `Complete`, `Partial`, `Missing`, `UnexpectedStructure`.
  - Added `BoundaryAssessment` (exact vs fuzzy tolerance $\le 2$ bp).
  - Redesigned cassette integration verification: requiring both left and right junctions (2 JCs $\to$ `Complete`, 1 JC $\to$ `Partial`, $\ge 3$ or aberrant $\to$ `UnexpectedStructure`).

### Phase C: Unified Domain Models & AuditResult
- Created `crates/prokadiff-classify/src/audit.rs`:
  - `OriginStatus`: `StarterBackground`, `Intended`, `PostEditDifferential`.
  - `SizeClass`: `Small`, `Structural`.
  - `IntendedRelation`: `None`, `Expected`, `PartOfExpectedEdit`, `UnexpectedAtOnTargetLocus`.
  - `GuideRelation`: `None`, `OnTarget`, `CandidateOffTarget`.
  - `ReviewPriority`: `Info`, `Review`, `HighAttention` (transparent review prioritization, never biological safety claims).
  - `EvidenceSummary`: RA, MC, JC flags, read depth, and coverage.
  - Single unified aggregate root: `AuditResult`.

### Phase D: Structured Machine-Readable Deliverables
- Implemented `crates/prokadiff-report/src/tables.rs`:
  - `edit_outcomes.tsv`: per-declared-edit verification audit (boundaries, sizes, diagnostics).
  - `post_edit_variants.tsv`: comprehensive multi-dimensional annotation matrix with stable IDs (`VAR_0001`...`VAR_NNNN`), gene annotations, and IS element metadata.
  - `provenance.tsv`: verifiable analysis provenance, tool versions, and oracle validation gates.
  - Maintained 100% backward compatibility for `unintended.tsv` and `summary.txt`.

### Phase E: Human-Readable Markdown Audit Report
- Implemented `crates/prokadiff-report/src/markdown.rs`:
  - Generates deterministic, reproducible `report.md` structured across 10 standard sections:
    1. Executive Summary
    2. Sample & Experiment Context
    3. Intended Edit Assessment
    4. Genome-wide Differential Summary
    5. Structural and Mobile-element Events
    6. Candidate Guide-dependent Off-target Events (with mandatory scientific caution banner)
    7. Distal Small Variants
    8. Analysis Limitations & Cautions
    9. Technical Provenance
  - Implemented snapshot fixture tests in `crates/prokadiff-report/tests/audit_report_fixtures.rs` (clean edit, partial cassette, IS1 insertion, candidate off-target, unvalidated CFD gate).

### Phase F: Structural & Mobile Element Audit Enhancement
- Extended GenBank feature parsing in `crates/prokadiff-evidence/src/fasta.rs` (`parse_genbank_features`) for CDS, tRNA, rRNA, gene.
- Implemented `find_gene_annotation` in `crates/prokadiff-classify/src/audit.rs` supporting overlapping CDS/tRNA/rRNA and intergenic flanking gene distance calculation (`+X bp from geneA, -Y bp to geneB`).
- Implemented IS family classification (`classify_is_family`) recognizing IS1, IS2, IS3 (IS150/IS911), IS4 (IS186/IS10/IS50), IS5 (IS903/IS1182), IS30, IS110, Tn3, Tn7 families.
- Extracted Target Site Duplications (TSD) directly from GD duplication sizes / repeat sequences.
- Rendered gene and IS family contexts in `report.md` and `post_edit_variants.tsv`.

### Phase G: Multi-Clone Batch Mode Specification
- Authored `docs/batch_mode_spec.md`:
  - Defined one-starter-versus-$N$-clones batch architecture.
  - Specified manifest table schema (`clone_id`, `r1`, `r2`, `intended_tsv`, `spacer`, `pam`).
  - Defined data models for `RecurrenceStatus` (`SharedAll`, `SharedSubset`, `CloneSpecific`) and cohort aggregation.
  - Outlined `batch_report.md` layout adhering to non-prescriptive policy (Section 66).

---

## 3. New Data Structures

```rust
// In prokadiff-classify::audit
pub struct AuditResult {
    pub sample: SampleMetadata,
    pub intended_edits: Vec<IntendedEditAssessment>,
    pub variants: Vec<AnnotatedVariant>,
    pub provenance: AnalysisProvenance,
}

pub struct AnnotatedVariant {
    pub variant_id: String,
    pub entry: GdEntry,
    pub origin_status: OriginStatus,
    pub size_class: SizeClass,
    pub intended_relation: IntendedRelation,
    pub guide_relation: GuideRelation,
    pub mobile_element_relation: Option<MobileElementAnnotation>,
    pub repeat_relation: Option<RepeatAnnotation>,
    pub gene_annotation: Option<GeneAnnotation>,
    pub evidence: EvidenceSummary,
    pub review_priority: ReviewPriority,
    pub legacy_class: Option<MutationClass>,
}

pub struct AnnotatedFeature {
    pub seq_id: String,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub feature_type: String,
    pub locus_tag: Option<String>,
    pub gene_name: Option<String>,
    pub product: Option<String>,
}
```

---

## 4. Output Artifacts Inventory

When running `prokadiff`:
```text
outdir/
├── report.md                  # Executive human-readable markdown audit report
├── summary.txt                # Legacy/shell key-value summary
├── edit_outcomes.tsv          # Intended edit verification matrix
├── post_edit_variants.tsv     # Unified multi-dimensional variant annotations
├── unintended.tsv             # Backward-compatible 3-class differential mutations
├── offtarget_sites.tsv        # Scanned guide-homologous reference sites
├── mutation_offtarget_links.tsv # Associations between variants and candidate sites
├── starter.gd                 # Starter strain breseq-compatible GD
├── edited.gd                  # Edited clone breseq-compatible GD
└── provenance.tsv             # Verifiable technical execution provenance
```

---

## 5. Sample `report.md` Excerpt

```markdown
# ProkaDiff Genome Audit Report

## 1. Executive Summary

- **Intended Edit Outcome: PASS** (1/1 declared edits confirmed)
- **Total post-edit differential variants:** 1
- **Structural and mobile-element events:** 1
- **Candidate guide-dependent off-target events:** 0
- **Distal / collateral small variants:** 0

> **Audit Review Recommendation:** 1 high-attention structural or near-target event(s) and 0 review-priority variant(s) require manual inspection.
> *ProkaDiff provides observational audit priorities, not biological safety claims. Final strain suitability depends on specific experimental requirements.*

## 5. Structural and Mobile-element Events

| Variant ID | Type | Locus | Description | Evidence | Priority |
| :--- | :---: | :--- | :--- | :---: | :---: |
| **VAR_0001** | `MOB` | NC_000913.3:3445210 | Mobile element IS1 insertion (TSD=9bp) | `RA=0;MC=0;JC=1` | **HIGH_ATTENTION** |

## 8. Analysis Limitations & Cautions

- **Short-Read Structural Resolution:** Variant and junction calling relies on short Illumina read mappings. Extremely repetitive genomic regions or large balanced inversions may not be fully resolved.
- **Observational Nature of Pairwise WGS:** Post-edit differential variants cannot be automatically attributed to the editing procedure without appropriate experimental controls.
- **Candidate Guide Associations:** Spatial proximity and sequence homology to predicted guide sites indicate candidate association, but do not establish nuclease cleavage causality without experimental controls.
- **CFD Scoring Unavailable:** CFD scoring was disabled or marked unavailable because external-oracle parity validation was not enabled for this run.
```

---

## 6. Backward Compatibility Assurance

- `unintended.tsv` remains completely unchanged in columns, schema, and values.
- `summary.txt` retains all historical fields (`starter_vs_ref`, `edited_vs_ref`, `post_subtract`, `unintended_mutations`, etc.).
- CLI flags (`--starter`, `--edited`, `--ref`, `--intended`, `--editor`, `--spacer`, `--pam`) operate identically while seamlessly producing the enhanced deliverables.

---

## 7. Verification & Test Suite Status

Executed on compute node partition `qcpu_18i`:
- `cargo fmt --all --check`: **PASSED**
- `cargo clippy --workspace --all-targets -- -D warnings`: **PASSED (0 warnings)**
- `cargo test --workspace`: **ALL PASSED (238 tests across all crates, 0 failures, 1 ignored cluster-only test)**

---

## 8. Known Limitations & Remaining Roadmap

1. **CFD / Hsu Scoring:** Intentionally gated as `DISABLED_UNVALIDATED_ORACLE` until empirical external oracle validation is achieved against published datasets.
2. **Multi-Clone Batch Mode Execution:** Full pipeline execution defined in `docs/batch_mode_spec.md`; ready for runtime CLI implementation in subsequent sprints.
3. **External DNA Integration:** Donor / vector backbone integration tracking defined in data models, awaiting dedicated plasmid junction modules.
