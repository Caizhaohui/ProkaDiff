# ProkaDiff Genome Audit Report

## 1. Executive Summary

**Intended Edit Outcome: ALERT** — 1 intended edit(s) exhibited unexpected/aberrant genomic structures at on-target loci.

### Key Findings Breakdown
- **Total post-edit differential variants:** 1
- **Mobile element (IS) insertions:** 0
- **Genomic rearrangements & large structural deletions:** 1
- **Candidate guide-dependent off-target events:** 1
- **Distal / collateral small variants:** 0

> **Audit Review Recommendation:** 1 high-attention structural or near-target event(s) and 0 review-priority variant(s) require manual inspection.
> *ProkaDiff provides observational audit priorities, not biological safety claims. Final strain suitability depends on specific experimental requirements.*

## 2. Sample & Experiment Context

| Parameter | Declared Value |
| :--- | :--- |
| Starter strain WGS | `starter.fq` |
| Edited strain WGS | `edited.fq` |
| Reference genome | `ref.fa` |
| Editor | `cas9` |
| Guide spacer | `ACGT` |
| PAM profile | `NGG` |
| Parallel threads | `1` |

## 3. Intended Edit Assessment

| Edit ID | Kind | Locus | Status | Left Boundary | Right Boundary | Exp / Obs Size | Diagnostics |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :--- |
| **edit_one** | `snp` | chr:10-10 | **UNEXPECTED_STRUCTURE** | NA | NA | NA / NA | none |

## 4. Genome-wide Differential Summary

| Priority Category | Count | Primary Causes / Characteristics |
| :--- | :---: | :--- |
| **HIGH_ATTENTION** | **1** | Large structural rearrangements, IS transposon insertions, or aberrant target structures |
| **REVIEW** | **0** | Small point mutations, indels, or coding variants requiring review |
| **INFO** | **0** | Confirmed expected edits and neutral background changes |

## 5. Structural and Mobile-element Events

| Variant ID | Type | Locus | Description | Evidence | Priority |
| :--- | :---: | :--- | :--- | :---: | :---: |
| **VAR_0002** | `JC` | chr:20 | Novel sequence junction linked to chr:30 | `RA=NA;MC=NA;JC=NA` | **HIGH_ATTENTION** |

## 6. Candidate Guide-dependent Off-target Events

| Variant ID | Locus | Change | Mismatches | PAM | Distance to Site | Review Priority |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **VAR_0002** | chr:20 | `JC` | 1 | `AGG` | 3 bp | **HIGH_ATTENTION** |

> **Scientific Caution:** Spatial proximity and sequence homology to predicted guide sites indicate candidate association, but do not establish nuclease cleavage causality without experimental controls.

## 7. Distal Small Variants

No distal small variants were detected.

## 8. Analysis Limitations & Cautions

- **Short-Read Structural Resolution:** Variant and junction calling relies on short Illumina read mappings. Extremely repetitive genomic regions or large balanced inversions may not be fully resolved.
- **Observational Association vs. Biological Causality:** Post-edit differential variants cannot be automatically attributed as direct enzymatic consequences of editing without independent controls or multiple clones.
- **CFD Scoring Unavailable:** CFD off-target cleavage scores were omitted because external-oracle parity validation was not enabled for this analysis.

## 9. Technical Provenance

| Parameter | Provenance Detail |
| :--- | :--- |
| ProkaDiff Engine Version | `test` |
| Git Commit SHA | `test` |
| Offtarget Engine Status | `FIXTURE_VALIDATED` |
| CFD Scoring Gate | `DISABLED` |
| Bulge Search Status | `EXACT_UNGAPPED` |
| Analysis Timestamp | `2026-09-24T00:00:00Z` |

