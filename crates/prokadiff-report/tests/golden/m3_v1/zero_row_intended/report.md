# ProkaDiff Genome Audit Report

## 1. Executive Summary

**Intended Edit Outcome:** *Not Specified* (Analysis run in untargeted differential mode).

### Key Findings Breakdown
- **Total post-edit differential variants:** 0
- **Mobile element (IS) insertions:** 0
- **Genomic rearrangements & large structural deletions:** 0
- **Candidate guide-dependent off-target events:** 0
- **Small variants with unassessed guide proximity:** 0

> **Audit Conclusion:** No additional post-edit differential variants were detected outside declared edits.

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

*No intended edits were declared. To verify targeting outcomes, supply an `--intended` TSV table.*

## 4. Genome-wide Differential Summary

No post-edit differential variants were detected in the edited strain relative to the starter strain.

## 5. Structural and Mobile-element Events

No post-edit differential structural or mobile-element events were detected.

## 6. Candidate Guide-dependent Off-target Events

Candidate guide-site search was not performed; no off-target association can be assessed.

## 7. Small Variants With Unassessed Guide Proximity

No small post-edit differential variants were detected; guide-site proximity was not assessed.

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

