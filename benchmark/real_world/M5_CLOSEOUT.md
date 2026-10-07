# M5 Closeout — BL21 real-data validation

**Verdict: M5 COMPLETE (2026-10-07)**

Authoritative closure record for ProkaDiff milestone M5. Scientific review
concluded **M5 READY TO CLOSE**. This document records validated claims,
provenance caveats, and retained limitations. It does not start M6.

Results under `benchmark/real_world/results/` are gitignored. Counts and
job identities below summarize observed Slurm logs and product validation
artifacts; they are not versioned binary payloads inside `b558371`.

## A. Scope

M5 closed these components only:

| Component | Status |
| --- | --- |
| M5-SR1 Repeat-Ambiguous Junction Recovery | COMPLETE |
| M5-BL21-P1 single-sample reference-relative validation | COMPLETE |
| M5 PEER_COMPARATOR matrix (three directions × G1/G2) | COMPLETE |

Out of scope for this closure: PRJNA1088182 / PRJNA884016 public registry
promotion (M6+), matched biological parent redesign, CAST/IS110, long reads,
and any claim that Cas9 caused the BL21 structural events.

M1–M4 remain frozen. No crate changes accompany this documentation closure.

## B. Frozen scientific invariants

M5 evidence did not contradict the frozen M1–M4 contract:

- `DV1 EventId` is the authoritative public differential identity.
- MC/RA remain evidence-only; they do not receive EventIds.
- Fabricated exact JC/DEL/EventId from unresolved repeat diagnostics is
  forbidden.
- Candidate search states distinguish not-performed, zero-candidate, and
  performed-with-candidates.
- Association is not Cas cleavage causality.
- Intended assessment precedes masking.
- Public schema V2 / seven-product consistency remains the product contract.
- Hidden MOB companion handling and frozen subtraction semantics are
  unchanged.

## C. SR1 validation

**Scientific conclusion.** ProkaDiff preserves high-copy repeat-associated
junction evidence while refusing copy-specific structural calls when
repeat-copy identity is not independently resolved.

| Job | Role | Outcome |
| --- | --- | --- |
| 2815876 | Synthetic `qcpu_18i` SR1 suite | PASS — 20-copy control; 21-copy `SUPPORTED_AMBIGUOUS` / `CANDIDATE` with `exact_supported=false`; reorder determinism; RESOURCE_LIMIT fail-closed tests |
| 2815888 | Focused B21_3_1 vs `GCF_013167015.1` | PASS — lacZ-edge repeat evidence retained as `SUPPORTED_AMBIGUOUS`; multiple feasible copies; `exact_breakpoint_supported=false`; no fabricated exact target JC/DEL |

**Provenance qualification.** Jobs 2815876 and 2815888 were **not** clean
checkouts of `f30f93be5439d2f8eb32f79719fe8482713fcf5a`. Both recorded:

- `HEAD = 08d41edd7f6ce7f6d287fb19cebdedd4e3bcb09c`
- a dirty worktree that already contained the SR1 implementation under test

For job 2815876, the saved worktree `git diff HEAD` SHA256 was
`e6111f92e94cb6c5986d37a2dff07256db4a4f49f08bae223472b1dee10701e7`.

Commit `f30f93b` later **froze** the validated SR1 implementation. Closure
must not claim those two jobs executed a clean `f30f93b` checkout.

`MAX_STAGE2_FAMILY_GEOMETRIES=5000` was exercised by the synthetic overflow
fixture and failed closed. That is not a universal proof of sufficiency.

Closure does **not** claim: exact breakpoint, exact repeat copy, IS1A
identity, or Cas9 causality.

## D. Three-clone single-sample validation

**Scientific conclusion.** All three abnormal BL21 clones show
reference-relative structural abnormality overlapping lacZ. ProkaDiff
retained repeat-associated junction evidence but did not resolve the repeat
family to a single copy. For all three, `exact_breakpoint_supported=false`.

| Sample | Job | Notes |
| --- | --- | --- |
| B21_3_1 | 2815888 | Also the SR1 focused real-data job (dirty `08d41ed` worktree; see §C) |
| B21_4_1 | 2816157 | Explicitly recorded `SR1_COMMIT=f30f93b` |
| B21_2_1 | 2816170 | Explicitly recorded `SR1_COMMIT=f30f93b` |

Observed ProkaDiff patterns (see also
[M5_BL21_P1_RECONCILIATION.md](M5_BL21_P1_RECONCILIATION.md)):

- lacZ-overlapping MC intervals
- `SUPPORTED_AMBIGUOUS` repeat-family evidence at the MC left edge
- multiple feasible copies
- no RESOURCE_LIMIT
- no EventId / DifferentialEvent minted from the repeat diagnostic sidecar

ProkaDiff, breseq, and curated interpretations remain **source-separated**.
breseq copy-specific JC calls and curated IS1/IS1A labels are external
evidence. They were not used to force a ProkaDiff copy choice or EventId.

ProkaDiff did **not** identify an exact repeat-copy breakpoint, an exact
IS1/IS1A copy, or Cas9 causality in these runs.

The bound `MAX_STAGE2_FAMILY_GEOMETRIES=5000` was **not** triggered in the
three BL21 datasets (`overflow_count=0`). Supported claim only: the current
bound was not triggered in the tested BL21 datasets.

## E. Six-cell PEER_COMPARATOR matrix

**Scientific conclusion.** Comparator-relative intended assessment remains
internally consistent under `PEER_COMPARATOR` semantics, even though all
three samples independently contain target-locus abnormalities.

Execution / validation commit:
`b558371a85367bcbb1b97a63647aff3d663f2bd7`.

| Job | Direction | Guides | Slurm | Exit |
| --- | --- | --- | --- | --- |
| 2888235 | B21_3_1 vs B21_4_1 | G1 + G2 | COMPLETED | 0:0 |
| 2888236 | B21_4_1 vs B21_3_1 | G1 + G2 | COMPLETED | 0:0 |
| 2888237 | B21_2_1 vs B21_3_1 | G1 + G2 | COMPLETED | 0:0 |

All six guide cells validated **PASS**. For every cell:

- M2 status = `MISSING`
- intended-related DifferentialEvents = 0
- `matched_event_ids_dv1` empty
- `unexpected_event_ids_dv1` empty
- `false_complete = 0`
- `gate_passed = true`
- candidate search = `PERFORMED_WITH_CANDIDATES`
- association count = 0
- seven M4 products present and product validation PASS

Guide search actually executed:

- G1: 1 candidate site, including a 0-mismatch candidate
- G2: 2 candidate sites, including 0- and 4-mismatch candidates

Guide-site presence and mismatch count do **not** establish cleavage or
editing causality. Zero associations are valid for these `MISSING` products.

### Directional DifferentialEvent counts

| Direction | DifferentialEvents |
| ---: | ---: |
| B21_3_1 vs B21_4_1 | 32 |
| B21_4_1 vs B21_3_1 | 44 |
| B21_2_1 vs B21_3_1 | 54 |

Directional asymmetry is expected and acceptable because subtraction is
directional. A−B and B−A symmetry is not required.

`MISSING` after peer subtraction does **not** retract the single-sample
finding that each clone is abnormal relative to `GCF_013167015.1`.

## F. EventId / product validation

Within the current six-cell matrix only:

- public EventIds use `DV1_<32 hex>`
- no numeric GD IDs used as public identity
- no MC/RA EventIds
- no unresolved association EventIds
- `mutation_offtarget_links.tsv` contains no data rows
- `hidden_mob_companions = 0`
- seven M4 product files are internally consistent
  (`summary.txt`, `report.md`, `edit_outcomes.tsv`,
  `post_edit_variants.tsv`, `unintended.tsv`, `offtarget_sites.tsv`,
  `mutation_offtarget_links.tsv`)

Do not generalize these integrity observations beyond the tested matrix.

## G. Provenance caveats

- Matrix `run_identity.tsv` files record `worktree_dirty=true`. The
  dirty-file list was **not** preserved.
- The current working tree relative to HEAD contains only the intentionally
  untracked file
  `benchmark/real_world/scripts/run_m5_bl21_p1_evidence.sbatch`.
  That does **not** prove it was the only dirty path at matrix runtime.
- Result directories are gitignored and are **not** part of commit
  `b558371`. This closeout summarizes observed validation results and Slurm
  provenance.
- SR1 synthetic / B21_3_1 focused jobs ran on dirty `08d41ed` trees; see §C.
  Jobs 2816157 and 2816170 explicitly recorded SR1 commit `f30f93b`.

History is not rewritten: earlier plan text that said the peer matrix had
not been submitted is superseded by this closeout and by the updated plan
status, not by altering past job logs.

## H. Scientific limitations retained after M5

1. The BL21 experiment lacks a verified biological parent.
2. `PEER_COMPARATOR` is a computational comparator, not a parent-child
   relationship.
3. Cas9 causality is not established.
4. ProkaDiff did not independently resolve the exact repeat copy in these
   BL21 cases.
5. breseq copy-specific JC calls and curated IS1/IS1A interpretations are
   external evidence sources.
6. Guide-site presence and mismatch count do not establish cleavage or
   editing causality.
7. `MAX_STAGE2_FAMILY_GEOMETRIES=5000` is not universally validated.
8. Public benchmark registry work for PRJNA1088182 / PRJNA884016 remains
   outside M5 closure.

## I. Final M5 verdict

**M5 COMPLETE (2026-10-07).**

M5 no longer blocks starting M6. This closeout does **not** begin M6 and
does not claim M6 readiness beyond the removal of the M5 gate.
