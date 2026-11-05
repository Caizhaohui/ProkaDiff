# M7 stage execution plan

**Governance:** `PROKADIFF_MULTI_AI_V1_3` · **Compute authorization:** `M7_STAGE_COMPUTE_AUTHORIZED`  
**Dedicated worktree:** `/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff-m7-worktree` (`m7-autonomous-dev` @ `1a0c7b5`)  
**Main worktree:** preserve M6-P2A dirty files; do not mix commits.

## Resource classification by workstream

| Workstream | Phase | Workload class | Partition (preferred) | Typical resources | Notes |
| --- | --- | --- | --- | --- | --- |
| **P0a/P0b** | Done | Governance + pytest | Login (historical) / compute preferred | — | Fail-closed quarantine complete |
| **P1 Negative concordance smoke** | Next | WGS differential smoke | **`qcpu_18i`** | 8–16 CPU, 32–48 GB, ≤8 h | SY22A1 vs SY221 vs `CP106998.1`; expect 0 chromosomal calls |
| **P1 FASTQ staging** | Queued | Download + provenance | **`qcpu_18i`** or data node policy | Bounded I/O | Read-only on source; manifest + SHA-256 |
| **Benchmark scripts / unit tests** | Ongoing | pytest, cargo (workspace) | **`qfree`** or active **`qcpu_18i`** alloc | 4–8 CPU, 16–32 GB, ≤2 h | No breseq on login |
| **Comparator / breseq oracle** | Later M7 | Structural compare | **`qcpu_18i`** | 8–16 CPU, 32–48 GB | breseq 0.40.x oracle only |
| **Registry / scoring validation** | Ongoing | Deterministic scripts | Compute allocation | Light CPU | Registry stays `PENDING`; no auto-promotion |
| **Interactive dev** | Ongoing | OpenCode / Grok / Cursor | **`qfree`** if nodes available else existing **`qcpu_18i`** interactive | 4 CPU, 12 GB, 2 h | Controller orchestrates |

## Estimated compute cost (order-of-magnitude)

| Item | CPU-hours (est.) | Cumulative M7 tracker |
| --- | ---: | --- |
| M7-P0b regression on compute (smoke) | &lt;0.01 | ~0.01 (verified run 2026-10-28) |
| M7-P1 negative smoke (planned) | 4–16 | TBD after `sacct` |
| Full M7 remaining (planned envelope) | 50–200 | Review at **1,024 CPU-h** threshold |

Track actuals from **`sacct`** into HANDOFF; do not invent CPU-hours.

## Partition selection algorithm (summary)

1. Read-only `sinfo` / `scontrol show partition` / `squeue -u $USER`.
2. Confirm group membership (`hpcuser` → **qfree** candidate; **qtib** → **qcpu_18i** candidate).
3. Pick cheapest eligible partition for workload class (table above).
4. Enforce standing envelope (≤16 CPU, ≤48 GB, ≤8 h, ≤2 concurrent Controller jobs).
5. Record job ID, partition, and allocation in HANDOFF before treating results as valid.

## Scientific gates (unchanged)

- M7-P0a/P0b frozen; PRJNA884016 **no** positive chromosomal MOB/JC truth.
- Slurm **COMPLETED** ≠ biological validation; Luna + deterministic tests for acceptance.

## Full M7 submilestone graph (remaining)

| Phase | ID | Deliverable | Acceptance |
| --- | --- | --- | --- |
| P1a | `m7_p1_fastq_stage` | FASTQ + provenance for SRR21706776/95 | Slurm COMPLETE, checksums in TSV |
| P1b | `m7_p1_negative_concordance` | ProkaDiff differential smoke | 0 calls in unintended/post_edit + empty residual GD |
| P2 | `m7_p2_breseq_oracle` | Optional breseq 0.40.x comparator run | Same-job evidence under results/ |
| P3 | `m7_p3_compare_metrics` | compare_* scripts vs zero truth | NOT_EVALUABLE for sensitivity; precision rules |
| P4 | `m7_p4_independent_review` | Codex Sol/Terra on results + diff | ACCEPT |
| P5 | `m7_p5_luna_validation` | Luna + deterministic gates | PASS or documented blocker |
| Close | `m7_closeout` | `M7_CLOSEOUT.md` | COMPLETE or COMPLETE_WITH_LIMITATIONS |

**Original objective (DEVELOPMENT_PLAN):** PRJNA884016 assembly-truth reconciliation with ProkaDiff calls and structural comparison — **not** positive SV sensitivity on chromosome.

## Next executable task

**Task ID:** `m7_p1_fastq_stage` — Slurm job **3049255** RUNNING (`SRR21706776`).
