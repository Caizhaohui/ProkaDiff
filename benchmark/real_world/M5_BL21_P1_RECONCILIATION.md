# M5-BL21-P1 three-clone reconciliation

Reference for every sample: `GCF_013167015.1`, contig `NZ_CP053602.1`.
lacZ CDS: `complement(334064..337138)`.

These runs are reference-relative single-sample evidence. They are not matched-parent comparisons and not peer-comparator analyses. SR1 implementation freeze commit `f30f93be5439d2f8eb32f79719fe8482713fcf5a`. Execution/documentation commit for this reconciliation `dc977ce7d90954ef1185feb4ddfa6700d99dd34a`. Jobs 2816157 and 2816170 explicitly recorded `SR1_COMMIT=f30f93b`. Job 2815888 (also the SR1 focused real-data run) logged HEAD `08d41ed` with a dirty worktree; see [M5_CLOSEOUT.md](M5_CLOSEOUT.md) §C for the provenance qualification.

Molecule counts below are per diagnostic row. Rows at the same anchor are not added together.

## PROKADIFF

| Field | B21_3_1 | B21_4_1 | B21_2_1 |
| --- | --- | --- | --- |
| Job | 2815888 | 2816157 | 2816170 |
| Output | `B21_3_1_evidence_m5_sr1` | `B21_4_1_evidence_m5_bl21_p1` | `B21_2_1_evidence_m5_bl21_p1` |
| MC interval | 331956–352202 | 307991–352201 | 307061–352201 |
| Exact target JC/DEL | none | none | none |
| Repeat anchor | 331955/− | 307991/− | 307061/− |
| Repeat state | SUPPORTED_AMBIGUOUS | SUPPORTED_AMBIGUOUS | SUPPORTED_AMBIGUOUS |
| Placement / feasible copies | primary 29/29 and 28/28; Stage-2 11/11 and 5/5 | primary 29/29 and 28/28 | primary 29/29 and 28/28; Stage-2 16/16 |
| Distinct molecules (plus/minus) | 17 (13/4); 92 (12/80); Stage-2 3 (2/1) and 4 (3/1) | 11 (5/6); 81 (12/69) | 11 (9/2); 103 (15/88); Stage-2 8 (2/6) |
| Reciprocal evidence | false | false | false on these target rows |
| exact_breakpoint_supported | false | false | false |
| RESOURCE_LIMIT | none | none | none |
| Wall time | 43:27 | 33:09 | 54:17 |
| Peak RSS | ~11.0 GiB | ~11.6 GiB | ~11.2 GiB |

`placement_set_complete=true` on these target rows. No EventId or DifferentialEvent is produced from the repeat diagnostic.

## BRESEQ

Existing oracle files under `analysis_bl21_cas9_wgs/results/sv/<sample>/augmented_full/breseq/output/output.gd`. Not rerun for this record.

| Field | B21_3_1 | B21_4_1 | B21_2_1 |
| --- | --- | --- | --- |
| MC | 331956–352268 | 307992–352239 | 307062–352239 |
| JC | 331955/−1 → 351541/+1 | 307991/−1 → 351541/+1 | 264860/−1 → 307061/−1 |
| JC support | `new_junction_read_count=385` | `new_junction_read_count=389` | `new_junction_read_count=484` |
| UN | 331956–352164 | 307992–352164 | 307062–352164 |

## CURATED

Source: `benchmark/real_world/truth/bl21_curated_truth.tsv`. This is manual interpretation, not a ProkaDiff call.

| Sample | Interpretation |
| --- | --- |
| B21_3_1 | 19.6 kb deletion removing lacZ and ending at downstream IS1A (`HO396_RS01630`) |
| B21_4_1 | 43.5 kb deletion removing lacZ and ending at downstream IS1A (`HO396_RS01630`) |
| B21_2_1 | 45.2 kb complex rearrangement removing lacZ between an upstream IS1 and a downstream IS1 |

## Claim boundaries

| Claim | B21_3_1 | B21_4_1 | B21_2_1 |
| --- | --- | --- | --- |
| Abnormal target-locus coverage loss exists | SUPPORTED BY PROKADIFF | SUPPORTED BY PROKADIFF | SUPPORTED BY PROKADIFF |
| Repeat-associated junction evidence exists at the MC left edge | SUPPORTED BY PROKADIFF | SUPPORTED BY PROKADIFF | SUPPORTED BY PROKADIFF |
| The repeat family is unresolved and more than one copy remains feasible | SUPPORTED BY PROKADIFF | SUPPORTED BY PROKADIFF | SUPPORTED BY PROKADIFF |
| Copy-specific junction coordinate and read support | SUPPORTED ONLY BY BRESEQ | SUPPORTED ONLY BY BRESEQ | SUPPORTED ONLY BY BRESEQ |
| IS1/IS1A identity and the named kilobase rearrangement | SUPPORTED ONLY BY CURATED INTERPRETATION | SUPPORTED ONLY BY CURATED INTERPRETATION | SUPPORTED ONLY BY CURATED INTERPRETATION |
| Exact repeat-copy breakpoint called by ProkaDiff | NOT SUPPORTED | NOT SUPPORTED | NOT SUPPORTED |
| Exact IS1A copy resolved by ProkaDiff | NOT SUPPORTED | NOT SUPPORTED | NOT SUPPORTED |
| Cas9 caused the structural event | NOT SUPPORTED | NOT SUPPORTED | NOT SUPPORTED |
| Every observed genome-wide change was editing-induced | NOT SUPPORTED | NOT SUPPORTED | NOT SUPPORTED |

breseq coordinates and curated IS labels were not used to choose a ProkaDiff copy, change a diagnostic state, or create an EventId.

## ProkaDiff versus breseq

Both tools detect abnormal structure overlapping lacZ. breseq reports a copy-specific, repeat-redundant JC. ProkaDiff keeps family-level evidence because its own placements do not reduce the feasible copy set to one copy. That difference is not recorded here as a ProkaDiff false negative.

| Sample | Class |
| --- | --- |
| B21_3_1 | Biological agreement on the lacZ-overlapping missing interval, with a representation difference and a repeat-copy resolution difference. The breseq other side is 351541/+1. ProkaDiff leaves 5–29 feasible copies and emits no exact JC. |
| B21_4_1 | Same class. Missing-interval left edges agree at 307991/307992. breseq names 351541/+1. ProkaDiff primary rows leave 28 and 29 feasible copies. |
| B21_2_1 | Same class, with a different breseq geometry: 264860/−1 → 307061/−1. ProkaDiff supports the 307061/− anchor as an unresolved family. One Stage-2 interval begins at 351541 and still contains 16 feasible copies; that interval is not a call of the breseq junction. |

No sample is an unexplained discrepancy. No ProkaDiff product contains an exact target JC that lacks support.

## Resource empirical summary

| Sample | Diagnostic rows | CANDIDATE | SUPPORTED_AMBIGUOUS | RESOLVED | RESOURCE_LIMIT | Max placement | Wall | Peak RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| B21_3_1 | 3619 | 3612 | 7 | 0 | 0 | 124 | 43:27 | ~11.0 GiB |
| B21_4_1 | 2817 | 2815 | 2 | 0 | 0 | 176 | 33:09 | ~11.6 GiB |
| B21_2_1 | 4923 | 4908 | 15 | 0 | 0 | 84 | 54:17 | ~11.2 GiB |

`overflow_count` is 0 in all three. The current geometry bound of 5000 was not triggered in these three datasets. This record does not claim that 5000 is sufficient for repeat-rich bacterial genomes in general.

## Peer-comparator matrix purpose

The single-sample question, whether each clone shows a lacZ abnormality relative to `GCF_013167015.1`, is already answered. The peer matrix does not retest that question.

Historical note: when this reconciliation was first written (2026-09-28), the
peer matrix had not yet been submitted. That statement is superseded.

The `PEER_COMPARATOR` matrix was executed and validated (jobs 2888235,
2888236, 2888237; commit `b558371`). Authoritative matrix counts, EventId
checks, provenance caveats, and the M5 COMPLETE verdict are recorded in
[M5_CLOSEOUT.md](M5_CLOSEOUT.md).

Completed directions (each with independent G1 and G2):

- B21_3_1 / B21_4_1
- B21_4_1 / B21_3_1
- B21_2_1 / B21_3_1

`B21_4_1` and `B21_2_1` remain computational comparators in the manifest.
`biological_parent_verified` remains false. Shared lacZ abnormalities can be
removed by subtraction; all six validated cells reported comparator-relative
`MISSING`. The matrix must not require `UNEXPECTED_STRUCTURE`.

## Frozen peer-matrix acceptance

These acceptance rules remain the frozen contract used by the completed
matrix. For each peer run, `MISSING` is valid when all of the following hold:

- no intended-related DifferentialEvent remains after subtraction
- `matched_event_ids_dv1` is empty
- `unexpected_event_ids_dv1` is empty
- no product row contradicts that status

The run is still rejected if any of the following occur:

- a false `COMPLETE`
- an EventId reference that does not resolve to a real differential event
- a fabricated candidate-site measurement
- wording that states Cas9 causality or an exact repeat-copy breakpoint without ProkaDiff evidence
- a differential event that the inputs support and the product silently drops

G1 and G2 are separate runs. One guide's status does not fill in the other guide.

Single-sample reconciliation in this file remains reference-relative evidence
only. M5 as a whole is **COMPLETE**; see [M5_CLOSEOUT.md](M5_CLOSEOUT.md).
M6 has not started.
