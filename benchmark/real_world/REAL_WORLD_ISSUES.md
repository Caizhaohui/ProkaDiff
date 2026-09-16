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
- **Status:** RESOLVED

