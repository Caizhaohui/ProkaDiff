# Case Study: Post-Edit Genome Audit of Cas9-Edited *Escherichia coli* BL21(DE3)

## 1. Executive Summary

This case study documents the real-world validation of **ProkaDiff** on deep whole-genome sequencing (WGS, ~370×–400× coverage) of *Escherichia coli* BL21(DE3) strains subjected to CRISPR-Cas9 genome editing targeting the `lacZ` locus.

In the laboratory, three resequenced clones (`B21_2_1`, `B21_3_1`, `B21_4_1`) presented a puzzling diagnostic phenotype: **complete PCR failure (no amplification band) at the target locus**, while control genomic loci amplified normally.

Reference-relative ProkaDiff evidence addresses whether each resequenced clone differs from `GCF_013167015.1` at `lacZ`. It does not establish a verified parent-child outcome, and it does not establish that Cas9 caused the structure. The M5-BL21-P1 record is [`benchmark/real_world/M5_BL21_P1_RECONCILIATION.md`](../benchmark/real_world/M5_BL21_P1_RECONCILIATION.md).

### Key Findings
- **PROKADIFF (M5-BL21-P1, 2026-09-28):** `B21_3_1` (job 2815888), `B21_4_1` (job 2816157), and `B21_2_1` (job 2816170) each show lacZ-overlapping missing coverage. ProkaDiff detected repeat-associated structural evidence at the left edge of that interval and retained it as `SUPPORTED_AMBIGUOUS`. Multiple repeat copies remain feasible, `exact_breakpoint_supported=false`, and the diagnostic does not emit an exact target JC/DEL or an EventId.
- **Historical peer comparison:** `B21_3_1` was also compared with `B21_4_1` used only as a computational starter (`B21_3_1_diff_rw005`, job 2679261). That report's intended-edit status is `UNEXPECTED_STRUCTURE` for the declared 860 bp deletion. `B21_4_1` is not a verified biological parent. The status is comparator-relative. It is not an exact breakpoint call, and a later peer matrix does not require `UNEXPECTED_STRUCTURE`.
- **BRESEQ and CURATED, kept separate:** breseq reports copy-specific junctions. Curated interpretation describes IS-associated deletions or rearrangements of about 19.6 kb, 43.5 kb, and 45.2 kb. ProkaDiff does not independently identify the exact IS copy.
- **Reference-relative variants:** differences from `GCF_013167015.1` are not, by themselves, editing-induced mutations. Subtracting one abnormal clone from another does not isolate a parent-child post-edit set.
- **Scientific Neutrality:**
  Distal small variants stay observational. This case study does not treat SOS-stress mutagenesis, Cas9 causality, or an exact IS-mediated breakpoint as a ProkaDiff result.

---

## 2. Experimental Design & Locus Topology

| Feature | Coordinates (`NZ_CP053602.1`) | Description |
| :--- | :--- | :--- |
| Target Gene | `334,064 – 337,138` (reverse strand) | `lacZ` beta-galactosidase |
| Spacer 1 (N20-1) | `334,902 – 334,921` (PAM: `CGG`) | Predicted cleavage at 334,918 \| 334,919 |
| Spacer 2 (N20-2) | `335,642 – 335,661` (PAM: `CGG`) | Predicted cleavage at 335,658 \| 335,659 |
| Left Homology Arm (LHA) | `334,576 – 334,875` (300 bp) | Plasmid donor arm |
| Right Homology Arm (RHA) | `335,736 – 336,035` (300 bp) | Plasmid donor arm |
| Intended Clean Deletion | `334,876 – 335,735` (860 bp) | Seamless excision of inter-arm region |
| Flanking PCR Primers | `334,376–334,397` (F), `336,242–336,262` (R) | WT: 1,887 bp; Intended edit: 1,027 bp |

---

## 3. Observed Genomic Events & Structural Resolution

> **Source note (M5-BL21-P1):** junction coordinates and read counts below are
> **BRESEQ** or **CURATED**. ProkaDiff MC intervals and repeat-diagnostic states are
> **PROKADIFF**. ProkaDiff detects repeat-associated structural evidence and does not
> reduce the feasible repeat-copy set to one copy, so it does not emit those exact
> junctions. That difference is a representation and repeat-copy resolution difference,
> not a ProkaDiff false negative.

### Clone `B21_3_1`
- **PROKADIFF:** MC `NZ_CP053602.1:331,956–352,202` (job 2815888). Repeat anchor `331,955/−`, state `SUPPORTED_AMBIGUOUS`, feasible copies 29 and 28 on the two primary rows, `exact_breakpoint_supported=false`. No exact target JC/DEL.
- **BRESEQ:** JC `331,955 (−1)` → `351,541 (+1)`, 385 reads. MC `331,956–352,268`.
- **CURATED:** 19.6 kb deletion removing `lacZ` and interpreted as ending at downstream IS1A. ProkaDiff does not make that copy or mechanism call.

### Clone `B21_4_1`
- **PROKADIFF:** MC `NZ_CP053602.1:307,991–352,201` (job 2816157). Repeat anchor `307,991/−`, state `SUPPORTED_AMBIGUOUS`, feasible copies 29 and 28 on the two primary rows, `exact_breakpoint_supported=false`. No exact target JC/DEL.
- **BRESEQ:** JC `307,991 (−1)` → `351,541 (+1)`, 389 reads. MC `307,992–352,239`.
- **CURATED:** 43.5 kb deletion removing `lacZ` and interpreted as ending at downstream IS1A. ProkaDiff does not make that copy or mechanism call.

### Clone `B21_2_1`
- **PROKADIFF:** MC `NZ_CP053602.1:307,061–352,201` (job 2816170). Repeat anchor `307,061/−`, state `SUPPORTED_AMBIGUOUS`, feasible copies 29, 28, and a Stage-2 row with 16 copies, `exact_breakpoint_supported=false`. No exact target JC/DEL.
- **BRESEQ:** JC `264,860 (−1)` → `307,061 (−1)`, 484 reads. MC `307,062–352,239`.
- **CURATED:** 45.2 kb complex rearrangement removing `lacZ` between flanking IS1 elements. ProkaDiff does not resolve that junction or assign the IS copies.

---

## 4. Release Gate Verification

- **Criterion:** `false_complete == 0` on curated real-world ground truth.
- **Evaluation:**
  - Standard caller without boundary validation might misinterpret disappearance of wild-type reads as intended deletion.
  - ProkaDiff checks whether an intended interval still has supporting evidence. On these clones that evidence is missing coverage plus unresolved repeat-family diagnostics, not a verified pair of exact junction coordinates.
  - Result: an earlier differential audit run on `B21_3_1` (edited) vs `B21_4_1` (used
    as starter) returned **`MISSING`**, not `UnexpectedStructure`
    (`benchmark/real_world/results/bl21_user/B21_3_1_diff/report.md`, `summary.txt`).
    Root cause: two independent bugs. RW-004 — the classifier's subtract/mask pipeline
    (`prokadiff-classify::classify`) filtered out `MC` (missing-coverage) records before
    they ever reached the intended-edit assessment. RW-005 — GenBank-derived contig ids
    dropped the `.1` version suffix that the `--intended` TSV uses, so even after fixing
    RW-004 the assessment's `seq_id` lookup still failed to match on this dataset. Both
    fixes are now implemented (`crates/prokadiff-classify/src/classify.rs`,
    `crates/prokadiff-evidence/src/fasta.rs` + `fasta/genbank.rs`) and **verified via a
    full end-to-end CLI re-run** (Slurm job 2679261, 2026-09-16):
    `benchmark/real_world/results/bl21_user/B21_3_1_diff_rw005/report.md` now reports
    `Intended Edit Outcome: ALERT`, with `edit_1` assessed as `UNEXPECTED_STRUCTURE`.
    `false_complete == 0` holds for that historical peer run. The report status is
    comparator-relative and does not identify an exact on-target breakpoint. Both
    RW-004 and RW-005 are `RESOLVED` in `REAL_WORLD_ISSUES.md`.
    `B21_2_1` has reference-relative ProkaDiff evidence (job 2816170) and no peer
    release-gate run. Any `UNEXPECTED_STRUCTURE` label for this clone remains a curated
    interpretation, not a ProkaDiff differential status.

---

## 5. Benchmarking vs Oracle `breseq 0.40.2`

> `B21_2_1` reference-relative evidence is job 2816170 (§3). The genome-wide JC/MC
> counts below are the earlier P1-3 comparison, not the M5-BL21-P1 repeat-diagnostic
> inventory. No `--starter`/`--edited` differential audit has been run with `B21_2_1`.

**Correction (M5-SR1 / M5-BL21-P1):** an earlier draft said ProkaDiff called
`331,955 / 351,541` with exact breseq concordance. That JC is a breseq/curated
coordinate. Current ProkaDiff output still has no exact target JC. After SR1 the
repeat-associated reads are retained as `SUPPORTED_AMBIGUOUS` family evidence rather
than discarded for having more than 20 placements. The missing exact JC is the
unresolved-copy representation, not evidence that the reads disappeared. RW-004 remains
the historical record of the intended-edit classifier bug; it is not a requirement that
ProkaDiff reproduce breseq's copy-specific junction.

| Metric | ProkaDiff Evidence Engine | `breseq 0.40.2` Oracle | Concordance |
| :--- | :--- | :--- | :--- |
| `B21_3_1` on-target Junction | family-level `SUPPORTED_AMBIGUOUS` at 331,955/−; no exact JC | `331,955 / 351,541` (385 reads) | **Repeat-copy resolution difference** |
| `B21_4_1` on-target Junction | family-level `SUPPORTED_AMBIGUOUS` at 307,991/−; no exact JC | `307,991 / 351,541` (389 reads) | **Repeat-copy resolution difference** |
| Background SNP Filter | Reference-relative calls; clone subtraction is not a verified parent-child filter | Requires manual `gdtools SUBTRACT` | **Different comparison** |
| Target Locus Status | Historical peer run `B21_3_1` vs `B21_4_1` reported `UNEXPECTED_STRUCTURE`; single-sample P1 does not emit that status | Raw genome diff records only | **Peer status, not an exact call** |

Full-genome comparison against the `B21_3_1` breseq oracle
(`benchmark/real_world/results/bl21_user/B21_3_1_evidence/variant_metrics.tsv`) confirms
the evidence engine's JC/MC output is unreliable well beyond just this one locus — it
currently emits a large number of additional JC/MC calls elsewhere in the genome that do
not correspond to any oracle call, and separately fails to reproduce the on-target call
itself:

### Initial Baseline (Pre-Fix, 2026-09-16)

| Sample | Type | ProkaDiff calls | breseq calls | TP | FP | FN | Precision | Recall |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `B21_3_1` | JC | 146 | 15 | 6 | **140** | 9 | 0.0411 | 0.4000 |
| `B21_3_1` | MC | 75 | 6 | 2 | **73** | 4 | 0.0267 | 0.3333 |
| `B21_4_1` | JC | 168 | 11 | 4 | **164** | 7 | 0.0238 | 0.3636 |
| `B21_4_1` | MC | 79 | 8 | 3 | **76** | 5 | 0.0380 | 0.3750 |
| `B21_2_1` | JC | 134 | 18 | 6 | **128** | 12 | 0.0448 | 0.3333 |
| `B21_2_1` | MC | 74 | 7 | 3 | **71** | 4 | 0.0405 | 0.4286 |

### Post-P1-3 Secondary Filter Results (2026-09-17/18)

After implementing the depth-aware secondary filter (`req_reads = max(3, round(depth * 0.05))`) in `prokadiff-evidence` and filtering breseq `reject=` candidates in evaluation harness:

| Sample | Type | ProkaDiff calls | breseq truth | TP | FP | FN | Precision | Recall | Change vs Baseline |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `B21_3_1` | JC | 35 | 10 | 6 | **29** | 4 | **0.1714** | **0.6000** | FP: 140 → 29 (**-79.3%**), Precision: 4.1% → 17.1% |
| `B21_3_1` | MC | 75 | 6 | 2 | **73** | 4 | 0.0267 | 0.3333 | 73/73 FP explained by breseq UN |
| `B21_4_1` | JC | 42 | 10 | 4 | **38** | 6 | **0.0952** | **0.4000** | FP: 164 → 38 (**-76.8%**), Precision: 2.4% → 9.5% |
| `B21_4_1` | MC | 79 | 8 | 3 | **76** | 5 | 0.0380 | 0.3750 | 76/76 FP explained by breseq UN |
| `B21_2_1` | JC | 43 | 12 | 5 | **38** | 7 | **0.1163** | **0.4167** | FP: 128 → 38 (**-70.3%**), Precision: 4.5% → 11.6% |
| `B21_2_1` | MC | 74 | 7 | 3 | **71** | 4 | 0.0405 | 0.4286 | 71/71 FP explained by breseq UN |

For comparison, SNP calling remains high-accuracy (precision 0.94–0.98 / recall ≥0.99) on each of these samples (`benchmark/real_world/results/bl21_user/{B21_3_1,B21_4_1,B21_2_1}_evidence_p1_3/variant_metrics.tsv`).
The depth-aware secondary filter is the historical RW-002 result: it removed more than 75% of genome-wide JC calls in that comparison without dropping the breseq matches counted as true positives in the table above. Those tables score exact JC coordinate agreement. They do not show that ProkaDiff resolved the lacZ repeat copy. For MC, every false positive in that scoring overlaps a breseq `UN` interval — 73/73 (`B21_3_1`), 76/76 (`B21_4_1`), 71/71 (`B21_2_1`) — which RW-003 treats as a metric-fairness gap. The historical peer report's `UNEXPECTED_STRUCTURE` is a comparator-relative status with `false_complete == 0`, not an exact breakpoint identification.

---

## 6. Recommendations for Engineering Workflows

1. **Do not rely solely on negative PCR:** PCR non-amplification must not be assumed to be PCR failure; gross chromosomal deletions up to 45 kb can eliminate primer binding sites entirely.
2. **Include junction-spanning diagnostic primers:** If a follow-up assay uses `331,955 / 351,541`, that coordinate pair is a breseq/curated example for `B21_3_1`, not a ProkaDiff exact breakpoint.
3. **Keep IS interpretation sourced:** curated notes describe IS-associated rearrangements. ProkaDiff's current evidence supports unresolved repeat-family structure at the target, not a demonstrated Cas9-driven IS1A mechanism.
