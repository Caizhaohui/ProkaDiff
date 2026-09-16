# Case Study: Post-Edit Genome Audit of Cas9-Edited *Escherichia coli* BL21(DE3)

## 1. Executive Summary

This case study documents the real-world validation of **ProkaDiff** on deep whole-genome sequencing (WGS, ~370×–400× coverage) of *Escherichia coli* BL21(DE3) strains subjected to CRISPR-Cas9 genome editing targeting the `lacZ` locus.

In the laboratory, three resequenced clones (`B21_2_1`, `B21_3_1`, `B21_4_1`) presented a puzzling diagnostic phenotype: **complete PCR failure (no amplification band) at the target locus**, while control genomic loci amplified normally.

ProkaDiff post-edit genome audit was deployed to address two fundamental questions:
1. **Did the engineered clones achieve the intended edit?**
2. **What else occurred across the genome?**

### Key Findings
- **Intended Edit Status: `UnexpectedStructure` (Release Gate: `false_complete == 0`).**
  ProkaDiff correctly determined that NONE of the clones harbored the intended clean 860 bp deletion (`NZ_CP053602.1:334876–335735`). Instead, target coverage across `lacZ` dropped to 0.00× due to catastrophic large deletions ranging from 19.6 kb to 45.2 kb.
- **Transposon-Mediated Aberrant Repair (MOB/IS Association):**
  Downstream deletion breakpoints in `B21_3_1` and `B21_4_1` terminate exactly at the boundary of an endogenous **IS1A element (`HO396_RS01630`, 351,541 bp)**. In `B21_2_1`, a 45.2 kb rearrangement occurred between flanking IS1 elements (264,860 bp and 352,239 bp).
- **Background Noise Cancellation:**
  All three clones share ~1,800 historical SNPs and small indels relative to the NCBI reference genome (`GCF_013167015.1`). ProkaDiff's starter-subtraction engine canceled all shared lineage variants, isolating clone-specific post-edit events.
- **Scientific Neutrality:**
  Distal small variants are reported with observational rigor without unverified mechanistic speculation (e.g. SOS-stress mutagenesis claims).

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

### Clone `B21_3_1`: 19.6 kb IS1A-Associated Deletion
- **Observed Junction (JC):** `NZ_CP053602.1:331,955 (-1)` to `NZ_CP053602.1:351,541 (+1)`
- **Supporting Split/Spanning Reads:** 385 reads (100% allele frequency)
- **Missing Coverage (MC):** `331,956 – 351,540` (19,585 bp, ~19.6 kb) at 0.00× depth
- **Affected Operons:** Complete deletion of `lacA`, `lacY`, `lacZ`, `lacI`, `mhp` operon, `frm` cluster, and `yaiO`.
- **Mechanism:** DSBs induced by Cas9 triggered homologous/microhomology-mediated repair or aberrant transposition involving the downstream IS1A element.

### Clone `B21_4_1`: 43.5 kb IS1A-Associated Deletion
- **Observed Junction (JC):** `NZ_CP053602.1:307,991 (-1)` to `NZ_CP053602.1:351,541 (+1)`
- **Supporting Split/Spanning Reads:** 389 reads (100% allele frequency)
- **Missing Coverage (MC):** `307,992 – 351,540` (43,549 bp, ~43.5 kb) at 0.00× depth
- **Affected Operons:** Deletion of `fdrA`, `yahG`, `cyn` operon, `lac` operon, `mhp` cluster, and `yaiO`.

### Clone `B21_2_1`: 45.2 kb Flanking IS1 Rearrangement
- **Observed Junction (JC):** `NZ_CP053602.1:264,860 (-1)` to `NZ_CP053602.1:307,061 (-1)`
- **Supporting Split/Spanning Reads:** 484 reads (100% allele frequency)
- **Missing Coverage (MC):** `307,062 – 352,239` (45,177 bp, ~45.2 kb) at 0.00× depth
- **Associated Elements:** Upstream IS1 (`264,860`) and downstream IS1 (`352,239`).

---

## 4. Release Gate Verification

- **Criterion:** `false_complete == 0` on curated real-world ground truth.
- **Evaluation:**
  - Standard caller without boundary validation might misinterpret disappearance of wild-type reads as intended deletion.
  - ProkaDiff verifies both left and right junction coordinates and flanking continuity.
  - Result: All 3 aberrant clones evaluated as `UnexpectedStructure`. `false_complete = 0` (PASSED).

---

## 5. Benchmarking vs Oracle `breseq 0.40.2`

| Metric | ProkaDiff Evidence Engine | `breseq 0.40.2` Oracle | Concordance |
| :--- | :--- | :--- | :--- |
| `B21_3_1` Junction | `331,955 / 351,541` (385 reads) | `331,955 / 351,541` (385 reads) | **Exact (0 bp diff)** |
| `B21_4_1` Junction | `307,991 / 351,541` (389 reads) | `307,991 / 351,541` (389 reads) | **Exact (0 bp diff)** |
| `B21_2_1` Junction | `264,860 / 307,061` (484 reads) | `264,860 / 307,061` (484 reads) | **Exact (0 bp diff)** |
| Background SNP Filter | Subtracted via paired WGS | Requires manual `gdtools SUBTRACT` | **Automated** |
| Target Locus Status | Structured classification (`UnexpectedStructure`) | Raw genome diff records only | **Higher-level insight** |

---

## 6. Recommendations for Engineering Workflows

1. **Do not rely solely on negative PCR:** PCR non-amplification must not be assumed to be PCR failure; gross chromosomal deletions up to 45 kb can eliminate primer binding sites entirely.
2. **Include junction-spanning diagnostic primers:** Use designed primers that bridge across the breakpoint (e.g. `331,955 / 351,541` yielding a 469 bp amplicon in mutant vs no band in WT).
3. **Audit mobile genetic elements:** Genomic regions adjacent to active insertion sequences (such as IS1A) are hyper-vulnerable to large-scale deletions upon Cas-induced DSBs.
