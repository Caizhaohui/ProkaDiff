# M6-P3c.2: Pilot Execution and Orthogonal Structural Reconciliation Report (SRR28349244)

## 1. Executive Summary

As authorized by task `M6-P3C2-GROK-PILOT-EXECUTION`, this report documents the execution of the independent `breseq 0.40.2` comparator pilot for sample `Cas9_1` (run accession `SRR28349244`, BioProject `PRJNA1088182`) on Slurm partition `qcpu_18i`. All operations strictly adhered to clean-room development rules (AGENTS.md Rules 4, 5, 11, and 12).

The breseq job executed successfully to completion in **12 minutes 48 seconds** (MaxRSS: 1.55 GB) without errors. Small-variant mutation calls (7 SNPs, 2 small INS) exhibited **100% coordinate and base identity** between ProkaDiff and breseq. Structural mutation promotion demonstrated robust biological concordance: 1/1 large deletions (82 bp DEL at coordinate 3,815,859) and 3/3 mobile element insertions (IS5 at 889,647, IS1 at 1,873,031, IS1 at 3,799,683) were mutually identified. Under the decoupled 4-attribute MOB concordance model, two MOBs match with exact TSD duplication lengths (`EXACT_MOB_MATCH`), while MOB 147 at 1,873,031 displays a 1 bp microhomology slip (`dup=8` vs `dup=9`), classified strictly as `CONCORDANT_MOB_1BP_TSD_RESOLVED`.

Orthogonal 4-dimensional classification of the 49 ProkaDiff junction candidates revealed:
- **5 junctions** were corroborated with exact dual-breakend coordinate identity ($D=0$, `EXACT_MATCH`).
- **6 junctions** were corroborated at the repeat-family level (`CORROBORATED_AT_FAMILY_LEVEL_UNRESOLVED_AT_LOCUS`, Section 3.J Case 6), corresponding to the constituent junction edges of the 3 mobile element insertions.
- **1 junction** is a near-miss boundary candidate (breseq 39 <-> ProkaDiff 123, $D=10$ bp, breakend 1 delta=6 bp) arising from a 16 bp microhomology window.
- **37 junctions** were uncalled by breseq (`NOT_CORROBORATED`), comprising 4 multi-copy repeat chimeras with massive read depth suppressed by breseq repeat filters and 33 low/moderate support candidates below breseq's frequency thresholds.
- All 17 breseq JCs and all 49 ProkaDiff JCs are completely accounted for in a bidirectional bijective accounting matrix with zero duplicate assignments.

Greedy bipartite matching was evaluated against maximum-cardinality matching. For `SRR28349244`, all candidate pairings are strictly 1:1 with zero multiplicity, making greedy matching equivalent to the global maximum. The theoretical starvation vulnerability of greedy matching was formally demonstrated via a dedicated unit test in `test_m6_p3c_reconcile.py`.

---

## 2. Sample Metadata & Resource Allocation

### A. Run Identity & Metadata
*Note on Biological Identity: Per project standards, `Cas9_1` is described strictly by its submitted sample name (`Cas9_1`). It is not characterized as a proven clonal isolate unless independently verified in wet-lab documentation.*

| Attribute | Specification | Verification / Provenance |
| :--- | :--- | :--- |
| **Run Accession** | `SRR28349244` | Verified in NCBI SRA / PRJNA1088182 |
| **Sample Name** | `Cas9_1` | Re-sequenced engineered *E. coli* pool |
| **BioSample / Experiment** | `SAMN40463616` / `SRX23955191` | Matched with `run_provenance.tsv` |
| **Organism / Strain** | *Escherichia coli* K-12 MG1655 | Reference strain |
| **Reference Contig** | `NC_000913.3` (4,641,652 bp) | File SHA256: `313752bfe...` |
| **Reference Sequence SHA256**| `e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0` | Gate verified |
| **FASTQ R1 Checksum** | `df654a7884e4daad915afdd4042b3dab29090226acbabbceed0028c9b9c7890f` | Matched with staged provenance |
| **FASTQ R2 Checksum** | `934e6b76160e60e145fed5ed7eb480f10c1857cdf4212d9a367f2574b120741d` | Matched with staged provenance |
| **Total Read Spots** | 1,600,722 paired spots (3,201,444 reads) | Paired Illumina HiSeq |

### B. Slurm Job Execution & Telemetry Comparison
Both tools were benchmarked on Slurm partition `qcpu_18i` using 8 CPU cores and 32 GB RAM allocation.

| Metric | breseq 0.40.2 (Oracle) | ProkaDiff Evidence (P1b Baseline) | Ratio / Difference |
| :--- | :---: | :---: | :---: |
| **Slurm Job ID** | `2997962` | `2894180` | Single-node runs |
| **Host Node** | `bnode2` | `bnode1` | Slurm partition `qcpu_18i` |
| **Exit Status** | `0` (Success) | `0` (Success) | Clean execution |
| **Wall Clock Time** | **12m 48.46s** (768.5 s) | **1h 35m 20s** (5,720 s) | breseq is $7.44\times$ faster |
| **User CPU Time** | 1,566.26 s | 45,084.36 s | ProkaDiff consumed $28.8\times$ CPU |
| **System CPU Time** | 137.23 s | 72.17 s | breseq had more I/O context |
| **CPU Utilization** | 221% | 789% | ProkaDiff parallelized across 8 cores |
| **Peak Memory (MaxRSS)** | **1,548,752 KB** (~1.55 GB) | **2,918,432 KB** (~2.92 GB) | breseq used $47\%$ less peak memory |

*Telemetry Insight:* ProkaDiff's P1b baseline spends significant time on exhaustively sorting and tracking multi-overlapping split reads across every single candidate window. breseq applies an initial hash-based candidate junction filter, reducing the downstream graph construction overhead.

---

## 3. Schema Validation & Divergence Findings

Analysis of breseq's raw Genome Diff (`output.gd`) and directory structure revealed several key differences between the callers:

1. **Contig Accession Normalization:**
   - **breseq Behavior:** breseq parses GenBank annotation headers using the LOCUS field (`NC_000913`) rather than the versioned accession (`NC_000913.3`).
   - **ProkaDiff Behavior:** ProkaDiff uses the FASTA record header (`NC_000913.3`).
   - **Resolution:** `reconcile_m6_p3c.py` implements accession-version normalization via `contig_match()`, ensuring robust contig identity across both naming conventions.
2. **Junction Evidence Attributes:**
   - **breseq Format:** breseq records junction depth in `new_junction_read_count` (total) and `coverage_plus` / `coverage_minus` (strand counts).
   - **ProkaDiff Format:** ProkaDiff uses `pd_support_reads`, `pd_plus_reads`, and `pd_minus_reads`.
   - **Resolution:** `reconcile_m6_p3c.py` property getters support both conventions transparently.
3. **Internal Filtering and Reject Flags:**
   - breseq retains non-passing junction candidates in `output.gd` marked with `reject=FREQUENCY_CUTOFF` or `reject=COVERAGE_EVENNESS_SKEW`. ProkaDiff omits sub-threshold candidates during the initial graph filtering phase.

---

## 4. Variant Call Concordance

### A. Small Variants (SNPs and Indels)
All 9 small variants detected by ProkaDiff matched breseq predictions with **100% coordinate and base concordance**:

| ProkaDiff Record | Variant Type | Position | Ref | Alt | breseq Record | Concordance Status |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `SNP 1` | SNP | 223,648 | A | G | `SNP 1` | **EXACT MATCH** |
| `SNP 2` | SNP | 297,931 | C | T | `SNP 2` | **EXACT MATCH** |
| `SNP 3` | SNP | 301,690 | C | T | `SNP 3` | **EXACT MATCH** |
| `INS 4` | INS | 2,018,530 | . | C | `INS 6` | **EXACT MATCH** |
| `SNP 5` | SNP | 2,177,839 | T | C | `SNP 7` | **EXACT MATCH** |
| `SNP 6` | SNP | 2,177,842 | A | G | `SNP 8` | **EXACT MATCH** |
| `SNP 7` | SNP | 2,456,631 | A | G | `SNP 9` | **EXACT MATCH** |
| `SNP 8` | SNP | 3,911,857 | T | C | `SNP 12` | **EXACT MATCH** |
| `INS 9` | INS | 4,296,381 | . | GC | `INS 13` | **EXACT MATCH** |

### B. Structural Mutations (Deletions & Mobile Element Insertions)
High-level structural mutations promoted from junction and missing coverage evidence were evaluated under the decoupled 4-attribute MOB concordance model (manifest: [`prjna1088182_m6_p3c_mob_concordance.tsv`](manifests/prjna1088182_m6_p3c_mob_concordance.tsv)):

| Event Type | ProkaDiff ID | breseq ID | Position | Element / Size | Strand | TSD / Dup | Concordance Classification | Notes |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- | :--- |
| **Large Deletion** | `DEL 73` | `DEL 11` | 3,815,859 | 82 bp DEL | . | . | **EXACT MATCH** | Identical 82 bp span (3,815,859..3,815,940) |
| **IS Insertion** | `MOB 148` | `MOB 4` | 889,647 | IS5 | -1 | 4 bp vs 4 bp | **EXACT_MOB_MATCH** | Exact site, exact family (IS5), exact strand, exact TSD |
| **IS Insertion** | `MOB 147` | `MOB 5` | 1,873,031 | IS1 | -1 | 8 bp vs 9 bp | **CONCORDANT_MOB_1BP_TSD_RESOLVED** | Exact site & family (IS1); 1 bp microhomology TSD slip |
| **IS Insertion** | `MOB 146` | `MOB 10` | 3,799,683 | IS1 | -1 | 8 bp vs 8 bp | **EXACT_MOB_MATCH** | Exact site, exact family (IS1), exact strand, exact TSD |

*Scientific Note on MOB 147 (1 bp TSD Microhomology Slip):*
Both callers place the insertion anchor at coordinate 1,873,031 and identify the inserted element as IS1 on the negative strand. However, breseq reports a 9 bp target site duplication (`dup=9`, constituent breakend at 1,873,039), whereas ProkaDiff reports an 8 bp target site duplication (`dup=8`, constituent breakend at 1,873,038). This 1 bp difference arises from aligner resolution of split-read overhangs across the microhomologous boundary (`|dup_A - dup_B| = 1 bp`). In strict accordance with the decoupled concordance model, this event is classified as `CONCORDANT_MOB_1BP_TSD_RESOLVED` rather than an exact match.

---

## 5. Junction Concordance & 4D Evidence Classification

### A. Summary of ProkaDiff Candidate Junctions ($N = 49$)

```
+--------------------------------------------------------------------------+
| Total ProkaDiff Candidate Junctions: 49                                 |
+--------------------------------------------------------------------------+
  |
  +--> Dimension 1: Genomic Placement
  |      - UNIQUE_ANCHOR:   38 (77.6%)
  |      - REPEAT_ANCHORED:  6 (12.2%)
  |      - REPEAT_AMBIGUOUS: 5 (10.2%)
  |
  +--> Dimension 2: Read Evidence
  |      - STRONG:   5 (10.2%)
  |      - LIMITED: 44 (89.8%)
  |
  +--> Dimension 3: Computational Corroboration
  |      - CORROBORATED (Exact Coordinate Match, D=0):           5 (10.2%)
  |      - CORROBORATED_AT_FAMILY_LEVEL_UNRESOLVED_AT_LOCUS:      6 (12.2%)
  |      - NOT_CORROBORATED (breseq uncalled / near-miss):       38 (77.6%)
  |
  +--> Dimension 4: Biological Validation
  |      - NOT_EVALUABLE: 49 (100.0%) [Requires orthogonal non-sequencing assay]
  |
  +--> Composite Classification
         - CANDIDATE: 49 (100.0%)
```

### B. Corroborated Structural Records ($D=0$)

| ProkaDiff ID | Breakend 1 | Breakend 2 | breseq ID | Match Type | Support Reads | Event Context |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **`97`** | `1 (+)` | `4,641,652 (-)` | `33` | `EXACT_MATCH` ($D=0$) | 70 | Circular genome origin-terminus junction |
| **`124`** | `1,207,805 (-)` | `1,209,602 (-)` | `40` | `EXACT_MATCH` ($D=0$) | 114 | Intra-chromosomal structural junction |
| **`134`** | `2,357,873 (-)` | `2,357,886 (+)` | `44` | `EXACT_MATCH` ($D=0$) | 16 | Small structural rearrangement candidate |
| **`140`** | `3,815,858 (-)` | `3,815,941 (+)` | `47` | `EXACT_MATCH` ($D=0$) | 105 | Constituent JC for 82 bp `DEL 73` / `DEL 11` |
| **`143`** | `4,171,803 (-)` | `4,213,162 (-)` | `48` | `EXACT_MATCH` ($D=0$) | 2,878 | High-depth chromosomal structural junction |

### C. Repeat-Family Corroboration (Section 3.J Case 6)
Six ProkaDiff junctions were corroborated at the repeat-family level (`CORROBORATED_AT_FAMILY_LEVEL_UNRESOLVED_AT_LOCUS`). These junctions represent the dual breakends of the 3 promoted mobile element insertions (`MOB 146`, `MOB 147`, `MOB 148`):
- Both callers identified the identical insertion anchor locus on the non-repeat flank.
- Because *E. coli* MG1655 contains multiple identical genomic copies of IS1 and IS5, the callers mapped split reads to different paralogous donor copies.
- Under Section 3.J Case 6, these events are recognized as **biologically equivalent structural insertions** despite divergence in internal repeat coordinates.

### D. Full Bidirectional Event-Accounting Specification & Matrix

Every junction called by either tool is fully and non-redundantly accounted for in the bidirectional matrix (manifest: [`prjna1088182_m6_p3c_pilot_bidirectional_jc_matrix.tsv`](manifests/prjna1088182_m6_p3c_pilot_bidirectional_jc_matrix.tsv)). The 17 breseq junctions and 49 ProkaDiff junctions partition into 5 mutually exclusive categories ($N = 54$ total rows):

```
+---------------------------------------------------------------------------------------+
| Bidirectional Bijective Accounting (54 Matrix Rows)                                   |
+---------------------------------------------------------------------------------------+
| Category               | ProkaDiff JCs (49) | breseq JCs (17) | Accounting Status     |
+------------------------+--------------------+-----------------+-----------------------+
| EXACT_MATCH            | 5 (97,124,134,     | 5 (33,40,44,    | Mutual Exact Match    |
|                        |    140,143)        |    47,48)       | (D = 0 bp)            |
+------------------------+--------------------+-----------------+-----------------------+
| REPEAT_FAMILY_MATCH    | 6 (98,100,106,     | 6 (45,43,42,    | Case 6 Paralogous     |
|                        |    120,121,132)    |    35,38,34)    | MOB constituent edges |
+------------------------+--------------------+-----------------+-----------------------+
| NEAR_MISS_BOUNDARY     | 1 (123)            | 1 (39)          | 16 bp microhomology   |
|                        |                    |                 | window slip (D=10 bp) |
+------------------------+--------------------+-----------------+-----------------------+
| BRESEQ_ONLY_REJECTED   | 0                  | 5 (36,37,41,    | Internal breseq skew/ |
|                        |                    |    46,49)       | frequency reject tags |
+------------------------+--------------------+-----------------+-----------------------+
| PROKADIFF_ONLY         | 37                 | 0               | 4 multi-copy chimeras |
|                        |                    |                 | + 33 sub-threshold    |
+------------------------+--------------------+-----------------+-----------------------+
| Total Unique Accounts  | 49 / 49 (100.0%)   | 17 / 17 (100.0%)| 0 Duplicate / Omitted |
+---------------------------------------------------------------------------------------+
```

#### 1. Near-Miss Boundary Analysis (ProkaDiff 123 <-> breseq 39)
- **Coordinates:** ProkaDiff JC 123 (`1207784(+) ~ 1209623(+)`, support=82) vs breseq JC 39 (`1207790(+) ~ 1209619(+)`, support=36).
- **Distance:** Breakend 1 delta = 6 bp, Breakend 2 delta = 4 bp ($D = 10$ bp).
- **Algorithmic Provenance:** breseq reports `alignment_overlap=16`, establishing that a 16 bp exact sequence microhomology exists across this breakpoint. Within a 16 bp homology window, different alignment algorithms place the boundary at different offsets. ProkaDiff aligned the split reads with a 6 bp offset relative to breseq. Because the contract tolerance requires $\le 5$ bp on *each* individual breakend, breakend 1 ($\Delta = 6$ bp) narrowly exceeds the threshold, accurately categorizing this event as a near-miss microhomology boundary artifact rather than a true topological discordance.

#### 2. breseq-Only Sub-Threshold Candidates ($N = 5$)
All 5 breseq-only junctions (JCs 36, 37, 41, 46, 49) represent low-frequency split reads ($N = 3$ to $5$ supporting reads) that failed breseq's internal quality filters:
- JCs 37, 46, 49 bear `reject=COVERAGE_EVENNESS_SKEW,FREQUENCY_CUTOFF`.
- JCs 36, 41 bear `reject=COVERAGE_EVENNESS_SKEW`.
These candidates are recorded in breseq's `output.gd` as rejected evidence items, and were filtered out by ProkaDiff's minimum support depth thresholds.

#### 3. ProkaDiff-Only Candidates ($N = 37$)
The 37 ProkaDiff-only junction candidates divide into two distinct algorithmic classes:
1. **Multi-Copy Repeat Chimeras ($N = 4$):**
   - JC 103 (14,956 reads) and JC 136 (6,489 reads) connecting 16S/23S ribosomal RNA operons (`rrsB` / `rrsE` / `rrnD`).
   - JC 139 (13,622 reads) internal to the `rrnE` operon.
   - JC 99 (590 reads) between IS1 repeats.
   These junctions possess massive read support but connect identical multi-copy repeats across the genome. breseq actively suppresses these using its coverage distribution calibration and repeat redundancy filters (`side_redundant=1`).
2. **Sub-Threshold / Marginal Candidates ($N = 33$):**
   Low to moderate read support (depth 4 to 48 reads) that fall below breseq's junction frequency cutoffs.

---

## 6. Greedy vs. Maximum-Cardinality Matching Assessment

### A. Assessment on Pilot Dataset SRR28349244
An empirical check for candidate multiplicity was performed across all $49 \times 17$ potential junction pairs:
- **Competing Candidate Multiplicity:** Exactly **0** junctions exhibited multiple valid comparator matches within tolerance $\le 5$ bp.
- **Ordering Sensitivity:** Every corroborated match in `SRR28349244` had an exact distance of $D=0$ and was strictly 1:1.
- **Conclusion:** On `SRR28349244`, greedy matching produces the identical assignment as maximum-cardinality and minimum-weight matching.

### B. Formal Demonstration of Greedy Matching Starvation
To evaluate edge-case behavior on more complex cohorts, a synthetic unit test was added to [`test_m6_p3c_reconcile.py`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/scripts/test_m6_p3c_reconcile.py):
```python
def test_greedy_vs_maximum_cardinality_counterexample():
    # JA1 at 1000, JA2 at 1004; JB1 at 1001, JB2 at 997.
    # Distances: (JA1, JB1)=1, (JA1, JB2)=3, (JA2, JB1)=3, (JA2, JB2)=7 (invalid)
    # Greedy takes (JA1, JB1) first (min D=1), starving JA2 -> Cardinality = 1.
    # Maximum-cardinality pairs (JA1, JB2) and (JA2, JB1) -> Cardinality = 2.
```
- **Finding:** Greedy bipartite matching can yield a strictly smaller match cardinality ($|M_{greedy}| = 1 < |M_{max}| = 2$) when competing candidates cluster near the 5 bp microhomology boundary.
- **Recommendation:** While greedy matching is fully sufficient for `SRR28349244`, maximum-cardinality matching (e.g., Hopcroft-Karp) should be considered if dense candidate clusters are observed in high-complexity cohorts like `Cas9_110_111` (`SRR28349173`).

---

## 7. Limitations & Scientific Integrity Declarations

1. **Illumina Resequencing Chimeras:**
   Concordance between ProkaDiff and breseq demonstrates that both algorithms reliably detect chimeric reads present in the sequencing library. It does **not** prove that the junction exists in the physical bacterial chromosome. In-vitro chimeras formed during PCR amplification will be corroborated by both computational pipelines.
2. **Dimension 4 Invariance:**
   Because BioProject PRJNA1088182 lacks orthogonal PCR breakpoint Sanger sequencing or long-read Oxford Nanopore corroboration, all structural calls are strictly classified as `NOT_EVALUABLE` in Dimension 4 (Biological Validation).
3. **Storage & Cohort Hygiene:**
   All pilot outputs are safely isolated under `benchmark/real_world/results/prjna1088182/m6_p3c/SRR28349244/`. The baseline registry (`VALIDATION_REGISTRY.tsv`) remains strictly `PENDING`, and historical P2a/P2b/P3a manifests remain completely untouched.
