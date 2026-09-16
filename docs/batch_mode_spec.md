# ProkaDiff Multi-Clone Batch Mode Specification

## 1. Overview & Positioning

In microbial genome engineering workflows, experimentalists typically isolate and screen multiple edited colonies from a single transformation or mutagenesis procedure. 

Current pairwise auditing in ProkaDiff compares:
$$\text{Starter Strain WGS} \longleftrightarrow \text{Edited Strain WGS}$$

Batch mode generalizes this to the one-starter-versus-many-clones paradigm:
```text
                     Starter Strain WGS
                             │
     ┌───────────────┬───────┴───────┬───────────────┐
     ▼               ▼               ▼               ▼
 Edited Clone 01  Edited Clone 02  Edited Clone 03  Edited Clone N
```

The objective of Batch Mode is to:
1. Concurrently or sequentially audit all clones against the same starter background and reference genome.
2. Cross-compare post-edit differential variants across clones to identify **shared** vs **clone-specific** events.
3. Provide an executive **Batch Genome Audit Report** (`batch_report.md`) enabling experimentalists to survey editing completeness, structural integrity, and collateral changes across the entire cohort.

> **Strict Non-Prescriptive Policy (Taskbook §66 & §40):** ProkaDiff provides observational audit metrics, clone-specific variant tallies, and shared event patterns. ProkaDiff **never** automatically declares a clone "best", "safest", or "recommended". Strain selection decisions depend on user-specific phenotypes, expression titers, and downstream biological requirements.

---

## 2. Input Specification: Batch Manifest

Batch runs are configured via a simple TSV manifest file (e.g. `batch_manifest.tsv`), specifying the starter strain, shared reference, and all edited clones.

### Manifest Columns

| Column | Required | Description |
| :--- | :---: | :--- |
| `clone_id` | Yes | Unique alphanumeric identifier for each clone (e.g., `C01`, `Clone_del_A`) |
| `r1` | Yes | Path to R1 FASTQ (`.fq`, `.fq.gz`, `.fastq`, `.fastq.gz`) |
| `r2` | No | Path to R2 FASTQ for paired-end sequencing (leave empty for single-end) |
| `intended_tsv` | No | Optional per-clone intended edit TSV table (if omitted, uses global `--intended`) |
| `spacer` | No | Optional per-clone spacer sequence override |
| `pam` | No | Optional per-clone PAM profile override |

### Example `batch_manifest.tsv`

```tsv
clone_id	r1	r2	intended_tsv	spacer	pam
C01	data/C01_R1.fq.gz	data/C01_R2.fq.gz	intended.tsv	GAGTTCATCTACGCCGTGAA	NGG
C02	data/C02_R1.fq.gz	data/C02_R2.fq.gz	intended.tsv	GAGTTCATCTACGCCGTGAA	NGG
C03	data/C03_R1.fq.gz	data/C03_R2.fq.gz	intended.tsv	GAGTTCATCTACGCCGTGAA	NGG
C04	data/C04_R1.fq.gz	data/C04_R2.fq.gz	intended.tsv	GAGTTCATCTACGCCGTGAA	NGG
```

### CLI Invocations

```bash
# Basic batch execution
prokadiff \
  --starter starter_R1.fq.gz starter_R2.fq.gz \
  --ref reference.gbk \
  --batch-manifest batch_manifest.tsv \
  --editor cas9 \
  --outdir batch_audit/ \
  --threads 16
```

---

## 3. Data Structures & Domain Models

### 3.1 Clone Audit Record

```rust
/// Audit outcome for an individual clone within the batch.
#[derive(Clone, Debug)]
pub struct CloneAuditEntry {
    pub clone_id: String,
    pub audit_result: AuditResult,
    pub output_dir: PathBuf,
}
```

### 3.2 Shared Variant Recurrence

```rust
/// Recurrence classification of a post-edit variant across clones.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RecurrenceStatus {
    /// Variant observed only in a single clone.
    CloneSpecific,
    /// Variant observed in a subset of clones (count < total).
    SharedSubset { count: usize, total: usize },
    /// Variant observed across all evaluated clones.
    SharedAll { total: usize },
}

/// Unified multi-clone variant record.
#[derive(Clone, Debug)]
pub struct BatchVariantRecord {
    pub batch_variant_id: String,
    pub seq_id: String,
    pub position: u64,
    pub end: u64,
    pub gd_type: String,
    pub ref_allele: String,
    pub alt_allele: String,
    pub recurrence: RecurrenceStatus,
    pub observed_clones: Vec<String>,
    pub intended_relation: IntendedRelation,
    pub guide_relation: GuideRelation,
    pub mobile_element: Option<String>,
    pub gene_name: Option<String>,
    pub locus_tag: Option<String>,
}
```

### 3.3 Cohort Aggregate

```rust
/// Cohort-wide batch audit aggregation.
#[derive(Clone, Debug)]
pub struct BatchAuditResult {
    pub starter_names: Vec<String>,
    pub reference_names: Vec<String>,
    pub clones: Vec<CloneAuditEntry>,
    pub shared_variants: Vec<BatchVariantRecord>,
    pub run_timestamp: String,
    pub prokadiff_version: String,
}
```

---

## 4. Output Directory Layout

```text
batch_audit/
├── batch_report.md            # Human-readable cohort executive summary
├── batch_summary.tsv          # Tabular clone-by-clone metric matrix
├── shared_variants.tsv        # Multi-clone recurrence and variant occurrence
├── provenance.tsv             # Technical provenance of the batch execution
└── clones/
    ├── C01/
    │   ├── report.md
    │   ├── edit_outcomes.tsv
    │   ├── post_edit_variants.tsv
    │   ├── summary.txt
    │   ├── unintended.tsv
    │   ├── starter.gd
    │   └── edited.gd
    ├── C02/
    │   └── ...
    └── C03/
        └── ...
```

---

## 5. Machine-Readable Batch Outputs

### 5.1 `batch_summary.tsv`

One line per clone summarizing audit dimensions:

```tsv
clone_id	intended_status	intended_edits_complete	intended_edits_partial	intended_edits_missing	mob_count	structural_count	candidate_offtarget_count	distal_small_count	high_attention_count	review_count	info_count
C01	COMPLETE	1	0	0	0	0	0	2	0	2	1
C02	COMPLETE	1	0	0	1	0	0	1	1	1	1
C03	PARTIAL	0	1	0	0	2	1	4	3	4	0
C04	COMPLETE	1	0	0	0	0	0	0	0	0	1
```

### 5.2 `shared_variants.tsv`

Cross-clone recurrence analysis:

```tsv
batch_var_id	seq_id	position	end	gd_type	ref	alt	recurrence_status	clone_count	total_clones	clones	gene	guide_relation	mobile_element
BVAR_0001	NC_000913.3	1250001	1251200	DEL	1200bp	none	SHARED_ALL	4	4	C01,C02,C03,C04	geneA	ON_TARGET	NONE
BVAR_0002	NC_000913.3	3445210	3445210	MOB	none	IS1	CLONE_SPECIFIC	1	4	C02	intergenic	NONE	IS1(TSD=9bp)
BVAR_0003	NC_000913.3	2135820	2135820	SNP	A	G	CLONE_SPECIFIC	1	4	C03	recA	CANDIDATE_OFF_TARGET	NONE
BVAR_0004	NC_000913.3	4201100	4201100	SNP	C	T	SHARED_SUBSET	3	4	C01,C02,C04	dnaK	NONE	NONE
```

---

## 6. Human-Readable Cohort Report (`batch_report.md`)

### Format Specification

```markdown
# Engineering Batch Genome Audit Report

## 1. Cohort Executive Summary

- **Starter Strain:** `starter_R1.fq.gz`, `starter_R2.fq.gz`
- **Total Evaluated Clones:** 4
- **Intended Edit Completion:**
  - **Complete (all declared edits confirmed):** 3 clones (75.0%)
  - **Partial (aberrant / incomplete integration):** 1 clone (25.0%)
  - **Missing (no intended edit detected):** 0 clones (0.0%)

### Differential Event Distribution Across Clones

| Clone ID | Intended Edit Status | MOB Events | Structural DEL / JC | Candidate Guide OT | Distal Small Variants | Review Priority |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **C01** | `COMPLETE` | 0 | 0 | 0 | 2 | 2 REVIEW |
| **C02** | `COMPLETE` | 1 | 0 | 0 | 1 | 1 HIGH_ATTENTION, 1 REVIEW |
| **C03** | `PARTIAL` | 0 | 2 | 1 | 4 | 3 HIGH_ATTENTION, 4 REVIEW |
| **C04** | `COMPLETE` | 0 | 0 | 0 | 0 | Clean Edit (0 Collateral) |

> **Audit Recommendation:** Clones exhibit heterogeneous collateral profiles. Clone `C04` represents a clean edit with zero observed differential variants outside the intended locus. Clone `C02` carries a novel IS1 transposon insertion requiring manual verification. Clone `C03` failed complete integration and exhibits multiple collateral structural rearrangements.
>
> *ProkaDiff reports observational differences and prioritizes manual review; it does not make biological fitness or safety claims.*

## 2. Recurrent & Shared Post-Edit Variants

Variants observed across multiple clones provide critical clues into parent-clone divergence:
- **Shared across all clones (4/4):** Indicates either an edit occurring at 100% efficiency (on-target) or pre-existing subclonal divergence in the starter culture prior to plating.
- **Shared across subsets (2/4 - 3/4):** Likely early culture drift or common collateral hotspots (e.g. common IS transposition targets or repair hotspots).
- **Clone-specific (1/4):** Independent post-transformation events.

| Variant ID | Locus | Type | Clones | Recurrence | Context / Gene | Review Priority |
| :--- | :--- | :---: | :--- | :---: | :--- | :---: |
| `BVAR_0001` | NC_000913.3:1250001 | DEL | C01, C02, C03, C04 | **4 / 4 (100%)** | Intended locus | INFO |
| `BVAR_0004` | NC_000913.3:4201100 | SNP | C01, C02, C04 | **3 / 4 (75%)** | dnaK (b0014) [CDS] | REVIEW |
| `BVAR_0002` | NC_000913.3:3445210 | MOB | C02 | **1 / 4 (25%)** | IS1 insertion | HIGH_ATTENTION |
| `BVAR_0003` | NC_000913.3:2135820 | SNP | C03 | **1 / 4 (25%)** | Candidate OT (2 mm, AGG) | REVIEW |

## 3. Individual Clone Audits

Detailed audit reports and variant tables for each individual clone are located in `clones/<clone_id>/report.md`.
```

---

## 7. Implementation Roadmap & Milestones

1. **Manifest Parser & Validator:** Implement parsing of `batch_manifest.tsv`, verify file paths exist, check unique `clone_id`.
2. **Batch Pipeline Orchestrator:** Sequential or pool-based execution of single-clone `run_sample` invocations, accumulating `AuditResult` instances.
3. **Cross-Clone Variant Aggregator:** Map coordinates and alleles to assign `batch_variant_id` and compute `RecurrenceStatus`.
4. **Batch Report Generator:** Render `batch_report.md`, `batch_summary.tsv`, and `shared_variants.tsv`.
5. **Fixtures & Testing:** Synthetic multi-clone fixture test verifying complete vs partial clones, shared variant counting, and deterministic markdown generation.
