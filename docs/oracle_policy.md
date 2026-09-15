# Oracle Policy for ProkaDiff

## 1. Scientific Positioning and Calling Oracle

ProkaDiff v0.x does not attempt to redefine bacterial short-read variant calling.

```text
breseq (0.40.2)
  │
  │ primary behavioral oracle
  ▼
ProkaDiff Rust Engine
  │
  ├── Reproduce breseq mutation calling behavior
  ├── Substantially accelerate execution via native Rust
  ├── Retain full GenomeDiff (.gd) format compatibility
  └── Add starter-vs-edited subtraction and CRISPR off-target interpretation
```

**Core Principle:**
- **Rust implementation != new calling standard.**
- **Rust implementation = accelerated, compatible implementation.**
- **breseq parity > algorithmic novelty.**

Differences between ProkaDiff and breseq are treated as unresolved discrepancies unless independently verified with extensive biological ground-truth data or experimental validation. We do not assume breseq is wrong when outputs diverge.

---

## 2. Pinned Oracle Environment

Variant calling sensitivity and junction resolution are intimately coupled to aligner heuristics. As documented in the breseq publication and release notes, differences in Bowtie2 versions can directly alter read mappings and downstream variant calls.

Therefore, the oracle environment is strictly pinned:
- **breseq**: `0.40.2`
- **bowtie2**: `2.5.4`
- **gdtools**: `0.40.2`

Oracle tools and versions are recorded in `testdata/VERSIONS.txt` and `benchmark/oracle/versions.tsv`. Upgrades to the oracle toolchain require explicit side-by-side benchmark verification before adoption.

---

## 3. The Three Parity Tiers

Parity evaluation between ProkaDiff and breseq must never collapse into a single binary pass/fail or count metric. Every benchmark must distinguish:

### Tier 1: Strict Parity
- Exact coordinate, mutation type (`SNP`, `SUB`, `INS`, `DEL`, `MOB`, `AMP`), and allele identity.
- Evaluated bidirectionally using `gdtools SUBTRACT`:
  - `gdtools SUBTRACT prokadiff.gd breseq.gd` → Over-calls (`over_call.gd`)
  - `gdtools SUBTRACT breseq.gd prokadiff.gd` → Under-calls (`under_call.gd`)
- This is the primary gate for point mutations, small indels, and standard mobilome calls.

### Tier 2: Normalized Parity
- 3' homopolymer and tandem-repeat sliding normalization (e.g., indel normalization per 1-based rightmost convention).
- Evaluated after canonical GenomeDiff left/right alignment rules are applied.

### Tier 3: Diagnostic Tolerant / Biological-Event Parity
- For split-read junctions (`JC`), mobilome breakpoints (`MOB`), or complex repeat events:
  - Evaluation allows coordinate tolerance (e.g., $\pm 5$ bp jitter from alternative read placement or repeat copy assignment).
- **Rule:** Diagnostic tolerance is strictly for debugging and root-cause analysis. A mismatch under Tier 1 must **never** be silently converted to a `PASS` merely because it matches under tolerant comparison. Both metrics must be reported simultaneously.

---

## 4. Discrepancy Handling Protocol

When ProkaDiff output differs from breseq:
1. **Preserve Inputs & Outputs:** Save input FASTQ/BAM, reference FASTA/GenBank, `breseq.gd`, and `prokadiff.gd`.
2. **Isolate Cause:** Determine whether discrepancy stems from RA (base quality/pileup/allele fraction), MC (coverage gap extension), JC (split seed/breakpoint calculation), or normalization/filtering.
3. **Construct Minimal Reproducible Fixture:** Create a minimal synthetic or extracted testcase in `testdata/parity/regressions/<issue_name>/`.
4. **Implement Conservative Fix:** Implement the minimal change required in ProkaDiff to align with the oracle.
5. **Add Regression Tests:** Add layer-0 unit tests and integration tests covering the fixture.
6. **Register in Discrepancy Registry:** Document the case in `benchmark/PARITY_ISSUES.md`.
7. **Verify Against Benchmark Matrix:** Re-run the full parity suite to ensure no regressions were introduced.

---

## 5. Prohibition of Heuristic Tweaking

Bioinformatic thresholds must not be adjusted based on intuition or subjective preference:
- Coverage thresholds (minimum depth, coverage-gap seeds)
- Allele frequency / consensus thresholds
- MAPQ and base quality cutoffs
- Junction split-read support and breakpoint overlap criteria

Any modification to calling thresholds must be backed by breseq oracle data or documented regression cases.

---

## 6. Scope Boundaries: Calling vs. Downstream Editing QC

ProkaDiff separates its scientific responsibilities into clear layers:
1. **Layer 1 (Variant Detection Engine):** breseq-compatible Bowtie2 alignment, RA/MC/JC evidence calling, GenomeDiff emission.
2. **Layer 2 (Starter-Edited Subtraction):** Direct subtraction of starter-strain background mutations ($M_{\text{edited}} - M_{\text{starter}} - M_{\text{intended}}$).
3. **Layer 3 (Gene-Editing Interpretation):**
   - Spatial association between called mutations and predicted guide target/homolog sites (`offtarget_sites.tsv`, `mutation_offtarget_links.tsv`).
   - Standard reference tools (Cas-OFFinder, FlashFry, CRISPRitz) serve as oracles for off-target search and scoring models.
   - Predictions are explicitly designated as *computational candidate sites / spatial associations*, never claimed as experimentally confirmed cleavage events without biological validation.
