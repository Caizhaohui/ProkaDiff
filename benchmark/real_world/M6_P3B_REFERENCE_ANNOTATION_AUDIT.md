# M6-P3b Reference Annotation & MOB Evidence Baseline Audit

**Milestone:** `M6-P3b` (Reference Annotation & MOB Validation)
**Task ID:** `P3B_BASELINE_AUDIT`
**Base Commit:** `41dec1f2a3a5f4f0df4000cb080447fa414ce5dc`
**Date:** 2026-10-21
**Author:** Antigravity (Project Orchestrator)
**Status:** **AUDIT COMPLETE — READY FOR CLAUDE SCIENTIFIC DESIGN REVIEW**

---

## 1. Executive Summary

Phase M6-P2a evaluated 25 eligible paired WGS runs from PRJNA1088182 against reference genomes according to strain:
- **10 MG1655 runs** against `NC_000913.3` using GenBank input `data/references/MG1655.gbk`. Result: **MOB = 2 or 3** per run (22 MOB calls total).
- **15 BW25113 runs** against `CP009273.1` using plain FASTA input `data/references/m6_identity/CP009273.1.fasta`. Result: **MOB = 0** across all 15 runs, despite 35 to 72 JC calls per run (749 JC calls total).

This audit investigates the computational and algorithmic basis of this difference: whether the zero-MOB outcome in BW25113 is an artifact of reference file format (lack of repeat annotations in FASTA) or reflects genuine absence of mobile element transpositions in the sequencing reads.

---

## 2. Repository Baseline & Binary Identity

1. **Working Tree & Commit State:**
   - Active branch: `main`
   - Latest commit: [`41dec1f`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff) (`benchmark: validate M6 P3a complete read-level diagnostic set`)
   - Uncommitted working tree files: `stage_m6_p1a_fastq.sh` (mode support) and untracked P2a execution logs/manifests. Core engine (`crates/**`) and contracts are 100% clean.
2. **Production Binary Identity:**
   - Path: `target/release/prokadiff`
   - SHA256: `02077c761b1c9bd83cf1b626a3e2e955d8d2a1c63abb3e35852421e0dbc20bfd`
   - Matches the exact binary SHA256 recorded in [`prjna1088182_m6_p2a_run_results.tsv`](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/benchmark/real_world/manifests/prjna1088182_m6_p2a_run_results.tsv).

---

## 3. Reference Genome File Formats & Sequence Hashes

| Strain | Accession | File Format | Local Path | Sequence Length | Nucleotide SHA256 | Repeat Annotations Parsed |
| :--- | :---: | :---: | :--- | ---: | :--- | :---: |
| **MG1655** | `NC_000913.3` | GenBank | `data/references/MG1655.gbk` | 4,641,652 bp | `e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0` | **50** |
| **BW25113 (P2a)** | `CP009273.1` | FASTA | `data/references/m6_identity/CP009273.1.fasta` | 4,631,469 bp | `3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce` | **0** |
| **BW25113 (Deposited)** | `CP009273.1` | GenBank | NCBI Entrez `CP009273.1` (Grenier et al., 2014) | 4,631,469 bp | `3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce` | **401** |

*Verification:* Direct byte extraction from the NCBI `CP009273.1` GenBank record confirms **100% exact sequence identity** with `CP009273.1.fasta` (Length: 4,631,469 bp, SHA256: `3761072f...`).

---

## 4. Algorithmic Root Cause Analysis of MOB Promotion

Line-by-line inspection of the production evidence engine identified the exact mechanism:

### 4.1 Repeat Ingestion: `crates/prokadiff/src/main.rs` & `crates/prokadiff-evidence/src/engine/mod.rs`
1. `run_evidence` calls `load_repeats(&args.refs)` ([`main.rs` L156](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/crates/prokadiff/src/main.rs#L156)):
   ```rust
   fn load_repeats(refs: &[PathBuf]) -> Vec<prokadiff_evidence::fasta::RepeatRegion> {
       let mut repeats = Vec::new();
       for r in refs {
           if let Ok(reps) = prokadiff_evidence::fasta::parse_genbank_repeats(r) {
               repeats.extend(reps);
           }
       }
       repeats
   }
   ```
2. In `run_sample` ([`mod.rs` L84-L110](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/crates/prokadiff-evidence/src/engine/mod.rs#L84-L110)), if `opts.repeats` is empty, it attempts `parse_genbank_repeats(ref_fa)`. If `ref_fa` is a FASTA file and no `.gbk` exists in its parent directory, `opts.repeats` remains **empty (`Vec::new()`)**.

### 4.2 Annotation Parsing: `crates/prokadiff-evidence/src/fasta.rs`
1. `parse_genbank_repeats` ([`fasta.rs` L144-L259](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/crates/prokadiff-evidence/src/fasta.rs#L144-L259)) scans for feature lines starting with `repeat_region` or `mobile_element`.
2. For `mobile_element`, qualifier `/mobile_element="insertion sequence:..."` extracts the IS name. In `MG1655.gbk` and `CP009273.1.gbk`, NCBI uses qualifier `/mobile_element_type="..."`, which bypasses the qualifier parser, defaulting `name` to `"repeat"`.
3. In `MG1655.gbk`, 50 features are parsed with name `"repeat"`. In `CP009273.1.gbk`, 401 features are parsed (48 with name `"repeat"`, 353 with specific note names).

### 4.3 MOB Candidate Promotion: `crates/prokadiff-evidence/src/engine/emit.rs`
1. In `emit_from_pileup` ([`emit.rs` L402-L560](file:///hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/crates/prokadiff-evidence/src/engine/emit.rs#L402-L560)):
   ```rust
   for &(c_tgt1, n_tgt1, p_tgt1, m_tgt1, n_rep1, p_rep1, _m_rep1) in &j1_sides {
       let Some(rep1) = find_repeat_at(&opts.repeats, n_rep1, p_rep1) else {
           continue;
       };
       ...
       let Some(rep2) = find_repeat_at(&opts.repeats, n_rep2, p_rep2) else {
           continue;
       };
       if rep1.name != rep2.name {
           continue;
       }
   ```
2. **Deterministic Mathematical Consequence:**
   - When `opts.repeats` is empty (as occurred with `CP009273.1.fasta`), `find_repeat_at` **always returns `None`**.
   - As a result, candidate junction pairing for MOB formation is bypassed entirely.
   - **`candidate_mobs` is guaranteed to be empty**, and no MOB entries can ever be emitted in `output.gd`.
   - All mobile element junctions remain unpromoted as standalone `JC` records.

---

## 5. P2a Cohort Structural Event Summary

In PRJNA1088182 P2a:
- **MG1655 Cohort (10 paired runs):**
  - Average JC calls: 133 (range: 45 to 505)
  - MOB calls: 2 or 3 per run (Total: 22 MOB calls)
  - MOB parent JCs: Successfully absorbed and tagged with `mob_evidence=1`
- **BW25113 Cohort (15 paired runs):**
  - Average JC calls: 50 (range: 35 to 72; Total: 749 JC calls)
  - MOB calls: **0** across all 15 runs
  - MOB parent JCs: 0

---

## 6. Audit Conclusion & Handoff to Claude

1. **Input-Format Discrepancy Confirmed:** The observed `MOB=0` across 15 BW25113 runs is an inevitable mathematical consequence of passing a plain FASTA file to an engine whose MOB promotion logic strictly requires reference repeat annotations.
2. **Nucleotide Identity Confirmed:** An authentic, sequence-identical GenBank record for `CP009273.1` exists in NCBI Entrez (4,631,469 bp, SHA256 `3761072f...`).
3. **Hypothesis Formulation:** Providing `CP009273.1.gbk` provides `opts.repeats` to the engine without altering a single base pair of the reference genome sequence or the Bowtie2 alignment index.
4. **Handoff:** The baseline audit is complete. Claude is requested to review the scientific architecture and design a controlled A/B experiment.
