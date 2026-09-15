# ProkaDiff Parity Discrepancy Registry

This registry tracks discrepancies identified between ProkaDiff and breseq oracle outputs.

---

## CASE-001: Ambiguous N base mapped to A in RA pileup

- **dataset**: `ambiguous_n` (`testdata/parity/ambiguous_base`)
- **mutation**: N-containing reads inappropriately called as SNP A
- **breseq**: N is filtered or ignored as ambiguous evidence (no mutation called)
- **prokadiff (before)**: `base_index` fallback `_ => 0` counted N as A
- **category**: RA / pileup
- **status**: fixed
- **minimal_fixture**: `testdata/parity/ambiguous_base/`
- **fix_summary**: Updated `base_index(b: u8) -> Option<usize>` in `ra.rs` to ignore non-ACGT bases in allele counting.
- **parity_after_fix**: `strict_over = 0`, `strict_under = 0`

---

## CASE-002: Adjacent substitutions emitted as separate SNPs instead of SUB

- **dataset**: `sub_adjacent` (`testdata/parity/synthetic_sub`)
- **mutation**: Multi-base consecutive substitutions (e.g., AC -> TT)
- **breseq**: Merges adjacent 1-bp changes into single `SUB` mutation (`SUB <id> <parents> <seq_id> <pos> <size> <new_seq>`)
- **prokadiff (before)**: Emitted individual `SNP` entries
- **category**: normalization / mutation promotion / GenomeDiff
- **status**: fixed
- **minimal_fixture**: `testdata/parity/synthetic_sub/`
- **fix_summary**: Added `GdKind::Sub` in `prokadiff-gd` and merged adjacent consensus SNPs into `SUB` records in `emit.rs`.
- **parity_after_fix**: `strict_over = 0`, `strict_under = 0` (1 SUB vs 1 SUB)

---

## CASE-003: Real-data Clonal JC calling differences

- **dataset**: `clonal` (`REL606` / `SRR030257`)
- **mutation**: Structural junctions (`JC`)
- **breseq**: 34 accepted JCs, uses coverage-evenness skew test and read-position hash filtering
- **prokadiff**: 34 JCs emitted; 5 match breseq within $\pm 5$ bp; 29 over-called (unmatched) due to lack of skew filter in Stage 1
- **category**: JC
- **status**: open / documented (Stage 2/3 roadmap for coverage-evenness skew filtering)
- **minimal_fixture**: `benchmark/results/clonal_2433013/`
- **notes**: Red-line mutations (SNP, INS, DEL, MOB) match with 0 over-calls; JC differences isolated to unpromoted junctions.

---

## CASE-004: Mobile element insertion breakpoint coordinate jitter

- **dataset**: `synth_is_mob`
- **mutation**: `MOB` / `JC`
- **breseq**: Emits MOB with specific target-site duplication coordinate
- **prokadiff**: Emits MOB within $\pm 5$ bp clustering tolerance
- **category**: MOB / normalization
- **status**: tolerant-match
- **minimal_fixture**: `benchmark/results/synth_is_mob_2432998/`
- **parity_after_fix**: Strict mismatch (different representative coordinate), matched under Tier 3 tolerant matching ($\pm 5$ bp).
