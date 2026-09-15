# Off-Target Search Coordinate Conventions

## Summary Table

| Tool           | Start base | End semantics | ProkaDiff conversion |
|----------------|------------|---------------|----------------------|
| Cas-OFFinder 2 | 0-based    | inclusive (start + len - 1) | `start = ff_start + 1`, `end = ff_start + len` |
| FlashFry 1.15  | 0-based    | exclusive (half-open)       | `start = ff_start + 1`, `end = ff_stop`         |
| CRISPRitz 2.8.1| TBD        | TBD                         | TBD (BLOCKED_EXTERNAL_DEPENDENCY)               |
| ProkaDiff      | **1-based** | **inclusive**              | canonical internal representation               |

## Cas-OFFinder 2.x Coordinates

Cas-OFFinder 2.x reports:
- `start`: 0-based index of the first nucleotide of the match on the reported strand
- The match length is always `len(spacer) + len(PAM)` (e.g., 20 + 3 = 23 for SpCas9 NGG)

**Conversion to ProkaDiff 1-based inclusive:**
```
prokadiff_start = cas_offinder_start + 1
prokadiff_end   = cas_offinder_start + len(target_seq)
```

**Verified** against controlled synthetic reference (`testdata/offtarget/cas_offinder/reference.fa`)  
Commit: `f46a64c`  12/12 test sites match expected coordinates.

## FlashFry 1.15 Coordinates

FlashFry 1.15 reports:
- `start`: 0-based inclusive (first nt of match)
- `stop`: 0-based exclusive (one past last nt of match, i.e. half-open interval)

**Conversion to ProkaDiff 1-based inclusive:**
```
prokadiff_start = flashfry_start + 1
prokadiff_end   = flashfry_stop      # half-open → 1-based inclusive: last nt at (stop-1), +1 for 1-based = stop
```

**Status**: Convention documented from FlashFry source/documentation.  
**NOT YET CONFIRMED** against real FlashFry 1.15 output with synthetic reference.  
Status: `BLOCKED_EXTERNAL_DEPENDENCY`

When real FlashFry output is available, run:
```bash
python3 benchmark/offtarget/flashfry/scripts/normalize_flashfry.py \
    benchmark/offtarget/flashfry/golden/raw_flashfry.tsv \
    benchmark/offtarget/flashfry/golden/normalized_flashfry.tsv
```
And verify against `testdata/offtarget/cas_offinder/reference.fa` known sites.

## CRISPRitz 2.8.1 Coordinates

Status: `BLOCKED_EXTERNAL_DEPENDENCY`  
Coordinate convention will be determined after real CRISPRitz output is collected.

## Implementation

Coordinate conversions are implemented in:
- `crates/prokadiff-offtarget/src/oracle/cas_offinder_v2.rs` — Cas-OFFinder v2
- `crates/prokadiff-offtarget/src/oracle/flashfry.rs` — FlashFry 1.15

All conversions output `OffTargetSite { start, end }` in **1-based inclusive** format.
