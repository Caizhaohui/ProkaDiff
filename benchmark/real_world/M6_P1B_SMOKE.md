# M6-P1b two-run smoke

Status: **M6-P1b PASS — ready for M6-P2 planning**. M6 remains **ACTIVE**.
M6-P2 has not started. Tier B remains **NOT_EVALUABLE**.

## Scope

Reference-relative publication-event recovery only. Not a parent audit, not
Cas9 causality, and not a specificity or F1 estimate.

Executed paired runs:

| Sample | Run | Publication truth | Reference |
| --- | --- | --- | --- |
| Cas9_1 | SRR28349244 | `301690 C>T`, `COORDINATE_VERIFIED` | `NC_000913.3` |
| Cas9_2 | SRR28349243 | `NONE` | `NC_000913.3` |

`SRR28357593` and `SRR28357592` were not executed. They remain
`SINGLE_BIOLOGICAL_READ_PER_SPOT` and are not part of this paired smoke.

Provenance commits: M6-P0 `06f6259`, P1a read-structure `7025051` (this job's
`HEAD`). Workflow: production `prokadiff evidence` against
`data/references/MG1655.gbk`. No `--starter`. No breseq.

## Gates before calling

Local chromosome length 4,641,652, sequence SHA256
`e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0`.
Both runs resolved to `NC_000913.3`. Compressed and uncompressed FASTQ
SHA256 values matched the P1a provenance, and R1/R2 counts were 1,600,722
and 1,514,507.

## Execution

One Slurm job, `2897505`, partition `qcpu_18i`, 8 CPU, 24G requested, node
`bnode1`, exit 0, elapsed 4:42:34. Bowtie2 2.5.4. Three evidence runs were
sequential, which is why the job was long: each sample spent about 1.5 hours
in sensitive unmatched-read realignment and the junction second pass.

| Run | Wall | Max RSS (`/usr/bin/time`) | Output |
| --- | --- | --- | --- |
| SRR28349244 | 1:35:20 | 2,918,432K | `results/.../m6_p1b/SRR28349244/output.gd` |
| SRR28349243 | 1:33:46 | 2,852,576K | `results/.../m6_p1b/SRR28349243/output.gd` |
| SRR28349244 rerun | 1:32:50 | 2,918,640K | `results/.../m6_p1b/SRR28349244_rerun/output.gd` |

The job-level sacct MaxRSS was 2,920,080K. Per-run figures above are the
ones used here.

## Cas9_1 recovery

Publication evidence and the ProkaDiff call are separate.

Publication: `301690 C>T` on the `U00096.3` / `NC_000913.3` chromosome.
The local reference base at 1-based 301690 is `C`.

ProkaDiff evidence GD:

`SNP 3 . NC_000913.3 301690 T pd_ra_depth=86 pd_ra_frequency=1 pd_ra_support_reads=86`

Classification: **RECOVERED_EXACT**. The numeric `3` is a GD record id, not a
DV1 EventId. This path does not create DifferentialEvents, so EventId checks
are **NOT_EVALUABLE**. No caller threshold was changed.

## Cas9_2 inventory

`NONE` means the publication table reports no mutation. It is not a
false-positive set. `301690` is absent. Measured GD records, excluding MC
from the mutation total: SNP 6, SUB 0, INS 1, DEL 1, JC 48, MOB 2, plus 88
MC evidence records. Label: **ADDITIONAL_CALL_INVENTORY**.

## Products

`output.gd` and `repeat_ambiguous_junctions.tsv` were generated and are the
evaluable reference-relative evidence. The seven M4 differential audit files
were not produced. They are **NOT_EVALUABLE_BY_TIER_A** because they need a
comparator or an intended-edit table. That is not a failure. M2 status was
not emitted.

## Determinism

The Cas9_1 rerun `output.gd` is byte-identical to the first run, including
the `301690 T` record. BAM files were not compared.

## Limitations

No biological parent. No breseq 0.40.2 oracle in this smoke. The publication
remains `PUBLICATION_REPORTED`. Two single-read SRA runs were not executed.
M6 is not closed.
