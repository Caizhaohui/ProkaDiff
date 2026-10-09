# M6-P3a targeted read-level diagnostic plan

Status: **M6-P3a PLAN READY FOR EXECUTION**. M6 remains **ACTIVE**.
This document does not execute alignment. No Slurm job is submitted here.
P2a and P2b tables are unchanged.

The target list is a pointer into the P2b diagnostic table. It is not a new
publication-truth file.

## Frozen P2b state

`manifests/prjna1088182_m6_p2b_event_diagnostics.tsv` has ten
`UNRESOLVED` rows, each once. Each still has:

- `p2a_recovery=NOT_RECOVERED`
- `p2a_diagnostic=UNKNOWN`
- `diagnosis=UNRESOLVED`
- `bam_artifact=ABSENT`
- depth, ALT count, and allele fraction `NOT_MEASURED`

`NOT_MEASURED` is not zero. `3112908 +AAT` stays `NOT_EVALUABLE` and is not
a diagnostic target.

## Initial runs

Checked against the P2b table, the relationship execution reference, and the
P2a FASTQ provenance. Sample names were not taken from this plan's prompt alone.

| Priority | Run | Sample | Reference | Allele | Spots | Gzip bytes R1+R2 |
| --- | --- | --- | --- | --- | ---: | ---: |
| 1 | `SRR28349204` | Cas9_3 _and2others | `NC_000913.3` | `1571257 G>T` | 4,051,149 | 304,839,611 + 320,910,509 |
| 1 | `SRR28349176` | Cas9_104_105 | `NC_000913.3` | `2776013 T>C` | 3,433,823 | 275,848,964 + 290,135,662 |
| 1 | `SRR28349176` | Cas9_104_105 | `NC_000913.3` | `4035824 T>C` | same run | same files |
| 1 | `SRR28349230` | I-SceI_7_8_9 | `CP009273.1` | `2872949 C>A` | 3,230,652 | 248,344,601 + 263,641,290 |

The other six unrecovered alleles stay in
`manifests/prjna1088182_m6_p3a_targets.tsv` at priority 2. They are not in
the first alignment batch. All ten are documented mixed gDNA. Mixture
percentages are not allele fractions and are not a classification input.

Local references, already hashed in P2a:

- MG1655 execution file `data/references/MG1655.gbk`, chromosome SHA256
  `e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0`,
  length 4,641,652.
- BW25113 execution file `data/references/m6_identity/CP009273.1.fasta`,
  chromosome SHA256
  `3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce`,
  length 4,631,469. Header id is `CP009273.1`.

## Alignment to reproduce

RA pileup reads the primary BAM. In `align_to_bam`, that BAM is the
coordinate-sorted conversion of stage-1 SAM only. Stage-2 unmatched
realignment writes separate SAM files and is not merged into that BAM.
The junction pass is a different `AlignKind`. A diagnostic that merged
stage 2, or that ran the junction seed settings, would be a new alignment
policy. The first diagnostic does not do that.

Production stage 1, Bowtie2 2.5.4, from `bowtie2_align_args(PrimaryStage1)`:

```text
bowtie2-build --threads 8 <reference.fasta> <work>/bt2ref
bowtie2 --local --no-unal -p 8 -x <work>/bt2ref \
  -1 <run>_1.fastq.gz -2 <run>_2.fastq.gz \
  --un-conc <work>/unconc.fastq --un <work>/unpaired.fastq \
  -S <work>/aligned.sam
```

No read trimming and no base-quality filter is applied before Bowtie2.
Paired files stay in R1-then-R2 order. `--no-unal` omits unmapped reads
from the SAM. Stage 1 does not set `-L`, `--score-min`, or `-k`; those
belong to stage 2 (`-L 9`, `--score-min L,6,0.2`, `-k 200`) and to the
junction pass (`-L 10`, `--score-min L,1,0.70`, `-k 2000`).

`bowtie2-build` needs FASTA. For MG1655, production `materialize_ref`
writes the GenBank ORIGIN to FASTA because `MG1655.gbk` is not already a
FASTA. The diagnostic does the same extraction, checks the chromosome
SHA256 above, and names the record `NC_000913.3`. For BW25113 the existing
FASTA is the index input. That is the production reference path, not a
second reference.

Sorting matches `sam_to_sorted_bam`: noodles reads the stage-1 SAM into
memory and sorts by reference id, then alignment start. Production does
not write a `.bai`. The diagnostic may stream the BAM once and keep only
reads that overlap the target base. An index is optional for that stream.

Temporary files per run: Bowtie2 index, `aligned.sam`, unmatched FASTQ,
and `aligned.bam`. The in-memory SAM sort is the memory risk. P2a full
evidence, which also ran stage 2 and junctions, stayed under about 11 GB
RSS for libraries of this size. The diagnostic requests 32 GB so the
stage-1 SAM can be held without assuming that peak.

## Measurements

These counts are independent diagnostic measurements. They are not a claim
that the script reproduced every internal RA filter. The caller drops reads
with MAPQ below 10 (`UNIQUE_MAPQ`), sets MAPQ to 0 when the XS score is at
least the AS score, trims one aligned base from each end of a read
(`END_TRIM = 1`), and can mask homopolymer ends. RA stores base and strand
only. It does not filter on base quality. A Phred cutoff is therefore an
independent measurement, not caller-eligible depth.

For each target coordinate, on the stage-1 BAM:

| Count | Definition |
| --- | --- |
| Raw depth | Primary alignments (not secondary, not supplementary, not unmapped) whose CIGAR covers the 1-based base |
| REF / ALT / other | Base at that coordinate among the raw set. Deletion covering the base is other, not REF or ALT |
| MAPQ ≥ 10 depth | Same, after the XS≥AS rule sets MAPQ to 0 |
| Caller-shaped depth | MAPQ ≥ 10 set, excluding a base that is the first or last aligned base of that read |
| Q20 depth | Raw set with Phred ≥ 20 at the base. Independent of the caller |
| Strand counts | Plus and minus among the caller-shaped set |
| Base-quality distribution | Phred histogram for the raw base at the coordinate |
| Mapping-quality distribution | MAPQ histogram for raw alignments covering the coordinate |
| Clip or indel adjacent | A CIGAR I, D, or S operation within 5 bp of the coordinate. This is a flag, not a recovery window |
| Multi-mapping | Alignments with MAPQ < 10, or XS≥AS, covering the coordinate |

Allele fraction uses caller-shaped ALT count divided by caller-shaped depth.
Raw depth is reported beside it and is not substituted for it.

## Classification

Apply only after the counts exist. Do not use published mixture percentages.

| Class | Required measurement |
| --- | --- |
| `REFERENCE_DISAGREEMENT` | Chromosome base at the coordinate differs from the published reference allele |
| `INSUFFICIENT_COVERAGE` | Caller-shaped depth < 5 and raw ALT count is 0 |
| `LOW_QUALITY_SUPPORT` | Raw ALT count > 0, Q20 ALT count is 0, and Q20 depth ≥ 5 |
| `MAPPING_AMBIGUITY` | Caller-shaped ALT count is 0, low-MAPQ ALT count > 0, and caller-shaped depth < 5 |
| `UNRESOLVED` | Caller-shaped depth < 20 and ALT count is 0, or depth < 5 while some ALT remains, or the extraction failed |
| `ALT_ABSENT_WITH_ADEQUATE_DEPTH` | Caller-shaped depth ≥ 20 and caller-shaped ALT count is 0 |
| `LOW_ALT_FREQUENCY` | Caller-shaped depth ≥ 5 and ALT fraction < 0.75 |
| `ALT_PRESENT_ABOVE_CONSENSUS_THRESHOLD` | Caller-shaped depth ≥ 5 and ALT fraction ≥ 0.75 |

`ALT_PRESENT_ABOVE_CONSENSUS_THRESHOLD` means the allele is visible in this
pileup. It does not, by itself, identify a caller defect. The 0.75 and 5
values are the existing RA defaults used as labels for the measurement.
They are not new caller thresholds. Depth from 5 through 19 with zero ALT
stays `UNRESOLVED` rather than an absence claim.

## Execution plan, not started

Partition `qcpu_18i`. One task per priority-1 run. 8 CPU, 32 GB, 4 hour
wall time. At most two tasks at once. FASTQ paths and checksums are the
P2a provenance rows. The command is `bowtie2` as written above, then a
locus extract. It is not `prokadiff evidence`.

Keep `aligned.bam` until the locus table has one row per priority-1 allele
and the reference hash check passed. Do not delete the BAM in the same step
as the extract. Outputs stay under gitignored
`benchmark/real_world/results/prjna1088182/m6_p3a/<run>/`. Expected small
files: `locus_pileup.tsv`, `command.txt`, and the reference SHA256 check.
BAM, SAM, and the Bowtie2 index stay untracked.

Scratch per run: index under 20 MB, SAM and BAM on the order of a few
gigabytes for these 3–4 million spot libraries, unmatched FASTQ smaller
than the inputs. 40 GB free on the run directory is enough headroom.
The three gzip pairs are already on disk (about 0.5–0.6 GB each).

## Risks

Stage-1 BAM is the production RA input, but homopolymer masking and
deletion right-alignment are not in the first extract. Caller-shaped depth
can therefore differ from the internal column. The plan labels that gap.

In-memory SAM sorting can exceed the P2a RSS if the SAM is larger than the
final BAM. The 32 GB request is the mitigation. A failed extract leaves the
BAM in place.

Priority 2 alleles are not explained by the first three runs.

## Verdict

**M6-P3a PLAN READY FOR EXECUTION**

**M6 overall: ACTIVE**
