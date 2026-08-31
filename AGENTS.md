# Agent rules for ProkDiff

Follow [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md). Project name is **ProkDiff**
(crates: `prokdiff*`；二进制一律 `prokdiff`).

## Rule 1

Execute **one milestone at a time**. Do not implement later milestones early.

## Rule 2

Before starting a milestone:

1. Read `DEVELOPMENT_PLAN.md`
2. Read [docs/schema.md](docs/schema.md) and [docs/parity.md](docs/parity.md) if touching I/O or oracle tests
3. Read current code and tests
4. Run baseline tests (`cargo test --workspace` when the workspace exists)

## Rule 3

Each milestone must ship **implementation + tests + documentation** together.

## Rule 4

**Engine is path B only.** Call system `bowtie2` for alignment; implement RA / MC / JC in Rust (`prokdiff-evidence`). **Never** shell out to `breseq` as a runtime engine. breseq is **oracle only** for output parity and for wall-clock / resource comparison on compute nodes.

## Rule 5

**Clean-room vs breseq.** Implement from the published [breseq methods](https://breseq.barricklab.org/0.40.1/methods/). **Do not copy breseq source code** (GPL-2.0). Oracle comparison must cover **(1) output** (`.gd` + `gdtools SUBTRACT`), **(2) wall-clock runtime**, and **(3) compute resources** (peak RSS / CPU). Never claim parity or speed without measured numbers on a compute node.

## Rule 6

CLI **requires both** `--starter` and `--edited` FASTQ. Missing either → non-zero exit with a message that starter-strain WGS is mandatory. **Never** report all differences vs NCBI as unintended edits.

## Rule 7

Classification order is locked: **(3) structural → (1) near_homolog → (2) scattered_snv**. One mutation → one class. **Never** disable class (2) because “this editor does not make DSBs.”

## Rule 8

First-period `--editor`: `cas9` | `cas12a` | `dsb` only. `cast` / `is110` must error and point to the roadmap. Hypothesis columns are optional annotations, not facts. Never hard-code “CAST / IS110 = no genome-wide mutagenesis.”

## Rule 9

BAM I/O uses **`noodles`**. Do not add `rust-htslib` unless `DEVELOPMENT_PLAN.md` explicitly changes this. Every new dependency needs a written rationale.

## Rule 10

**Do not** put Widney/Copley SRA (PRJNA1088182) in testdata. **Do not** commit full *E. coli* deep FASTQ. Use `testdata/generate.sh` or git-lfs for synthetic reads. Private paired WGS stays in `testdata/private/` (gitignored).

## Rule 11

**Do not run test workloads on the login node.** Login node is only for edit / `cargo fmt` / `cargo clippy` / compile checks. Any run that invokes bowtie2, breseq, FASTQ I/O, synthetic genome fixtures, layer 1–3, or ProkDiff↔breseq timing must be submitted to Slurm partition **`qcpu_18i`** with **adequate** `cpus-per-task` / memory (enough for a fair side-by-side vs breseq; do not under-request). Pure in-memory layer-0 unit tests (no external binaries) may run in CI or briefly on login; prefer queue when unsure.

## Rule 12

Pin **breseq 0.40.x** + matching bowtie2 in `testdata/VERSIONS.txt` for oracle runs. Do not retarget 0.35.4 for literature reproduction.

## Rule 13

`--starter` / `--edited` may repeat. **Two consecutive files under the same flag = one PE pair** (R1 then R2); a single file = SE. `--intended` is optional (omit → all post-subtract calls are unintended).

## Rule 14

Never modify user resequencing source data. Real-data tests are **read-only**. Outputs go under `--outdir` or scratch.

## Rule 15

If real-data paths or breseq/bowtie2 binaries are missing: stop layer 2/3 work, continue pure unit work. Never fabricate parity, wall-clock, or resource numbers. Write measured comparisons under `benchmark/results/` (see [docs/parity.md](docs/parity.md)).
