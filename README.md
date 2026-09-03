# ProkaDiff

**Prokaryotic genome diff for gene-editing QC — starter vs edited WGS**

[中文说明](README.zh.md)

Rust tool for whole-genome unintended-mutation QC after prokaryotic gene editing. Both the **starter strain** and the **edited strain** are sequenced (Illumina WGS). An NCBI / close genome is a **coordinate skeleton only**. Unintended mutations are those present in the edited sample relative to the starter, then labeled into three observational classes.

This is not mammalian GUIDE-seq, and it is not “every SNP versus MG1655.”

## Positioning


| | |
| --- | --- |
| Name | **ProkaDiff** (`prokadiff` CLI) |
| Reads (v0.1) | Illumina short reads, clonal consensus |
| Engine | External **Bowtie2** + Rust RA / MC / JC (clean-room vs [breseq methods](https://breseq.barricklab.org/0.40.1/methods/)). System `breseq` is **oracle only**, never the runtime engine. |
| Diff | `M_edited \ M_parent \ {intended}` |
| Editors (v0.1) | `--editor cas9 \| cas12a \| dsb` (`cast` / `is110` rejected with a roadmap message) |

## Dual-sample requirement

`--starter` and `--edited` FASTQ are **both required**. Missing either exits non-zero with a message that starter-strain WGS is mandatory. Do not treat all differences versus NCBI as unintended edits.

- Repeatable flags: **two consecutive files under the same flag = one PE pair** (R1 then R2); a single file = SE.
- `--intended` is optional; if omitted, every post-subtract call is unintended.
- `--ref` is repeatable (chromosome + plasmids + donor backbone).

## Three classes (observational)

Assignment order is locked: **structural (3) → near_homolog (1) → scattered_snv (2)**. One mutation → one class. Class (2) is **never** disabled because “this editor does not make DSBs.”

1. **`near_homolog`** — near a spacer homolog + that editor’s PAM (Cas9 default NGG; Cas12a default TTTV; ≤4 mismatches; default 50 bp). Annotation, not sole evidence. `--editor dsb` does not assign this class.
2. **`scattered_snv`** — remaining SNPs / short indels. Reported for every editor. Optional hypothesis `sos_widney2014` unless `--no-hypothesis`.
3. **`structural`** — MOB / JC, large deletions (>2 bp), AMP / CON.

## CLI

Product (both strains required; classify after subtract):

```bash
prokadiff \
  --starter starter_R1.fastq.gz --starter starter_R2.fastq.gz \
  --edited  edited_R1.fastq.gz  --edited  edited_R2.fastq.gz \
  --ref     NC_000913.3.fa \
  --intended intended.tsv \
  --editor cas9 \
  --spacer  NNNNNNNNNNNNNNNNNNNN \
  --threads 8 \
  --outdir  results/prokadiff
```

Writes `starter.gd`, `edited.gd`, `unintended.tsv`, and `summary.txt`.

`summary.txt` records whether `--intended` was supplied (`intended_provided`), how many rows were declared (`intended_declared`), how many post-subtract mutations were masked (`intended_observed`), and a roll-up `intended_status` (`all_observed` / `partial` / `none_observed` / `NA`). SNP `ref` in the TSV is the reference base; JC rows also have `side2_seq_id` / `side2_position`. Full contract: [docs/schema.md](docs/schema.md).

Single-sample engine (oracle / parity jobs):

```bash
prokadiff evidence \
  --ref genome.fa \
  --fastq sample_R1.fastq --fastq sample_R2.fastq \
  --threads 8 \
  --outdir results/sample
```

I/O contract: [docs/schema.md](docs/schema.md). Oracle rules (output + wall-clock + peak RSS): [docs/parity.md](docs/parity.md).

## What v0.1 does not do

- CAST / bRNA-IS110-specific homolog and donor interpretation
- Population polymorphism, long reads, GUIDE-seq
- Shelling out to system `breseq` as the engine
- Built-in assembly; rewriting Bowtie2 / HMMER / ESM
- Shipping full deep *E. coli* FASTQ or Widney/Copley SRA (PRJNA1088182) in git
- Reporting NCBI-only diffs as gene-editing off-targets

## Status

Stage 4 report layer: `summary.txt` includes declared / observed / status for intended edits; `unintended.tsv` writes SNP reference bases and JC side-2 coordinates. Engine MVP (`prokadiff evidence`) is unchanged. Layer-0 tests are in-memory. Layer-1 / layer-2 engine jobs must run on Slurm **`qcpu_18i`**.

## License

MIT. Methods: [breseq 0.40.1](https://breseq.barricklab.org/0.40.1/methods/). Scientific framing: Widney et al. 2024 bioRxiv [10.1101/2024.03.19.584922](https://doi.org/10.1101/2024.03.19.584922).
