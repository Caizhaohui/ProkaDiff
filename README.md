# ProkDiff

**Prokaryotic genome diff for gene-editing QC — starter vs edited WGS**

[中文说明](README.zh.md)

Rust tool for whole-genome unintended-mutation QC after prokaryotic gene editing. Both the **starter strain** and the **edited strain** are sequenced (Illumina WGS). An NCBI / close genome is a **coordinate skeleton only**. Unintended mutations are those present in the edited sample relative to the starter, then labeled into three observational classes.

This is not mammalian GUIDE-seq, and it is not “every SNP versus MG1655.”

## Positioning


| | |
| --- | --- |
| Name | **ProkDiff** (`prokdiff` CLI) |
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

1. **`near_homolog`** — near a spacer homolog + that editor’s PAM (Cas9 default NGG; Cas12a default TTTV). Annotation, not sole evidence.
2. **`scattered_snv`** — remaining SNPs / short indels. Reported for every editor.
3. **`structural`** — MOB / JC, large deletions, plasmid-related events.

Classification lands in a later milestone. The engine MVP writes Genome Diff evidence only.

## CLI

Product (both strains required; engine runs twice; no classify yet):

```bash
prokdiff \
  --starter starter_R1.fastq.gz --starter starter_R2.fastq.gz \
  --edited  edited_R1.fastq.gz  --edited  edited_R2.fastq.gz \
  --ref     NC_000913.3.fa \
  --editor cas9 \
  --threads 8 \
  --outdir  results/prokdiff
```

Single-sample engine (oracle / parity jobs):

```bash
prokdiff evidence \
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

Engine MVP (stage 2): Bowtie2 wrapper, noodles BAM I/O, RA / MC / JC, `prokdiff evidence`. Product classify (the three classes above) is not implemented yet.

Layer-0 unit tests are in-memory. Synthetic FASTQ generation (`testdata/generate.sh`) and breseq parity jobs must run on Slurm partition **`qcpu_18i`**, not a login node.

## License

MIT. Methods: [breseq 0.40.1](https://breseq.barricklab.org/0.40.1/methods/). Scientific framing: Widney et al. 2024 bioRxiv [10.1101/2024.03.19.584922](https://doi.org/10.1101/2024.03.19.584922).
