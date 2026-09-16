# ProkaDiff development environment

The project development environment is the Conda environment `prokadiff` at
`/hpcfs/fhome/caizhh/.conda/envs/prokadiff` on the shared HPC system. Its
toolchain is pinned in [`environment.yml`](../environment.yml): Rust 1.98.1,
rust-analyzer 2026.08.24, rust-src 1.98.1, Python 3.12, Bowtie2 2.5.4,
breseq 0.40.2, and wgsim 1.0. `breseq` and `gdtools` are oracle-only tools;
ProkaDiff never invokes breseq at runtime.

## Run the Rust checks

Use the environment explicitly so shell startup files cannot silently switch
back to the base environment:

```bash
ENV_PREFIX=/hpcfs/fhome/caizhh/.conda/envs/prokadiff
mamba run -p "$ENV_PREFIX" cargo fmt --all --check
mamba run -p "$ENV_PREFIX" cargo clippy --workspace --all-targets -- -D warnings
mamba run -p "$ENV_PREFIX" cargo test --workspace
```

The default workspace test command does not run FASTQ/Bowtie2 tests. Tests that
invoke Bowtie2, breseq, gdtools, wgsim, or process sequencing data must still be
submitted to Slurm partition `qcpu_18i` with adequate resources; use the
existing `scripts/submit_*.sh` wrappers for those workloads.

## Activate interactively

```bash
mamba activate /hpcfs/fhome/caizhh/.conda/envs/prokadiff
```

The project-scoped `.codex/lsp-client.json` routes Rust LSP startup through the
same environment, so `rust-analyzer` and its matching Rust toolchain are used
for diagnostics, definitions, references, and rename checks.
