# Benchmark results

Measured ProkDiff ↔ breseq comparisons (output / wall-clock / peak RSS) go here.

Do not invent numbers. Jobs run on Slurm partition `qcpu_18i`. Generated CSVs and job trees are gitignored except files explicitly allowlisted in `.gitignore`.

## Tracked notes

- [clonal_leftover.md](clonal_leftover.md) — job **2369858** leftover classification (breseq **0.39.0**, node `bnode5.tibhpc.net`). Layer 2 will be **re-run after breseq 0.40.x**; do not resubmit 2369858 for this engine pass.

## Job pointers (local only)

- `scripts/submit_parity_clonal.sh` → Slurm job **2369858** → `clonal_2369858/` and `clonal.csv` (gitignored)
- `scripts/submit_parity_example6.sh` → Slurm job **2369859** → `example6_<jobid>/` and `example6.csv` (gitignored)

Oracle env recorded in `testdata/VERSIONS.txt` is breseq **0.39.0** + bowtie2 2.5.4. Do not label those runs as 0.40.x. Do not write speedup claims in specs.
