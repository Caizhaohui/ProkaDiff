# Benchmark results

Measured ProkDiff ↔ breseq comparisons (output / wall-clock / peak RSS) go here.

Do not invent numbers. Jobs run on Slurm partition `qcpu_18i`. Generated CSVs are gitignored except files explicitly allowlisted in `.gitignore`.

Layer-2 engine jobs were submitted (not yet finished):

- `scripts/submit_parity_clonal.sh` → Slurm job **2369858** → `clonal_<jobid>/` and `clonal.csv`
- `scripts/submit_parity_example6.sh` → Slurm job **2369859** → `example6_<jobid>/` and `example6.csv`

Do not fill speedup claims until those tables exist. Oracle env is breseq **0.39.0** (see `testdata/VERSIONS.txt`); do not label the run as 0.40.x.
