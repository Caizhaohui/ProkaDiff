# Benchmark evidence

Only compact Markdown records in [`verified/`](verified/) are publishable benchmark evidence. A record must identify one `qcpu_18i` Slurm job and node, input workload, aligned thread count, tool versions, bidirectional comparator success, and wall-clock plus peak RSS for both tools.

`scripts/submit_parity_clonal.sh` creates `verified/clonal_<jobid>.md` only after `layer2_gd_compare.py` exits successfully and both `/usr/bin/time -v` files are complete. It also requires `SLURM_JOB_PARTITION=qcpu_18i` and checks breseq and Bowtie2 against `testdata/VERSIONS.txt`; a failed comparison, missing field, wrong partition, or version mismatch leaves no record.

Raw job trees, generated CSVs, logs, and incomplete historical measurements remain ignored. They may be used locally for diagnosis, but cannot support claims about parity, recall, resource use, or speed.
