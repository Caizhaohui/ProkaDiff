#!/usr/bin/env bash
# Inspect deposited SRA read structure. No ProkaDiff, Bowtie2, or breseq.
set -euo pipefail
REPO_DIR="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff"
BIN="$REPO_DIR/data/tools/sra-tools/bin"
OUT="$REPO_DIR/benchmark/real_world/results/m6_p1a/read_structure"
mkdir -p "$OUT"
export PATH="$BIN:$PATH"

runs=(SRR28349244 SRR28357593 SRR28349243 SRR28357592)
for run in "${runs[@]}"; do
    sra="$REPO_DIR/data/prjna1088182/tmp/$run/$run/$run.sra"
    test -s "$sra"
    {
        echo "=== $run sra-stat --quick ==="
        sra-stat --quick "$sra"
        echo "=== $run sra-stat -x ==="
        sra-stat -x "$sra"
        echo "=== $run SEQUENCE columns ==="
        vdb-dump -T SEQUENCE -o "$sra" || true
        echo "=== $run first 5 SEQUENCE rows ==="
        vdb-dump -T SEQUENCE -R 1-5 -C READ_LEN,READ_TYPE,READ_START,SPOT_LEN -f tab "$sra" || \
            vdb-dump -T SEQUENCE -R 1-5 -f json "$sra" | head -c 4000
    } > "$OUT/${run}.structure.txt" 2>&1

    exp="$OUT/${run}_split_experiment"
    rm -rf "$exp"
    mkdir -p "$exp/split-files" "$exp/split-3"
    fastq-dump -X 8 --split-files -O "$exp/split-files" "$sra" > "$exp/split-files.log" 2>&1 || true
    fastq-dump -X 8 --split-3 -O "$exp/split-3" "$sra" > "$exp/split-3.log" 2>&1 || true
    {
        echo "=== split-files names ==="
        find "$exp/split-files" -type f -printf '%f %s\n' | sort
        echo "=== split-3 names ==="
        find "$exp/split-3" -type f -printf '%f %s\n' | sort
        echo "=== first headers split-files ==="
        for f in "$exp"/split-files/*.fastq; do
            [[ -f "$f" ]] || continue
            echo "-- $(basename "$f")"
            head -n 4 "$f"
        done
    } > "$OUT/${run}.export.txt"
done
echo "inspection complete"
