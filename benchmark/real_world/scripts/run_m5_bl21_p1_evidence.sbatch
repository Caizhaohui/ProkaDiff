#!/usr/bin/env bash
#SBATCH --job-name=pkd_m5_bl21_p1
#SBATCH --partition=qcpu_18i
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=benchmark/real_world/results/bl21_user/slurm_m5_bl21_p1_%j.log
#SBATCH --error=benchmark/real_world/results/bl21_user/slurm_m5_bl21_p1_%j.err

# Reference-relative single-sample evidence. Not a parent or peer comparison.
set -euo pipefail
export PATH="/hpcfs/fhome/caizhh/.conda/envs/prokadiff/bin:$PATH"
export RUST_BACKTRACE=1

SAMPLE="${1:?sample id required}"
case "$SAMPLE" in
  B21_4_1|B21_2_1) ;;
  *) echo "unsupported sample for M5-BL21-P1: $SAMPLE" >&2; exit 2 ;;
esac

REPO_DIR="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff"
cd "$REPO_DIR"
BL21_DATA="/hpcfs/fhome/caizhh/19_BL21_edited"
REF="$BL21_DATA/BL21/ncbi_dataset/data/GCF_013167015.1/genomic.gbff"
R1="$BL21_DATA/TCN202608300034/00.CleanData/${SAMPLE}/${SAMPLE}_1.clean.fq.gz"
R2="$BL21_DATA/TCN202608300034/00.CleanData/${SAMPLE}/${SAMPLE}_2.clean.fq.gz"
OUT_DIR="$REPO_DIR/benchmark/real_world/results/bl21_user/${SAMPLE}_evidence_m5_bl21_p1"

if [[ -e "$OUT_DIR/output.gd" || -e "$OUT_DIR/repeat_ambiguous_junctions.tsv" ]]; then
  echo "refusing to overwrite existing evidence in $OUT_DIR" >&2
  exit 3
fi
mkdir -p "$OUT_DIR"
test -s "$REF"
test -s "$R1"
test -s "$R2"

echo "=== M5-BL21-P1 ${SAMPLE} reference-relative evidence ==="
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-manual}"
echo "NODE: $(hostname)"
echo "CPUS: ${SLURM_CPUS_PER_TASK:-16}"
echo "MEM_REQUESTED: 48G"
echo "TIME_LIMIT: 02:00:00"
echo "GIT_HEAD: $(git rev-parse HEAD)"
echo "SR1_COMMIT: f30f93be5439d2f8eb32f79719fe8482713fcf5a"
echo "REFERENCE: $REF"
echo "REFERENCE_ACCESSION: GCF_013167015.1"
echo "R1: $R1"
echo "R2: $R2"
echo "OUT_DIR: $OUT_DIR"
echo "COMPARISON_ROLE: REFERENCE_RELATIVE_SINGLE_SAMPLE"

echo "=== BUILD RELEASE BINARY FROM CURRENT SOURCE ==="
cargo build --release -p prokadiff --bin prokadiff
PROKADIFF="$REPO_DIR/target/release/prokadiff"

/usr/bin/time -v -o "$OUT_DIR/telemetry.txt" \
  "$PROKADIFF" evidence \
    --ref "$REF" \
    --fastq "$R1" "$R2" \
    --threads 16 \
    --outdir "$OUT_DIR"

python3 - "$OUT_DIR" << 'PY'
import csv
import sys
from collections import Counter
from pathlib import Path

out = Path(sys.argv[1])
sidecar = out / "repeat_ambiguous_junctions.tsv"
rows = list(csv.DictReader(sidecar.open(), delimiter="\t"))
print(f"repeat_diagnostic_rows: {len(rows)}")
print("states:", dict(Counter(row["state"] for row in rows)))
if rows:
    print("max_placement_count:", max(int(row["placement_count"] or 0) for row in rows))
    print("max_overflow_count:", max(int(row["overflow_count"] or 0) for row in rows))
    print("max_repeat_family_count_field:", max(int(row["repeat_family_count"] or 0) for row in rows))
print("evidence complete")
PY
