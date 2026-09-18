#!/usr/bin/env bash
#SBATCH --job-name=pkd_b21_4_p13
#SBATCH --partition=qcpu_18i
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=benchmark/real_world/results/bl21_user/slurm_p1_3_b21_4_%j.log
#SBATCH --error=benchmark/real_world/results/bl21_user/slurm_p1_3_b21_4_%j.err

set -euo pipefail
export PATH="/hpcfs/fhome/caizhh/.conda/envs/prokadiff/bin:$PATH"
export RUST_BACKTRACE=1

REPO_DIR="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff"
cd "$REPO_DIR"

PROKADIFF="$REPO_DIR/target/release/prokadiff"
BL21_DATA="/hpcfs/fhome/caizhh/19_BL21_edited"
REF="$BL21_DATA/BL21/ncbi_dataset/data/GCF_013167015.1/genomic.gbff"
B21_4_R1="$BL21_DATA/TCN202608300034/00.CleanData/B21_4_1/B21_4_1_1.clean.fq.gz"
B21_4_R2="$BL21_DATA/TCN202608300034/00.CleanData/B21_4_1/B21_4_1_2.clean.fq.gz"
ORACLE_B21_4="$BL21_DATA/analysis_bl21_cas9_wgs/results/sv/B21_4_1/augmented_full/breseq/output/output.gd"

OUT_DIR="$REPO_DIR/benchmark/real_world/results/bl21_user/B21_4_1_evidence_p1_3"
mkdir -p "$OUT_DIR"

echo "=== B21_4_1 P1-3 verification run ==="
echo "Node: $(hostname)"
echo "Date: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
echo "ProkaDiff Commit: $(git rev-parse HEAD)"

/usr/bin/time -v -o "$OUT_DIR/telemetry.txt" \
    "$PROKADIFF" evidence \
        --ref "$REF" \
        --fastq "$B21_4_R1" "$B21_4_R2" \
        --threads 16 \
        --outdir "$OUT_DIR"

echo "evidence complete."

python3 "$REPO_DIR/benchmark/real_world/scripts/compare_variant_calls.py" \
    --test-gd "$OUT_DIR/output.gd" \
    --truth-gd "$ORACLE_B21_4" \
    --tsv | tee "$OUT_DIR/variant_metrics.tsv"
