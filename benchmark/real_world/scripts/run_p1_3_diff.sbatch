#!/usr/bin/env bash
#SBATCH --job-name=pkd_diff_p13
#SBATCH --partition=qcpu_18i
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=benchmark/real_world/results/bl21_user/slurm_p1_3_diff_%j.log
#SBATCH --error=benchmark/real_world/results/bl21_user/slurm_p1_3_diff_%j.err

set -euo pipefail
export PATH="/hpcfs/fhome/caizhh/.conda/envs/prokadiff/bin:$PATH"
export RUST_BACKTRACE=1

REPO_DIR="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff"
cd "$REPO_DIR"

PROKADIFF="$REPO_DIR/target/release/prokadiff"
BL21_DATA="/hpcfs/fhome/caizhh/19_BL21_edited"
REF="$BL21_DATA/BL21/ncbi_dataset/data/GCF_013167015.1/genomic.gbff"
INTENDED="$REPO_DIR/benchmark/real_world/truth/bl21_intended.tsv"

B21_3_R1="$BL21_DATA/TCN202608300034/00.CleanData/B21_3_1/B21_3_1_1.clean.fq.gz"
B21_3_R2="$BL21_DATA/TCN202608300034/00.CleanData/B21_3_1/B21_3_1_2.clean.fq.gz"

B21_4_R1="$BL21_DATA/TCN202608300034/00.CleanData/B21_4_1/B21_4_1_1.clean.fq.gz"
B21_4_R2="$BL21_DATA/TCN202608300034/00.CleanData/B21_4_1/B21_4_1_2.clean.fq.gz"

DIR_DIFF="$REPO_DIR/benchmark/real_world/results/bl21_user/B21_3_1_diff_p1_3"
mkdir -p "$DIR_DIFF"

echo "=== ProkaDiff Differential Audit P1-3 (B21_3_1 edited vs B21_4_1 starter) ==="
echo "Node: $(hostname)"
echo "Date: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
echo "ProkaDiff Commit: $(git rev-parse HEAD)"

/usr/bin/time -v -o "$DIR_DIFF/telemetry.txt" \
    "$PROKADIFF" \
        --starter "$B21_4_R1" "$B21_4_R2" \
        --edited "$B21_3_R1" "$B21_3_R2" \
        --ref "$REF" \
        --intended "$INTENDED" \
        --editor cas9 \
        --spacer "ATTTCGCTGGTGGTCAGATG" \
        --pam "NGG" \
        --threads 16 \
        --outdir "$DIR_DIFF"

echo "Differential audit complete."

python3 "$REPO_DIR/benchmark/real_world/scripts/compare_intended_edits.py" \
    --audit-json "$DIR_DIFF/edit_outcomes.tsv" \
    --expected-status UNEXPECTED_STRUCTURE \
    --tsv | tee "$DIR_DIFF/intended_metrics.tsv"

python3 "$REPO_DIR/benchmark/real_world/scripts/validate_reports.py" \
    "$DIR_DIFF/report.md" | tee "$DIR_DIFF/report_validation.txt"
