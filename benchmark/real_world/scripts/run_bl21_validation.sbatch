#!/usr/bin/env bash
#SBATCH --job-name=prokadiff_bl21
#SBATCH --partition=qcpu_18i
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=benchmark/real_world/results/bl21_user/slurm_%j.log
#SBATCH --error=benchmark/real_world/results/bl21_user/slurm_%j.err

set -euo pipefail

# Ensure environment
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

B21_2_R1="$BL21_DATA/TCN202608300034/00.CleanData/B21_2_1/B21_2_1_1.clean.fq.gz"
B21_2_R2="$BL21_DATA/TCN202608300034/00.CleanData/B21_2_1/B21_2_1_2.clean.fq.gz"

ORACLE_B21_3="$BL21_DATA/analysis_bl21_cas9_wgs/results/sv/B21_3_1/augmented_full/breseq/output/output.gd"
ORACLE_B21_4="$BL21_DATA/analysis_bl21_cas9_wgs/results/sv/B21_4_1/augmented_full/breseq/output/output.gd"

OUT_BASE="$REPO_DIR/benchmark/real_world/results/bl21_user"
mkdir -p "$OUT_BASE"

echo "=== Starting ProkaDiff BL21 Real-World Validation ==="
echo "Node: $(hostname)"
echo "Date: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
echo "ProkaDiff Commit: $(git rev-parse HEAD)"
echo "Bowtie2: $(bowtie2 --version | head -n 1)"
echo "Breseq Oracle: $(breseq --version | head -n 1)"

# -------------------------------------------------------------
# 1. Single-sample Evidence Engine on B21_3_1
# -------------------------------------------------------------
echo "--> [1/4] Running ProkaDiff Evidence Engine on B21_3_1..."
DIR_B21_3="$OUT_BASE/B21_3_1_evidence"
mkdir -p "$DIR_B21_3"

/usr/bin/time -v -o "$DIR_B21_3/telemetry.txt" \
    "$PROKADIFF" evidence \
        --ref "$REF" \
        --fastq "$B21_3_R1" "$B21_3_R2" \
        --threads 16 \
        --outdir "$DIR_B21_3"

echo "B21_3_1 evidence complete."

# -------------------------------------------------------------
# 2. Single-sample Evidence Engine on B21_4_1
# -------------------------------------------------------------
echo "--> [2/4] Running ProkaDiff Evidence Engine on B21_4_1..."
DIR_B21_4="$OUT_BASE/B21_4_1_evidence"
mkdir -p "$DIR_B21_4"

/usr/bin/time -v -o "$DIR_B21_4/telemetry.txt" \
    "$PROKADIFF" evidence \
        --ref "$REF" \
        --fastq "$B21_4_R1" "$B21_4_R2" \
        --threads 16 \
        --outdir "$DIR_B21_4"

echo "B21_4_1 evidence complete."

# -------------------------------------------------------------
# 3. Dual-strain Differential Audit: B21_3_1 (edited) vs B21_4_1 (starter)
# -------------------------------------------------------------
echo "--> [3/4] Running ProkaDiff Differential Audit (B21_3_1 edited vs B21_4_1 starter)..."
DIR_DIFF="$OUT_BASE/B21_3_1_diff"
mkdir -p "$DIR_DIFF"

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

# -------------------------------------------------------------
# 4. Evaluations & Release Gate Checks
# -------------------------------------------------------------
echo "--> [4/4] Running evaluation scripts..."

# Variant comparison vs breseq oracle
python3 "$REPO_DIR/benchmark/real_world/scripts/compare_variant_calls.py" \
    --test-gd "$DIR_B21_3/output.gd" \
    --truth-gd "$ORACLE_B21_3" \
    --tsv > "$DIR_B21_3/variant_metrics.tsv"

python3 "$REPO_DIR/benchmark/real_world/scripts/compare_variant_calls.py" \
    --test-gd "$DIR_B21_4/output.gd" \
    --truth-gd "$ORACLE_B21_4" \
    --tsv > "$DIR_B21_4/variant_metrics.tsv"

# Intended edit evaluation & false_complete gate check
python3 "$REPO_DIR/benchmark/real_world/scripts/compare_intended_edits.py" \
    --audit-json "$DIR_DIFF/audit.json" \
    --expected-status UNEXPECTED_STRUCTURE \
    --tsv > "$DIR_DIFF/intended_metrics.tsv"

# Human report validation
python3 "$REPO_DIR/benchmark/real_world/scripts/validate_reports.py" \
    "$DIR_DIFF/report.md" > "$DIR_DIFF/report_validation.txt"

# Extract telemetry stats
python3 -c "
import sys, re
def parse_time(fpath):
    wall = 0.0
    rss = 0.0
    with open(fpath) as f:
        for line in f:
            if 'Elapsed (wall clock) time' in line:
                m = re.search(r'(\d+):(\d+\.?\d*)', line)
                if m:
                    wall = int(m.group(1)) * 60 + float(m.group(2))
            elif 'Maximum resident set size' in line:
                m = re.search(r'(\d+)', line)
                if m:
                    rss = float(m.group(1)) / 1024.0 # MB
    return wall, rss

wall, rss = parse_time('$DIR_DIFF/telemetry.txt')
with open('$DIR_DIFF/summary_metrics.tsv', 'w') as out:
    out.write(f'median_runtime_sec\t{wall:.1f}\n')
    out.write(f'median_rss_mb\t{rss:.1f}\n')
"

echo "=== All Steps Completed Successfully ==="
