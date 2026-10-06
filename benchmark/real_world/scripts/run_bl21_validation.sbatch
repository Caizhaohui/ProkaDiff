#!/usr/bin/env bash
#SBATCH --job-name=prokadiff_bl21_m5
#SBATCH --partition=qcpu_18i
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=benchmark/real_world/results/bl21_user/slurm_%j.log
#SBATCH --error=benchmark/real_world/results/bl21_user/slurm_%j.err

set -euo pipefail
export PATH="/hpcfs/fhome/caizhh/.conda/envs/prokadiff/bin:$PATH"
export RUST_BACKTRACE=1

REPO_DIR="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff"
BL21_DATA="/hpcfs/fhome/caizhh/19_BL21_edited"
MANIFEST="$REPO_DIR/benchmark/real_world/manifests/bl21_user.tsv"
INTENDED="$REPO_DIR/benchmark/real_world/truth/bl21_intended.tsv"
PROKADIFF="$REPO_DIR/target/release/prokadiff"
OUT_BASE="$REPO_DIR/benchmark/real_world/results/bl21_user"
REF="$BL21_DATA/BL21/ncbi_dataset/data/GCF_013167015.1/genomic.gbff"

cd "$REPO_DIR"
mkdir -p "$OUT_BASE"
# shellcheck source=attempt_peer_guides.sh
source "$REPO_DIR/benchmark/real_world/scripts/attempt_peer_guides.sh"
python3 "$REPO_DIR/benchmark/real_world/scripts/validate_bl21_registry.py" --data-root "$BL21_DATA"
test -x "$PROKADIFF"

manifest_value() {
    local sample="$1" name="$2"
    awk -F '\t' -v sample="$sample" -v name="$name" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == name) column = i } NR > 1 && $3 == sample { if (!column) exit 2; print $column; found = 1 } END { if (!found) exit 1 }' "$MANIFEST"
}

echo "Node: $(hostname)"
echo "Date: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
echo "ProkaDiff Commit: $(git rev-parse HEAD)"
echo "Bowtie2: $(bowtie2 --version | head -n 1)"
echo "Breseq Oracle: $(breseq --version | head -n 1)"

for sample in B21_2_1 B21_3_1 B21_4_1; do
    r1="$(manifest_value "$sample" r1)"
    r2="$(manifest_value "$sample" r2)"
    evidence_dir="$OUT_BASE/${sample}_evidence"
    mkdir -p "$evidence_dir"
    /usr/bin/time -v -o "$evidence_dir/telemetry.txt" \
        "$PROKADIFF" evidence --ref "$REF" --fastq "$r1" "$r2" \
        --threads 16 --outdir "$evidence_dir"
    if ! python3 "$REPO_DIR/benchmark/real_world/scripts/validate_bl21_products.py" \
        --diagnostic-sidecar "$evidence_dir/repeat_ambiguous_junctions.tsv" \
        > "$evidence_dir/diagnostic_sidecar.txt"; then
        echo "repeat-ambiguous sidecar for $sample is missing or is not diagnostic-only" >&2
        run_status=1
    fi

    oracle="$BL21_DATA/analysis_bl21_cas9_wgs/results/sv/$sample/augmented_full/breseq/output/output.gd"
    residual_dir="$evidence_dir/oracle_residuals"
    bash "$REPO_DIR/benchmark/real_world/scripts/preserve_gdtools_residuals.sh" "$evidence_dir/output.gd" "$oracle" "$residual_dir" prokadiff_vs_breseq
    python3 "$REPO_DIR/benchmark/real_world/scripts/review_target_locus.py" \
        --sample "$sample" --target-id edit_1 --seq-id NZ_CP053602.1 \
        --start 334876 --end 335735 --prokadiff-gd "$evidence_dir/output.gd" \
        --breseq-gd "$oracle" \
        --curated-truth "$REPO_DIR/benchmark/real_world/truth/bl21_curated_truth.tsv" \
        --output "$evidence_dir/target_locus_review.tsv"
    python3 "$REPO_DIR/benchmark/real_world/scripts/compare_variant_calls.py" \
        --test-gd "$evidence_dir/output.gd" --truth-gd "$oracle" --tsv \
        > "$evidence_dir/variant_metrics.tsv"
done

for sample in B21_2_1 B21_3_1 B21_4_1; do
    starter="$(manifest_value "$sample" starter_id)"
    starter_r1="$(manifest_value "$starter" r1)"
    starter_r2="$(manifest_value "$starter" r2)"
    edited_r1="$(manifest_value "$sample" r1)"
    edited_r2="$(manifest_value "$sample" r2)"
    run_peer_guide() {
        local guide_id="$1" spacer product_dir failed=0
        case "$guide_id" in
            BL21_G1) spacer="$(manifest_value "$sample" spacer_1)" ;;
            BL21_G2) spacer="$(manifest_value "$sample" spacer_2)" ;;
            *) echo "unknown guide $guide_id" >&2; return 1 ;;
        esac
        product_dir="$OUT_BASE/${sample}_vs_${starter}/${guide_id}"
        mkdir -p "$(dirname "$product_dir")"
        if ! /usr/bin/time -v -o "$OUT_BASE/${sample}_vs_${starter}/${guide_id}.telemetry.txt" \
            "$PROKADIFF" \
            --starter "$starter_r1" "$starter_r2" \
            --edited "$edited_r1" "$edited_r2" \
            --ref "$REF" --intended "$INTENDED" --editor cas9 \
            --spacer "$spacer" --pam NGG --schema-version v2 \
            --threads 16 --outdir "$product_dir"; then
            return 1
        fi

        residual_dir="$product_dir/oracle_residuals"
        sample_oracle="$BL21_DATA/analysis_bl21_cas9_wgs/results/sv/$sample/augmented_full/breseq/output/output.gd"
        starter_oracle="$BL21_DATA/analysis_bl21_cas9_wgs/results/sv/$starter/augmented_full/breseq/output/output.gd"
        bash "$REPO_DIR/benchmark/real_world/scripts/preserve_gdtools_residuals.sh" "$product_dir/edited.gd" "$product_dir/starter.gd" "$residual_dir" prokadiff_pair || failed=1
        bash "$REPO_DIR/benchmark/real_world/scripts/preserve_gdtools_residuals.sh" "$sample_oracle" "$starter_oracle" "$residual_dir" breseq_pair || failed=1
        if [[ -f "$residual_dir/prokadiff_pair_A_minus_B.gd" && -f "$residual_dir/breseq_pair_A_minus_B.gd" ]]; then
            bash "$REPO_DIR/benchmark/real_world/scripts/preserve_gdtools_residuals.sh" "$residual_dir/prokadiff_pair_A_minus_B.gd" "$residual_dir/breseq_pair_A_minus_B.gd" "$residual_dir" pair_parity || failed=1
        fi

        if ! python3 "$REPO_DIR/benchmark/real_world/scripts/compare_intended_edits.py" \
            --edit-outcomes "$product_dir/edit_outcomes.tsv" --comparison-role PEER_COMPARATOR --tsv \
            > "$product_dir/intended_metrics.tsv"; then failed=1; fi
        if ! python3 "$REPO_DIR/benchmark/real_world/scripts/validate_reports.py" \
            --strict "$product_dir/report.md" > "$product_dir/report_validation.txt"; then failed=1; fi
        if ! python3 "$REPO_DIR/benchmark/real_world/scripts/validate_bl21_products.py" \
            --products "$product_dir" --sample "$sample" --guide-id "$guide_id" \
            --comparison-role PEER_COMPARATOR --target-review "$OUT_BASE/${sample}_evidence/target_locus_review.tsv" \
            --diagnostic-sidecar "$OUT_BASE/${sample}_evidence/repeat_ambiguous_junctions.tsv" \
            > "$product_dir/product_validation.tsv"; then failed=1; fi
        return "$failed"
    }
    attempt_peer_guides run_peer_guide
    if [[ "$GUIDE_FAILURES" -ne 0 ]]; then
        run_status=1
    fi
done

if [[ ${run_status:-0} -ne 0 ]]; then
    echo "BL21 validation failed for at least one comparison; all guide attempts completed." >&2
    exit 1
fi
echo "BL21 validation matrix completed: 3 PEER_COMPARATOR pairs x 2 independent guides."
