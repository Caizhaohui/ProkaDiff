#!/usr/bin/env bash
#SBATCH --job-name=pkd_diff_p13_m5
#SBATCH --partition=qcpu_18i
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=56G
#SBATCH --time=12:00:00
#SBATCH --output=benchmark/real_world/results/bl21_user/slurm_p1_3_diff_%j.log
#SBATCH --error=benchmark/real_world/results/bl21_user/slurm_p1_3_diff_%j.err

set -euo pipefail
export PATH="/hpcfs/fhome/caizhh/.conda/envs/prokadiff/bin:$PATH"
export RUST_BACKTRACE=1

REPO_DIR="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff"
BL21_DATA="/hpcfs/fhome/caizhh/19_BL21_edited"
MANIFEST="$REPO_DIR/benchmark/real_world/manifests/bl21_user.tsv"
INTENDED="$REPO_DIR/benchmark/real_world/truth/bl21_intended.tsv"
PROKADIFF="$REPO_DIR/target/release/prokadiff"
REF="$BL21_DATA/BL21/ncbi_dataset/data/GCF_013167015.1/genomic.gbff"
REF_FASTA="$BL21_DATA/BL21/ncbi_dataset/data/GCF_013167015.1/GCF_013167015.1_ASM1316701v1_genomic.fna"

EDITED="${1:?edited sample required}"
COMPARATOR="${2:?comparator sample required}"
case "${EDITED}:${COMPARATOR}" in
    B21_3_1:B21_4_1|B21_4_1:B21_3_1|B21_2_1:B21_3_1) ;;
    *)
        echo "unsupported PEER_COMPARATOR direction: ${EDITED} vs ${COMPARATOR}" >&2
        exit 2
        ;;
esac

cd "$REPO_DIR"
# shellcheck source=attempt_peer_guides.sh
source "$REPO_DIR/benchmark/real_world/scripts/attempt_peer_guides.sh"
OUT_BASE="$REPO_DIR/benchmark/real_world/results/bl21_user/m5_p3/${EDITED}_vs_${COMPARATOR}"
RUN_ID="${SLURM_JOB_ID:-manual_$(date -u +%Y%m%dT%H%M%SZ)}"
RUN_DIR="$OUT_BASE/validation_${RUN_ID}"
mkdir -p "$RUN_DIR"
registry_report="$(python3 "$REPO_DIR/benchmark/real_world/scripts/validate_bl21_registry.py" --data-root "$BL21_DATA" --identity-output "$RUN_DIR/reference_contig_identity.tsv")"
test -x "$PROKADIFF"

manifest_value() {
    local sample="$1" name="$2"
    awk -F '\t' -v sample="$sample" -v name="$name" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == name) column = i } NR > 1 && $3 == sample { if (!column) exit 2; print $column; found = 1 } END { if (!found) exit 1 }' "$MANIFEST"
}

starter_r1="$(manifest_value "$COMPARATOR" r1)"
starter_r2="$(manifest_value "$COMPARATOR" r2)"
edited_r1="$(manifest_value "$EDITED" r1)"
edited_r2="$(manifest_value "$EDITED" r2)"
manifest_version="$(manifest_value "$EDITED" registry_version)"
comparison_role="$(manifest_value "$EDITED" comparison_role)"
parent_verified="$(manifest_value "$EDITED" biological_parent_verified)"
reference_path="$(manifest_value "$EDITED" reference)"
manifest_comparator="$(manifest_value "$EDITED" starter_id)"
[[ "$manifest_version" == BL21_INPUTS_V2 ]]
[[ "$comparison_role" == PEER_COMPARATOR ]]
[[ "$parent_verified" == false ]]
if [[ "$manifest_comparator" != "$COMPARATOR" ]]; then
    echo "manifest starter_id ($manifest_comparator) does not match supplied comparator ($COMPARATOR); starter_id is a computational comparator, not a biological parent" >&2
    exit 2
fi

printf 'key\tvalue\n' > "$RUN_DIR/run_identity.tsv"
printf 'slurm_job_id\t%s\n' "$RUN_ID" >> "$RUN_DIR/run_identity.tsv"
printf 'prokadiff_commit\t%s\n' "$(git rev-parse HEAD)" >> "$RUN_DIR/run_identity.tsv"
printf 'worktree_dirty\t%s\n' "$(if [[ -n "$(git status --porcelain)" ]]; then echo true; else echo false; fi)" >> "$RUN_DIR/run_identity.tsv"
printf 'manifest_version\t%s\ncomparison_role\t%s\nbiological_parent_verified\t%s\n' "$manifest_version" "$comparison_role" "$parent_verified" >> "$RUN_DIR/run_identity.tsv"
printf 'reference_accession\tGCF_013167015.1\nreference_path\t%s\nreference_sha256\t%s\ngenbank_sha256\t%s\n' \
    "$reference_path" "$(sha256sum "$REF_FASTA" | awk '{print $1}')" "$(sha256sum "$REF" | awk '{print $1}')" >> "$RUN_DIR/run_identity.tsv"
printf 'contig_identity_validation\t%s\ncontig_identity_file\t%s\n' "$registry_report" "$RUN_DIR/reference_contig_identity.tsv" >> "$RUN_DIR/run_identity.tsv"
printf 'edited_sample\t%s\ncomparator_sample\t%s\ncomparison_type\tPEER_COMPARATOR\n' "$EDITED" "$COMPARATOR" >> "$RUN_DIR/run_identity.tsv"
printf 'edited_r1\t%s\nedited_r2\t%s\ncomparator_r1\t%s\ncomparator_r2\t%s\n' "$edited_r1" "$edited_r2" "$starter_r1" "$starter_r2" >> "$RUN_DIR/run_identity.tsv"
printf 'prokadiff_version\t%s\nbreseq_version\t%s\ngdtools_version\t%s\nbowtie2_version\t%s\n' \
    "$($PROKADIFF --version | head -n 1)" "$(breseq --version | head -n 1)" \
    "$(gdtools --version | head -n 1)" "$(bowtie2 --version | head -n 1)" >> "$RUN_DIR/run_identity.tsv"
printf 'slurm_cpus_per_task\t%s\nslurm_mem_per_node\t%s\nslurm_time_limit\t12:00:00\n' \
    "${SLURM_CPUS_PER_TASK:-16}" "${SLURM_MEM_PER_NODE:-56G}" >> "$RUN_DIR/run_identity.tsv"
for guide_id in BL21_G1 BL21_G2; do
    if [[ "$guide_id" == BL21_G1 ]]; then
        spacer="$(manifest_value "$EDITED" spacer_1)"
    else
        spacer="$(manifest_value "$EDITED" spacer_2)"
    fi
    printf 'guide_id\t%s\nspacer\t%s\nproduct_directory\t%s\n' "$guide_id" "$spacer" "$RUN_DIR/$guide_id" >> "$RUN_DIR/run_identity.tsv"
done

printf 'Node: %s\nDate: %s\nProkaDiff Commit: %s\n' "$(hostname)" "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$(git rev-parse HEAD)"
run_status=0
evidence_dir="$RUN_DIR/${EDITED}_evidence"
mkdir -p "$evidence_dir"
if ! /usr/bin/time -v -o "$RUN_DIR/evidence.telemetry.txt" \
    "$PROKADIFF" evidence --ref "$REF" --fastq "$edited_r1" "$edited_r2" \
    --threads 16 --outdir "$evidence_dir"; then
    run_status=1
fi
if ! python3 "$REPO_DIR/benchmark/real_world/scripts/validate_bl21_products.py" \
    --diagnostic-sidecar "$evidence_dir/repeat_ambiguous_junctions.tsv" \
    > "$evidence_dir/diagnostic_sidecar.txt"; then
    echo "repeat-ambiguous sidecar is missing or is not diagnostic-only" >&2
    run_status=1
fi

declare -A oracle_gd
for sample in "$EDITED" "$COMPARATOR"; do
    if [[ "$sample" == "$EDITED" ]]; then
        sample_r1="$edited_r1"
        sample_r2="$edited_r2"
    else
        sample_r1="$starter_r1"
        sample_r2="$starter_r2"
    fi
    breseq_dir="$RUN_DIR/${sample}_breseq"
    oracle_existing="$BL21_DATA/analysis_bl21_cas9_wgs/results/sv/$sample/augmented_full/breseq/output/output.gd"
    if [[ -s "$oracle_existing" ]]; then
        mkdir -p "$breseq_dir/output"
        cp "$oracle_existing" "$breseq_dir/output/output.gd"
        echo "Verified precomputed breseq oracle for $sample from $oracle_existing"
    else
        if ! /usr/bin/time -v -o "$RUN_DIR/${sample}.breseq.telemetry.txt" \
            breseq --reference "$REF_FASTA" --num-processors 16 --output "$breseq_dir" "$sample_r1" "$sample_r2"; then
            run_status=1
        fi
    fi
    oracle_gd["$sample"]="$breseq_dir/output/output.gd"
done

if [[ -s "$evidence_dir/output.gd" && -s "${oracle_gd[$EDITED]}" ]]; then
    bash "$REPO_DIR/benchmark/real_world/scripts/preserve_gdtools_residuals.sh" \
        "$evidence_dir/output.gd" "${oracle_gd[$EDITED]}" "$RUN_DIR/single_sample_residuals" prokadiff_edited_vs_breseq || run_status=1
    if ! python3 "$REPO_DIR/benchmark/real_world/scripts/review_target_locus.py" \
        --sample "$EDITED" --target-id edit_1 --seq-id NZ_CP053602.1 --start 334876 --end 335735 \
        --prokadiff-gd "$evidence_dir/output.gd" --breseq-gd "${oracle_gd[$EDITED]}" \
        --curated-truth "$REPO_DIR/benchmark/real_world/truth/bl21_curated_truth.tsv" \
        --output "$RUN_DIR/target_locus_review.tsv"; then run_status=1; fi
else
    run_status=1
fi

run_peer_guide() {
    local guide_id="$1" spacer product_dir failed=0
    case "$guide_id" in
        BL21_G1) spacer="$(manifest_value "$EDITED" spacer_1)" ;;
        BL21_G2) spacer="$(manifest_value "$EDITED" spacer_2)" ;;
        *) echo "unknown guide $guide_id" >&2; return 1 ;;
    esac
    product_dir="$RUN_DIR/$guide_id"
    mkdir -p "$product_dir"
    if ! /usr/bin/time -v -o "$RUN_DIR/${guide_id}.telemetry.txt" \
        "$PROKADIFF" \
        --starter "$starter_r1" "$starter_r2" \
        --edited "$edited_r1" "$edited_r2" \
        --ref "$REF" --intended "$INTENDED" --editor cas9 \
        --spacer "$spacer" --pam NGG --schema-version v2 \
        --threads 16 --outdir "$product_dir"; then
        return 1
    fi

    cp "$RUN_DIR/run_identity.tsv" "$product_dir/provenance_m5.tsv"
    residual_dir="$product_dir/oracle_residuals"
    bash "$REPO_DIR/benchmark/real_world/scripts/preserve_gdtools_residuals.sh" "$product_dir/edited.gd" "$product_dir/starter.gd" "$residual_dir" prokadiff_pair || failed=1
    if [[ -s "${oracle_gd[$EDITED]}" && -s "${oracle_gd[$COMPARATOR]}" ]]; then
        bash "$REPO_DIR/benchmark/real_world/scripts/preserve_gdtools_residuals.sh" "${oracle_gd[$EDITED]}" "${oracle_gd[$COMPARATOR]}" "$residual_dir" breseq_pair || failed=1
    else
        failed=1
    fi
    if [[ -f "$residual_dir/prokadiff_pair_A_minus_B.gd" && -f "$residual_dir/breseq_pair_A_minus_B.gd" ]]; then
        bash "$REPO_DIR/benchmark/real_world/scripts/preserve_gdtools_residuals.sh" "$residual_dir/prokadiff_pair_A_minus_B.gd" "$residual_dir/breseq_pair_A_minus_B.gd" "$residual_dir" pair_parity || failed=1
    fi

    if ! python3 "$REPO_DIR/benchmark/real_world/scripts/compare_intended_edits.py" \
        --edit-outcomes "$product_dir/edit_outcomes.tsv" --comparison-role PEER_COMPARATOR --tsv \
        > "$product_dir/intended_metrics.tsv"; then failed=1; fi
    if ! python3 "$REPO_DIR/benchmark/real_world/scripts/validate_reports.py" \
        --strict "$product_dir/report.md" > "$product_dir/report_validation.txt"; then failed=1; fi
    if ! python3 "$REPO_DIR/benchmark/real_world/scripts/validate_bl21_products.py" \
        --products "$product_dir" --sample "$EDITED" --guide-id "$guide_id" \
        --comparison-role PEER_COMPARATOR --target-review "$RUN_DIR/target_locus_review.tsv" \
        --diagnostic-sidecar "$evidence_dir/repeat_ambiguous_junctions.tsv" \
        > "$product_dir/product_validation.tsv"; then failed=1; fi
    return "$failed"
}
attempt_peer_guides run_peer_guide
if [[ "$GUIDE_FAILURES" -ne 0 ]]; then
    run_status=1
fi

if [[ -s "$evidence_dir/output.gd" && -s "${oracle_gd[$EDITED]}" && -s "$RUN_DIR/target_locus_review.tsv" \
    && -f "$RUN_DIR/BL21_G1/edit_outcomes.tsv" && -f "$RUN_DIR/BL21_G2/edit_outcomes.tsv" ]]; then
    if ! python3 "$REPO_DIR/benchmark/real_world/scripts/reconcile_target_locus.py" \
        --sample "$EDITED" --seq-id NZ_CP053602.1 --start 334876 --end 335735 \
        --prokadiff-gd "$evidence_dir/output.gd" --breseq-gd "${oracle_gd[$EDITED]}" \
        --curated-truth "$REPO_DIR/benchmark/real_world/truth/bl21_curated_truth.tsv" \
        --products "$RUN_DIR/BL21_G1" "$RUN_DIR/BL21_G2" \
        --output "$RUN_DIR/target_locus_reconciliation.tsv"; then run_status=1; fi
else
    run_status=1
fi

if [[ $run_status -ne 0 ]]; then
    echo "${EDITED} vs ${COMPARATOR} validation failed; both guide attempts were made when prerequisites permitted." >&2
    exit 1
fi
echo "${EDITED} vs ${COMPARATOR} PEER_COMPARATOR validation finished in $RUN_DIR."
