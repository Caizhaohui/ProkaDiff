#!/usr/bin/env bash
set -euo pipefail

repo_root=$(git rev-parse --show-toplevel)
cd "$repo_root"

require_file() {
    test -f "$1" || {
        printf 'missing required validation asset: %s\n' "$1" >&2
        exit 1
    }
}

require_file benchmark/real_world/VALIDATION_REGISTRY.tsv
require_file benchmark/real_world/manifests/prjna1088182_cohort.tsv
require_file benchmark/real_world/manifests/prjna1088182_relationships.tsv
require_file benchmark/real_world/manifests/prjna1088182_truth_coordinate_audit.tsv
require_file benchmark/real_world/manifests/prjna1088182_sra_runinfo.csv
require_file benchmark/real_world/truth/prjna1088182_publication_truth.tsv
require_file benchmark/real_world/manifests/prjna884016.tsv
require_file benchmark/real_world/manifests/prjna884016_sra_runinfo.csv
require_file benchmark/real_world/truth/prjna884016/SY22A1_truth.gd
require_file benchmark/real_world/truth/prjna884016/mob_truth.tsv
require_file benchmark/real_world/truth/prjna884016/jc_truth.tsv
require_file benchmark/real_world/truth/prjna884016/SY22A1_vs_SY221.rdiff

awk -F '\t' 'NR > 1 && $1 == "prjna1088182" {
    if ($4 != "false" || $5 != "false" || $6 != "false" || $7 != "false" || $8 != "false" || $9 != "PENDING") exit 1
    found = 1
} END { exit found ? 0 : 1 }' benchmark/real_world/VALIDATION_REGISTRY.tsv

awk -F '\t' 'NR > 1 && $1 == "prjna884016" {
    if ($4 != "false" || $5 != "false" || $6 != "false" || $7 != "false" || $8 != "false" || $9 != "PENDING") exit 1
    found = 1
} END { exit found ? 0 : 1 }' benchmark/real_world/VALIDATION_REGISTRY.tsv

test "$(awk 'END { print NR - 1 }' benchmark/real_world/manifests/prjna1088182_cohort.tsv)" -eq 25
test "$(awk 'END { print NR - 1 }' benchmark/real_world/truth/prjna1088182_publication_truth.tsv)" -eq 28
test "$(awk '$1 == "MOB" { count += 1 } END { print count + 0 }' benchmark/real_world/truth/prjna884016/SY22A1_truth.gd)" -eq 0
test "$(awk '$1 == "JC" { count += 1 } END { print count + 0 }' benchmark/real_world/truth/prjna884016/SY22A1_truth.gd)" -eq 0
test "$(awk 'NR > 1 { count += 1 } END { print count + 0 }' benchmark/real_world/truth/prjna884016/mob_truth.tsv)" -eq 0
test "$(awk 'NR > 1 { count += 1 } END { print count + 0 }' benchmark/real_world/truth/prjna884016/jc_truth.tsv)" -eq 0

# M7-P0b: sidecar outputs of the invalidated truncated SY22A1_vs_SY221 MUMmer run
# must not live in the active truth directory. Only the zeroed .rdiff remains.
for sidecar in delta report 1delta mdelta 1coords mcoords qdiff snps; do
    stale_sidecar="benchmark/real_world/truth/prjna884016/SY22A1_vs_SY221.${sidecar}"
    if test -e "$stale_sidecar"; then
        printf 'stale MUMmer sidecar present in active truth dir: %s\n' "$stale_sidecar" >&2
        exit 1
    fi
done

while IFS= read -r asset; do
    if git check-ignore -q "$asset"; then
        printf 'validation asset is ignored: %s\n' "$asset" >&2
        exit 1
    fi
    git ls-files --error-unmatch "$asset" >/dev/null || exit 1
done <<'ASSETS'
benchmark/real_world/manifests/prjna1088182_cohort.tsv
benchmark/real_world/manifests/prjna1088182_relationships.tsv
benchmark/real_world/manifests/prjna1088182_truth_coordinate_audit.tsv
benchmark/real_world/manifests/prjna1088182_sra_runinfo.csv
benchmark/real_world/truth/prjna1088182_publication_truth.tsv
benchmark/real_world/manifests/prjna884016_sra_runinfo.csv
benchmark/real_world/truth/prjna884016/SY22A1_truth.gd
benchmark/real_world/truth/prjna884016/mob_truth.tsv
benchmark/real_world/truth/prjna884016/jc_truth.tsv
benchmark/real_world/truth/prjna884016/SY22A1_vs_SY221.rdiff
ASSETS

python3 benchmark/real_world/scripts/m6_reference_policy.py --check

printf 'validation asset checks passed\n'
