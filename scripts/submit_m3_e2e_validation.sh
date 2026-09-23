#!/usr/bin/env bash
# Milestone 3 E2E Validation Job on Slurm partition qcpu_18i
#
#SBATCH -p qcpu_18i
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -J pd_m3_e2e
#SBATCH -o benchmark/logs/%x-%j.out
#SBATCH -e benchmark/logs/%x-%j.err

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "error: submit this script with sbatch (partition qcpu_18i); do not run it on a login node." >&2
  exit 1
fi

THREADS="${SLURM_CPUS_PER_TASK:-8}"
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"

CONDA_ENV="${PROKDIFF_CONDA_ENV:-/hpcfs/fhome/caizhh/.conda/envs/BactGenome}"
if [[ -f "${CONDA_ENV}/bin/activate" ]]; then
  # shellcheck disable=SC1091
  source "${CONDA_ENV}/bin/activate" "${CONDA_ENV}" 2>/dev/null || true
fi
export PATH="${CONDA_ENV}/bin:${PATH:-}"

mkdir -p "${ROOT}/benchmark/logs" "${ROOT}/benchmark/results" "${ROOT}/testdata/generated"

if ! command -v wgsim >/dev/null 2>&1; then
  echo "error: wgsim not found. Stop layer-1; do not invent results." >&2
  exit 1
fi

PROKDIFF="${ROOT}/target/release/prokadiff"
echo "building release prokadiff on compute node..."
cargo build --release -p prokadiff

echo "=== Generating synthetic dataset synth_cas9_near ==="
bash "${ROOT}/testdata/generate.sh" synth_cas9_near "${ROOT}/testdata/generated/synth_cas9_near"
GEN="${ROOT}/testdata/generated/synth_cas9_near"
JOBOUT="${ROOT}/benchmark/results/m3_e2e_${SLURM_JOB_ID}"
mkdir -p "${JOBOUT}/cas9" "${JOBOUT}/dsb"
SPACER="$(tr -d '[:space:]' < "${GEN}/spacer.txt")"

echo "=== Running ProkaDiff with Cas9 editor and spacer ==="
/usr/bin/time -v -o "${JOBOUT}/cas9.time" \
  "${PROKDIFF}" \
    --starter "${GEN}/starter_R1.fastq" --starter "${GEN}/starter_R2.fastq" \
    --edited "${GEN}/edited_R1.fastq" --edited "${GEN}/edited_R2.fastq" \
    --ref "${GEN}/ref.fa" \
    --editor cas9 \
    --spacer "${SPACER}" \
    --threads "${THREADS}" \
    --outdir "${JOBOUT}/cas9" \
  | tee "${JOBOUT}/cas9.stdout"

echo "=== Running ProkaDiff with DSB editor (no spacer) ==="
/usr/bin/time -v -o "${JOBOUT}/dsb.time" \
  "${PROKDIFF}" \
    --starter "${GEN}/starter_R1.fastq" --starter "${GEN}/starter_R2.fastq" \
    --edited "${GEN}/edited_R1.fastq" --edited "${GEN}/edited_R2.fastq" \
    --ref "${GEN}/ref.fa" \
    --editor dsb \
    --threads "${THREADS}" \
    --outdir "${JOBOUT}/dsb" \
  | tee "${JOBOUT}/dsb.stdout"

echo "=== Running M3 Assertion and Verification Suite ==="
python3 - "${GEN}" "${JOBOUT}" <<'PY'
import sys
import json
from pathlib import Path

gen = Path(sys.argv[1])
job = Path(sys.argv[2])
cas9_dir = job / "cas9"
dsb_dir = job / "dsb"

snp_pos = (gen / "snp.txt").read_text().strip().split("\t")[0]
print(f"Target synthetic SNP position: {snp_pos}")

# Helper to read TSV
def read_tsv(path: Path):
    if not path.exists():
        return []
    lines = path.read_text().splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line.strip()]

# 1. Check unintended.tsv
cas9_unintended = read_tsv(cas9_dir / "unintended.tsv")
dsb_unintended = read_tsv(dsb_dir / "unintended.tsv")

cas9_snps = [r for r in cas9_unintended if r["position"] == snp_pos and r["gd_type"] == "SNP"]
assert cas9_snps, f"FAIL: SNP at {snp_pos} not found in cas9 unintended.tsv"
cas9_snp = cas9_snps[0]
assert cas9_snp["class"] == "near_homolog", f"FAIL: cas9 class is {cas9_snp['class']}, expected near_homolog"
assert cas9_snp["pam_profile"] != "", "FAIL: cas9 pam_profile is empty"
assert cas9_snp["pam_profile"] != "NGG", f"FAIL: cas9 pam_profile has query pattern NGG instead of observed PAM: {cas9_snp['pam_profile']}"
assert int(cas9_snp["distance_to_site"]) <= 50, f"FAIL: distance_to_site {cas9_snp['distance_to_site']} > 50"
print(f"✓ cas9 unintended.tsv verified: class={cas9_snp['class']}, pam_profile={cas9_snp['pam_profile']}, dist={cas9_snp['distance_to_site']}, mm={cas9_snp['offtarget_mismatch']}")

dsb_snps = [r for r in dsb_unintended if r["position"] == snp_pos and r["gd_type"] == "SNP"]
assert dsb_snps, f"FAIL: SNP at {snp_pos} not found in dsb unintended.tsv"
dsb_snp = dsb_snps[0]
assert dsb_snp["class"] == "scattered_snv", f"FAIL: dsb class is {dsb_snp['class']}, expected scattered_snv"
assert dsb_snp["pam_profile"] == "", f"FAIL: dsb pam_profile should be empty, got {dsb_snp['pam_profile']}"
print(f"✓ dsb unintended.tsv verified: class={dsb_snp['class']}")

# 2. Check offtarget_sites.tsv
cas9_sites = read_tsv(cas9_dir / "offtarget_sites.tsv")
assert len(cas9_sites) > 0, "FAIL: cas9 offtarget_sites.tsv has no sites"
for site in cas9_sites:
    assert "SITE_UNKNOWN" not in site["site_id"], f"FAIL: SITE_UNKNOWN found in site: {site}"
    assert site["pam"] != "", f"FAIL: empty PAM in site: {site}"
print(f"✓ cas9 offtarget_sites.tsv verified: {len(cas9_sites)} sites, no SITE_UNKNOWN")

# 3. Check mutation_offtarget_links.tsv
cas9_links = read_tsv(cas9_dir / "mutation_offtarget_links.tsv")
assert len(cas9_links) > 0, "FAIL: cas9 mutation_offtarget_links.tsv has no links"
sites_by_id = {s["site_id"]: s for s in cas9_sites}
for link in cas9_links:
    assert "SITE_UNKNOWN" not in link["site_id"], f"FAIL: SITE_UNKNOWN found in link: {link}"
    assert link["site_id"].startswith("SITE_"), f"FAIL: invalid site_id {link['site_id']}"
    assert link["site_id"] in sites_by_id, f"FAIL: link references unknown site_id {link['site_id']}"
    matching_site = sites_by_id[link["site_id"]]
    assert link["pam"] == matching_site["pam"], f"FAIL: link pam {link['pam']} != site pam {matching_site['pam']}"
    assert link["mismatches"] == matching_site["mismatches"], f"FAIL: link mismatches {link['mismatches']} != site {matching_site['mismatches']}"
    assert int(link["distance_to_site"]) <= 50, f"FAIL: link distance {link['distance_to_site']} > 50"
    assert link["pam"] != "NGG", f"FAIL: link pam has query pattern NGG instead of observed PAM: {link['pam']}"
    # Verify no fabricated distance=0 / mismatch=0 association:
    mut_pos = int(link["mutation_position"])
    site_start = int(matching_site["start"])
    site_end = int(matching_site["end"])
    dist = int(link["distance_to_site"])
    if not (site_start <= mut_pos <= site_end):
        assert dist > 0, f"FAIL: fabricated distance=0 for non-overlapping site {link}"
print(f"✓ cas9 mutation_offtarget_links.tsv verified: {len(cas9_links)} links, valid site_ids, observed PAMs, and consistent site correspondence")

dsb_links = read_tsv(dsb_dir / "mutation_offtarget_links.tsv")
assert len(dsb_links) == 0, f"FAIL: dsb should have 0 links, got {len(dsb_links)}"
print(f"✓ dsb mutation_offtarget_links.tsv verified: 0 links")

# 4. Check TSV headers remain unchanged
assert (cas9_dir / "unintended.tsv").read_text().splitlines()[0] == \
    "seq_id\tposition\tend\tgd_type\tref\talt\tclass\teditor\tpam_profile\tofftarget_mismatch\tdistance_to_site\tside2_seq_id\tside2_position"
assert (cas9_dir / "offtarget_sites.tsv").read_text().splitlines()[0] == \
    "site_id\tseq_id\tstart\tend\tstrand\ttarget_seq\tpam\tmismatches\tbulge_type\tbulge_size\tsearch_backend\tcfd_score\thsu_score"
assert (cas9_dir / "mutation_offtarget_links.tsv").read_text().splitlines()[0] == \
    "mutation_id\tsite_id\tmutation_type\tmutation_position\tsite_start\tsite_end\tdistance_to_site\tmismatches\tpam\tcfd_score\tassociation_window"
print("✓ Output TSV headers verified unchanged")

# 5. Check report.md and provenance.tsv
cas9_report = (cas9_dir / "report.md").read_text()
assert "SITE_UNKNOWN" not in cas9_report, "FAIL: SITE_UNKNOWN found in cas9 report.md"
assert "candidate guide-dependent off-target" in cas9_report.lower(), "FAIL: candidate guide-dependent off-target not in cas9 report.md"
print(f"✓ cas9 report.md verified")

dsb_report = (dsb_dir / "report.md").read_text()
assert "SITE_UNKNOWN" not in dsb_report, "FAIL: SITE_UNKNOWN found in dsb report.md"
print(f"✓ dsb report.md verified")

cas9_prov = (cas9_dir / "provenance.tsv").read_text()
assert "VALIDATED" in cas9_prov or "fixture_validated" in cas9_prov.lower(), f"FAIL: offtarget search status in provenance: {cas9_prov}"
print(f"✓ cas9 provenance.tsv verified")

dsb_prov = (dsb_dir / "provenance.tsv").read_text()
assert "not_requested" in dsb_prov.lower(), f"FAIL: expected not_requested in dsb provenance: {dsb_prov}"
print(f"✓ dsb provenance.tsv verified")

# 6. Global text scan: Ensure SITE_UNKNOWN does NOT appear anywhere in output directories
for path in list(cas9_dir.glob("*")) + list(dsb_dir.glob("*")):
    if path.is_file():
        content = path.read_text(errors="ignore")
        assert "SITE_UNKNOWN" not in content, f"FAIL: SITE_UNKNOWN found in {path}"

print("✓ Global check passed: SITE_UNKNOWN does not appear in any emitted file!")
print("ALL M3 E2E ASSERTIONS PASSED SUCCESSFULLY!")
PY

echo "=== M3 Validation Job Completed Successfully: ${JOBOUT} ==="
