#!/usr/bin/env bash
# Layer-1 synth_parent_child: product subtract + intended mask on qcpu_18i.
#   mkdir -p benchmark/logs
#   sbatch scripts/submit_product_synth_parent_child.sh
#
#SBATCH -p qcpu_18i
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -J pd_parent_child
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
if [[ ! -x "${ROOT}/target/release/prokdiff" ]]; then
  echo "building release prokdiff on the compute node..."
  cargo build --release -p prokdiff
fi
PROKDIFF="${ROOT}/target/release/prokdiff"

bash "${ROOT}/testdata/generate.sh" synth_parent_child "${ROOT}/testdata/generated/synth_parent_child"
GEN="${ROOT}/testdata/generated/synth_parent_child"
JOBOUT="${ROOT}/benchmark/results/synth_parent_child_${SLURM_JOB_ID}"
mkdir -p "${JOBOUT}"

/usr/bin/time -v -o "${JOBOUT}/prokdiff.time" \
  "${PROKDIFF}" \
    --starter "${GEN}/starter_R1.fastq" --starter "${GEN}/starter_R2.fastq" \
    --edited "${GEN}/edited_R1.fastq" --edited "${GEN}/edited_R2.fastq" \
    --ref "${GEN}/ref.fa" \
    --editor dsb \
    --intended "${GEN}/intended.tsv" \
    --threads "${THREADS}" \
    --outdir "${JOBOUT}/prokdiff" \
  | tee "${JOBOUT}/prokdiff.stdout"

python3 - "${GEN}" "${JOBOUT}" <<'PY'
import sys
from pathlib import Path

gen = Path(sys.argv[1])
job = Path(sys.argv[2])
tsv = job / "prokdiff" / "unintended.tsv"
if not tsv.is_file():
    print("error: missing", tsv)
    sys.exit(1)

header = tsv.read_text().splitlines()[0].split("\t")
rows = []
for line in tsv.read_text().splitlines()[1:]:
    if not line.strip():
        continue
    rows.append(dict(zip(header, line.split("\t"))))

hist = set()
for line in (gen / "historical_snps.txt").read_text().splitlines()[1:]:
    if line.strip():
        hist.add(line.split("\t")[0])
intended = (gen / "intended.tsv").read_text().splitlines()[1].split("\t")[1]
extra_pos, _ref, extra_alt = (gen / "extra_snp.txt").read_text().strip().split("\t")

positions = {r["position"] for r in rows}
leaked = sorted(hist & positions)
if leaked:
    print("FAIL: historical SNPs in unintended:", leaked)
    sys.exit(2)
if intended in positions:
    print("FAIL: intended SNP still in unintended")
    sys.exit(2)
hits = [r for r in rows if r["position"] == extra_pos and r["gd_type"] == "SNP"]
if not hits:
    print("FAIL: extra SNP", extra_pos, "missing from unintended")
    sys.exit(2)
if hits[0]["class"] != "scattered_snv":
    print("FAIL: extra SNP class", hits[0]["class"], "expected scattered_snv")
    sys.exit(2)
if hits[0]["alt"] != extra_alt:
    print("FAIL: extra SNP alt", hits[0]["alt"], "expected", extra_alt)
    sys.exit(2)
print("synth_parent_child OK: historical subtracted, intended masked, extra SNP scattered_snv")
print(f"unintended_rows={len(rows)}")
PY

echo "done ${JOBOUT}"
