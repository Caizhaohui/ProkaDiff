#!/usr/bin/env bash
# Layer-1 synth_cas9_near: cas9 → near_homolog, dsb → scattered_snv.
#   mkdir -p benchmark/logs
#   sbatch scripts/submit_product_synth_cas9_near.sh
#
#SBATCH -p qcpu_18i
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -J pd_cas9_near
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

bash "${ROOT}/testdata/generate.sh" synth_cas9_near "${ROOT}/testdata/generated/synth_cas9_near"
GEN="${ROOT}/testdata/generated/synth_cas9_near"
JOBOUT="${ROOT}/benchmark/results/synth_cas9_near_${SLURM_JOB_ID}"
mkdir -p "${JOBOUT}/cas9" "${JOBOUT}/dsb"
SPACER="$(tr -d '[:space:]' < "${GEN}/spacer.txt")"

echo "=== cas9 ==="
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

echo "=== dsb ==="
/usr/bin/time -v -o "${JOBOUT}/dsb.time" \
  "${PROKDIFF}" \
    --starter "${GEN}/starter_R1.fastq" --starter "${GEN}/starter_R2.fastq" \
    --edited "${GEN}/edited_R1.fastq" --edited "${GEN}/edited_R2.fastq" \
    --ref "${GEN}/ref.fa" \
    --editor dsb \
    --threads "${THREADS}" \
    --outdir "${JOBOUT}/dsb" \
  | tee "${JOBOUT}/dsb.stdout"

python3 - "${GEN}" "${JOBOUT}" <<'PY'
import sys
from pathlib import Path

gen = Path(sys.argv[1])
job = Path(sys.argv[2])
snp_pos = (gen / "snp.txt").read_text().strip().split("\t")[0]


def snp_row(path: Path):
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:] if line.strip()]
    hits = [r for r in rows if r["position"] == snp_pos and r["gd_type"] == "SNP"]
    if not hits:
        raise SystemExit(f"FAIL: SNP {snp_pos} missing in {path}")
    return hits[0]


cas9 = snp_row(job / "cas9" / "unintended.tsv")
dsb = snp_row(job / "dsb" / "unintended.tsv")
if cas9["class"] != "near_homolog":
    raise SystemExit(f"FAIL: cas9 class {cas9['class']} expected near_homolog")
if dsb["class"] != "scattered_snv":
    raise SystemExit(f"FAIL: dsb class {dsb['class']} expected scattered_snv")
print("synth_cas9_near OK: cas9=near_homolog dsb=scattered_snv")
PY

echo "done ${JOBOUT}"
