#!/usr/bin/env bash
# Layer-2 Example 6 (A. baylyi cassette / JC): ProkDiff evidence vs official workshop
# output.gd AND same-node breseq. Mode = Example 6b (both refs as -r; ProkDiff has no
# --junction-only-reference). Submit with:
#   mkdir -p benchmark/logs testdata/layer2
#   sbatch scripts/submit_parity_example6.sh
#
# FASTQ is not posted next to the wiki results; this job tries known URLs and
# Advanced_Workshop_Files.zip. If reads are missing, exit non-zero and do not invent numbers.
#
#SBATCH -p qcpu_18i
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH -t 24:00:00
#SBATCH -J pd_ex6
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

DATA="${ROOT}/testdata/layer2/example6"
mkdir -p "${ROOT}/benchmark/logs" "${ROOT}/benchmark/results" "${DATA}"

missing=()
for bin in bowtie2 bowtie2-build breseq gdtools python3; do
  if ! command -v "${bin}" >/dev/null 2>&1; then
    missing+=("${bin}")
  fi
done
if ! command -v curl >/dev/null 2>&1 && ! command -v wget >/dev/null 2>&1; then
  missing+=("curl-or-wget")
fi
if (( ${#missing[@]} )); then
  echo "error: missing binaries: ${missing[*]}. Stop layer-2 example6; do not invent leftover / wall / RSS." >&2
  exit 1
fi
if [[ ! -x /usr/bin/time ]]; then
  echo "error: /usr/bin/time not found. Stop layer-2; do not invent timing numbers." >&2
  exit 1
fi

fetch_url() {
  local url="$1"
  local dest="$2"
  if command -v curl >/dev/null 2>&1; then
    curl -fL --retry 2 --retry-delay 3 --connect-timeout 30 -o "${dest}.partial" "${url}"
  else
    wget -O "${dest}.partial" "${url}"
  fi
}

# Try URL; 404 / no network is non-fatal here so later candidates can run.
try_download() {
  local url="$1"
  local dest="$2"
  if [[ -s "${dest}" ]]; then
    echo "skip download (exists): ${dest}"
    return 0
  fi
  echo "try ${url}"
  if fetch_url "${url}" "${dest}"; then
    mv "${dest}.partial" "${dest}"
    return 0
  fi
  rm -f "${dest}.partial"
  echo "not available: ${url}"
  return 1
}

require_download() {
  local url="$1"
  local dest="$2"
  if try_download "${url}" "${dest}"; then
    return 0
  fi
  echo "error: required download failed: ${url}" >&2
  echo "Stop layer-2 example6; do not invent leftover / wall / RSS." >&2
  exit 1
}

WIKI="https://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing"
GH="https://raw.githubusercontent.com/barricklab/Abaylyi-EE/master/reference"

require_download \
  "${GH}/Acinetobacter-baylyi-ADP1-WT.gff3" \
  "${DATA}/Acinetobacter-baylyi-ADP1-WT.gff3"
require_download \
  "${GH}/pBTK622_tdk-kanR_cassette_for%20Golden_Transformation.gbk" \
  "${DATA}/pBTK622_tdk-kanR_cassette_for_Golden_Transformation.gbk"
require_download \
  "${WIKI}/IntroWorkshop/Abaylyi_with_cassette_normal/output.gd" \
  "${DATA}/official_example6b.gd"

# Optional workshop zip (curated GDs; too small for genome FASTQ, but search anyway).
try_download "${WIKI}/Advanced_Workshop_Files.zip" "${DATA}/Advanced_Workshop_Files.zip" || true
if [[ -s "${DATA}/Advanced_Workshop_Files.zip" && ! -d "${DATA}/Advanced_Workshop_Files" ]]; then
  mkdir -p "${DATA}/Advanced_Workshop_Files"
  python3 - "${DATA}/Advanced_Workshop_Files.zip" "${DATA}/Advanced_Workshop_Files" <<'PY'
import sys, zipfile
zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])
PY
fi

R1="${DATA}/G2_CCGTCC_L007_R1_001.fastq.gz"
R2="${DATA}/G2_CCGTCC_L007_R2_001.fastq.gz"
FQ_CANDIDATES=(
  "${WIKI}/IntroWorkshop/G2_CCGTCC_L007_R1_001.fastq.gz"
  "${WIKI}/G2_CCGTCC_L007_R1_001.fastq.gz"
  "${WIKI}/IntroWorkshop/data/G2_CCGTCC_L007_R1_001.fastq.gz"
)
if [[ ! -s "${R1}" ]]; then
  for url in "${FQ_CANDIDATES[@]}"; do
    if try_download "${url}" "${R1}"; then
      r2url="${url/R1_001/R2_001}"
      try_download "${r2url}" "${R2}" || true
      break
    fi
  done
fi
if [[ -s "${R1}" && ! -s "${R2}" ]]; then
  try_download "${WIKI}/IntroWorkshop/G2_CCGTCC_L007_R2_001.fastq.gz" "${R2}" || true
  try_download "${WIKI}/G2_CCGTCC_L007_R2_001.fastq.gz" "${R2}" || true
fi

# Copies extracted from the zip, if any.
if [[ ! -s "${R1}" ]]; then
  found="$(find "${DATA}" -type f \( -name '*G2_CCGTCC*R1*' -o -name '*_R1_001.fastq.gz' \) | head -n 1 || true)"
  if [[ -n "${found}" ]]; then
    ln -sfn "${found}" "${R1}"
  fi
fi
if [[ ! -s "${R2}" ]]; then
  found="$(find "${DATA}" -type f \( -name '*G2_CCGTCC*R2*' -o -name '*_R2_001.fastq.gz' \) | head -n 1 || true)"
  if [[ -n "${found}" ]]; then
    ln -sfn "${found}" "${R2}"
  fi
fi

if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
  echo "error: Example 6 FASTQ not available after trying wiki URLs and Advanced_Workshop_Files.zip." >&2
  echo "Tried R1 names: G2_CCGTCC_L007_R1_001.fastq.gz under ${WIKI} and ${WIKI}/IntroWorkshop/" >&2
  echo "Official 6b output.gd is present at ${DATA}/official_example6b.gd but the engine cannot run without reads." >&2
  echo "Stop layer-2 example6; do not invent leftover / wall / RSS." >&2
  exit 1
fi

# GFF3 has ##FASTA (seq id ADP1-WT). ProkDiff evidence reads FASTA/GBK, not GFF3.
python3 - "${DATA}/Acinetobacter-baylyi-ADP1-WT.gff3" "${DATA}/ADP1-WT.fa" <<'PY'
import sys
from pathlib import Path
gff, dest = Path(sys.argv[1]), Path(sys.argv[2])
text = gff.read_text()
marker = "##FASTA"
idx = text.find(marker)
if idx < 0:
    raise SystemExit(f"error: {gff} has no ##FASTA; cannot build ProkDiff reference. Do not invent results.")
dest.write_text(text[idx + len(marker) :].lstrip())
print(f"wrote {dest}")
PY

python3 - "${DATA}/pBTK622_tdk-kanR_cassette_for_Golden_Transformation.gbk" "${DATA}/cassette.fa" <<'PY'
import sys
from pathlib import Path
gbk, dest = Path(sys.argv[1]), Path(sys.argv[2])
name = "pYTK001_with_tdk_kan"
seq = []
in_origin = False
for line in gbk.read_text().splitlines():
    t = line.strip()
    if t.upper().startswith("LOCUS") and len(t.split()) > 1:
        raw = t.split()[1]
        name = raw.replace("/", "_")
    if t.upper().startswith("ORIGIN"):
        in_origin = True
        continue
    if in_origin:
        if t == "//":
            break
        seq.extend(ch.upper() for ch in t if ch.isalpha())
if not seq:
    raise SystemExit(f"error: empty ORIGIN in {gbk}")
chunks = ["".join(seq[i : i + 80]) for i in range(0, len(seq), 80)]
dest.write_text(">" + name + "\n" + "\n".join(chunks) + "\n")
print(f"wrote {dest} name={name} len={len(seq)}")
PY

if [[ ! -x "${ROOT}/target/release/prokdiff" ]]; then
  echo "building release prokdiff on the compute node..."
  cargo build --release -p prokdiff
fi
PROKDIFF="${ROOT}/target/release/prokdiff"

JOBOUT="${ROOT}/benchmark/results/example6_${SLURM_JOB_ID}"
mkdir -p "${JOBOUT}/prokdiff" "${JOBOUT}/breseq"
TIME=(/usr/bin/time -v)
ORACLE_GD="${DATA}/official_example6b.gd"
GFF="${DATA}/Acinetobacter-baylyi-ADP1-WT.gff3"
CAS_GBK="${DATA}/pBTK622_tdk-kanR_cassette_for_Golden_Transformation.gbk"

echo "=== ProkDiff evidence (example6b: both refs) ==="
"${TIME[@]}" -o "${JOBOUT}/prokdiff.time" \
  "${PROKDIFF}" evidence \
    --ref "${DATA}/ADP1-WT.fa" \
    --ref "${DATA}/cassette.fa" \
    --fastq "${R1}" \
    --fastq "${R2}" \
    --threads "${THREADS}" \
    --outdir "${JOBOUT}/prokdiff" \
    --keep-bam \
  | tee "${JOBOUT}/prokdiff.stdout"

echo "=== breseq oracle example 6b (same node; env is breseq 0.40.2) ==="
"${TIME[@]}" -o "${JOBOUT}/breseq.time" \
  breseq -j "${THREADS}" -o "${JOBOUT}/breseq" \
    -r "${GFF}" \
    -r "${CAS_GBK}" \
    "${R1}" \
    "${R2}" \
  | tee "${JOBOUT}/breseq.stdout"

RUST_GD="${JOBOUT}/prokdiff/output.gd"
BRESEQ_GD="${JOBOUT}/breseq/output/output.gd"
if [[ ! -f "${RUST_GD}" ]]; then
  echo "error: ProkDiff did not write ${RUST_GD}" >&2
  exit 1
fi
if [[ ! -f "${BRESEQ_GD}" ]]; then
  echo "error: breseq did not write ${BRESEQ_GD}" >&2
  exit 1
fi

(
  cd "${JOBOUT}"
  gdtools SUBTRACT "${RUST_GD}" "${BRESEQ_GD}"
  mv -f output.gd subtract_rust_minus_breseq.gd
  gdtools SUBTRACT "${BRESEQ_GD}" "${RUST_GD}"
  mv -f output.gd subtract_breseq_minus_rust.gd
  gdtools SUBTRACT "${RUST_GD}" "${ORACLE_GD}"
  mv -f output.gd subtract_rust_minus_official.gd
  gdtools SUBTRACT "${ORACLE_GD}" "${RUST_GD}"
  mv -f output.gd subtract_official_minus_rust.gd
)

python3 "${ROOT}/scripts/layer2_gd_compare.py" \
  --jobout "${JOBOUT}" \
  --workload example6 \
  --threads "${THREADS}" \
  --prokdiff-gd "${RUST_GD}" \
  --breseq-gd "${BRESEQ_GD}" \
  --oracle-gd "${ORACLE_GD}" \
  --csv "${ROOT}/benchmark/results/example6.csv" \
  --jc-tol-bp 5

echo "done ${JOBOUT}"
