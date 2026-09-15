#!/usr/bin/env bash
# scripts/smoke_e2e.sh: Minimal Bowtie2 end-to-end smoke test suitable for GitHub Actions runner

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

if ! command -v bowtie2 >/dev/null 2>&1; then
  if [[ -d "/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin" ]]; then
    export PATH="/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin:${PATH}"
  fi
fi

echo "[smoke_e2e] Checking bowtie2 availability..."
if ! command -v bowtie2 >/dev/null 2>&1; then
  echo "Error: bowtie2 not found on PATH. Please install bowtie2 first." >&2
  exit 1
fi
bowtie2 --version | head -n 1

# Locate or build prokadiff binary
PROKADIFF_BIN="${REPO_ROOT}/target/release/prokadiff"
if [[ ! -x "${PROKADIFF_BIN}" ]]; then
  echo "[smoke_e2e] Building prokadiff release binary..."
  cargo build --release -p prokadiff
fi

TMP_DIR="$(mktemp -d -t prokadiff_smoke_XXXXXX)"
trap 'rm -rf "${TMP_DIR}"' EXIT

echo "[smoke_e2e] Generating synthetic fixture..."
FIXTURE_DIR="${TMP_DIR}/fixture"
mkdir -p "${FIXTURE_DIR}"
python3 "${REPO_ROOT}/testdata/parity/synthetic_sub/generate_sub.py"
cp "${REPO_ROOT}/testdata/parity/synthetic_sub/reference.fa" "${FIXTURE_DIR}/ref.fa"
cp "${REPO_ROOT}/testdata/parity/synthetic_sub/reads.fq" "${FIXTURE_DIR}/reads.fq"

OUTDIR="${TMP_DIR}/output"
mkdir -p "${OUTDIR}"

echo "[smoke_e2e] Running prokadiff evidence engine..."
"${PROKADIFF_BIN}" evidence \
  --ref "${FIXTURE_DIR}/ref.fa" \
  --fastq "${FIXTURE_DIR}/reads.fq" \
  --threads 2 \
  --outdir "${OUTDIR}"

if [[ ! -f "${OUTDIR}/output.gd" ]]; then
  echo "Error: prokadiff evidence failed to produce output.gd" >&2
  exit 1
fi

echo "[smoke_e2e] Verifying output.gd contains expected SUB mutation..."
if ! grep -q "SUB" "${OUTDIR}/output.gd"; then
  echo "Error: output.gd does not contain expected SUB mutation:" >&2
  cat "${OUTDIR}/output.gd" >&2
  exit 1
fi

echo "[smoke_e2e] Output GenomeDiff:"
cat "${OUTDIR}/output.gd"

echo "[smoke_e2e] PASS: Bowtie2 E2E smoke test succeeded."
