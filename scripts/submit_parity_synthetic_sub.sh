#!/usr/bin/env bash
#SBATCH -p qcpu_18i
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -J pd_sub
#SBATCH -o benchmark/logs/%x-%j.out
#SBATCH -e benchmark/logs/%x-%j.err

set -euo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"

CONDA_ENV="${PROKDIFF_CONDA_ENV:-/hpcfs/fhome/caizhh/.conda/envs/BactGenome}"
export PATH="${CONDA_ENV}/bin:${PATH:-}"

mkdir -p benchmark/logs benchmark/results

FIXTURE_DIR="${ROOT}/testdata/parity/synthetic_sub"
python3 "${FIXTURE_DIR}/generate_sub.py"

OUTDIR="${ROOT}/benchmark/results/synthetic_sub_${SLURM_JOB_ID}"
mkdir -p "${OUTDIR}/breseq_raw"

echo "=== Running breseq oracle on synthetic SUB ==="
breseq -r "${FIXTURE_DIR}/reference.fa" -j 4 -o "${OUTDIR}/breseq_raw" "${FIXTURE_DIR}/reads.fq"

cp "${OUTDIR}/breseq_raw/output/output.gd" "${OUTDIR}/breseq.gd"
cp "${OUTDIR}/breseq.gd" "${FIXTURE_DIR}/oracle.gd"
echo "=== Done breseq on synthetic SUB ==="
