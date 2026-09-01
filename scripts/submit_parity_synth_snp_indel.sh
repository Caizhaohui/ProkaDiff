#!/usr/bin/env bash
# Layer-1 synth_snp_indel: ProkDiff evidence vs breseq on qcpu_18i.
# Do not run on a login node. Submit with:
#   mkdir -p benchmark/logs
#   sbatch scripts/submit_parity_synth_snp_indel.sh
#
#SBATCH -p qcpu_18i
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -J pd_snp_indel
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

missing=()
for bin in wgsim bowtie2 bowtie2-build breseq gdtools python3; do
  if ! command -v "${bin}" >/dev/null 2>&1; then
    missing+=("${bin}")
  fi
done
if (( ${#missing[@]} )); then
  echo "error: missing binaries: ${missing[*]}. Stop layer-1; do not invent parity numbers." >&2
  exit 1
fi

echo "building release prokdiff (incremental; no-op if fresh)..."
cargo build --release -p prokdiff
PROKDIFF="${ROOT}/target/release/prokdiff"

bash "${ROOT}/testdata/generate.sh" "${ROOT}/testdata/generated/synth_snp_indel"
GEN="${ROOT}/testdata/generated/synth_snp_indel"
JOBOUT="${ROOT}/benchmark/results/synth_snp_indel_${SLURM_JOB_ID}"
mkdir -p "${JOBOUT}/prokdiff" "${JOBOUT}/breseq"

TIME=(/usr/bin/time -v)

echo "=== ProkDiff evidence ==="
"${TIME[@]}" -o "${JOBOUT}/prokdiff.time" \
  "${PROKDIFF}" evidence \
    --ref "${GEN}/ref.fa" \
    --fastq "${GEN}/edited_R1.fastq" \
    --fastq "${GEN}/edited_R2.fastq" \
    --threads "${THREADS}" \
    --outdir "${JOBOUT}/prokdiff" \
    --keep-bam \
  | tee "${JOBOUT}/prokdiff.stdout"

echo "=== breseq oracle ==="
"${TIME[@]}" -o "${JOBOUT}/breseq.time" \
  breseq -j "${THREADS}" -o "${JOBOUT}/breseq" \
    -r "${GEN}/ref.fa" \
    "${GEN}/edited_R1.fastq" \
    "${GEN}/edited_R2.fastq" \
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

# gdtools SUBTRACT writes "output.gd" in cwd; run inside the job dir.
(
  cd "${JOBOUT}"
  gdtools SUBTRACT "${RUST_GD}" "${BRESEQ_GD}"
  mv -f output.gd subtract_rust_minus_breseq.gd
  gdtools SUBTRACT "${BRESEQ_GD}" "${RUST_GD}"
  mv -f output.gd subtract_breseq_minus_rust.gd
)

python3 - "${JOBOUT}" "${GEN}/truth.gd" "${THREADS}" <<'PY'
import csv
import os
import re
import socket
import subprocess
import sys
from pathlib import Path

jobout = Path(sys.argv[1])
truth = Path(sys.argv[2])
threads = sys.argv[3]


def mut_keys(path: Path):
    keys = []
    if not path.exists():
        return keys
    for line in path.read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("="):
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        kind = parts[0]
        if kind not in {"SNP", "INS", "DEL"}:
            continue
        # skip evidence-only; mutation rows have id then parent then fields
        keys.append((kind, tuple(parts[3:6] if len(parts) >= 6 else parts[1:])))
    return keys


def leftover(path: Path):
    return mut_keys(path)


def parse_time(path: Path):
    wall_s = ""
    rss_kb = ""
    text = path.read_text(errors="replace")
    for line in text.splitlines():
        if "Elapsed (wall clock) time" in line:
            # e.g. 0:00:12 or 12.34
            m = re.search(r":\s*(\d+):(\d+):(\d+(?:\.\d+)?)$", line)
            if m:
                wall_s = str(int(m.group(1)) * 3600 + int(m.group(2)) * 60 + float(m.group(3)))
            else:
                m = re.search(r":\s*(\d+):(\d+(?:\.\d+)?)$", line)
                if m:
                    wall_s = str(int(m.group(1)) * 60 + float(m.group(2)))
        if "Maximum resident set size" in line:
            m = re.search(r":\s*(\d+)", line)
            if m:
                rss_kb = m.group(1)
    return wall_s, rss_kb


def ver(cmd):
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True, timeout=30)
        return out.strip().splitlines()[0][:120]
    except Exception as e:
        return f"unavailable:{e}"


node = socket.gethostname()
breseq_v = ver(["breseq", "--version"])
bt2_v = ver(["bowtie2", "--version"])
over = leftover(jobout / "subtract_rust_minus_breseq.gd")
under = leftover(jobout / "subtract_breseq_minus_rust.gd")
notes = (
    f"over_snp_ins_del={len(over)};under_snp_ins_del={len(under)};"
    f"over={over[:8]!r};under={under[:8]!r}"
)

csv_path = jobout.parent / "synth_snp_indel.csv"
new_file = not csv_path.exists()
with csv_path.open("a", newline="") as fh:
    w = csv.writer(fh)
    if new_file:
        w.writerow(
            [
                "workload",
                "tool",
                "threads",
                "wall_s",
                "max_rss_kb",
                "node",
                "job_id",
                "breseq_version",
                "bowtie2_version",
                "notes",
            ]
        )
    wall, rss = parse_time(jobout / "prokdiff.time")
    w.writerow(
        [
            "synth_snp_indel",
            "prokdiff",
            threads,
            wall,
            rss,
            node,
            os.environ.get("SLURM_JOB_ID", ""),
            breseq_v,
            bt2_v,
            notes,
        ]
    )
    wall, rss = parse_time(jobout / "breseq.time")
    w.writerow(
        [
            "synth_snp_indel",
            "breseq",
            threads,
            wall,
            rss,
            node,
            os.environ.get("SLURM_JOB_ID", ""),
            breseq_v,
            bt2_v,
            notes,
        ]
    )

print(f"wrote {csv_path}")
print(notes)
if over or under:
    print("WARNING: SNP/INS/DEL SUBTRACT leftover is non-empty (see notes).")
    sys.exit(2)
print("SNP/INS/DEL bidirectional SUBTRACT empty.")
PY

echo "done ${JOBOUT}"
