#!/usr/bin/env bash
#SBATCH --job-name=m6p1b_smoke
#SBATCH --partition=qcpu_18i
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --output=benchmark/real_world/results/prjna1088182/m6_p1b/slurm_%j.log
#SBATCH --error=benchmark/real_world/results/prjna1088182/m6_p1b/slurm_%j.err

# Reference-relative evidence only. Not a parent comparison and not M6-P2.
# 8 CPU / 24G / 6h: about 1.5 million paired spots, smaller than BL21 WGS.
set -euo pipefail
export PATH="/hpcfs/fhome/caizhh/.conda/envs/prokadiff/bin:$PATH"
export RUST_BACKTRACE=1

REPO_DIR="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff"
cd "$REPO_DIR"
REF="$REPO_DIR/data/references/MG1655.gbk"
BASE="$REPO_DIR/benchmark/real_world/results/prjna1088182/m6_p1b"
mkdir -p "$BASE"
echo "job ${SLURM_JOB_ID:-manual} host $(hostname) start $(date -Iseconds)"
echo "commit $(git rev-parse HEAD)"
echo "bowtie2 $(bowtie2 --version | head -1)"
echo "cpus ${SLURM_CPUS_PER_TASK:-8} mem 24G"

python3 - << 'PY'
import hashlib, sys
from pathlib import Path
sys.path.insert(0, "benchmark/real_world/scripts")
from m6_reference_policy import MG1655_EXECUTION, RELATIONSHIPS, read_tsv, require_execution_reference
repo = Path(".").resolve()
rows = read_tsv(RELATIONSHIPS)
for run in ("SRR28349244", "SRR28349243"):
    got = require_execution_reference(run, MG1655_EXECUTION, rows)
    if got != "NC_000913.3":
        raise SystemExit(f"{run} reference {got}")
seq = []
on = False
for line in (repo / "data/references/MG1655.gbk").open():
    if line.startswith("ORIGIN"):
        on = True
        continue
    if not on:
        continue
    if line.startswith("//"):
        break
    seq.append("".join(ch for ch in line if ch.isalpha()))
bases = "".join(seq).upper()
digest = hashlib.sha256(bases.encode()).hexdigest()
if len(bases) != 4641652 or digest != "e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0":
    raise SystemExit(f"reference mismatch {len(bases)} {digest}")
print("reference_gate_ok", digest)
PY

verify_fastq() {
    local run="$1" r1_sha="$2" r2_sha="$3" n="$4"
    local dir="$REPO_DIR/data/prjna1088182/fastq/$run"
    local gz1="$dir/${run}_1.fastq.gz" gz2="$dir/${run}_2.fastq.gz"
    local g1 g2 c1 c2
    g1="$(sha256sum "$gz1" | awk '{print $1}')"
    g2="$(sha256sum "$gz2" | awk '{print $1}')"
    [[ "$g1" == "$5" && "$g2" == "$6" ]] || { echo "gzip sha mismatch $run" >&2; exit 4; }
    c1="$(gzip -dc "$gz1" | sha256sum | awk '{print $1}')"
    c2="$(gzip -dc "$gz2" | sha256sum | awk '{print $1}')"
    [[ "$c1" == "$r1_sha" && "$c2" == "$r2_sha" ]] || { echo "content sha mismatch $run" >&2; exit 4; }
    local n1 n2
    n1="$(gzip -dc "$gz1" | awk 'END{printf "%d", NR/4}')"
    n2="$(gzip -dc "$gz2" | awk 'END{printf "%d", NR/4}')"
    [[ "$n1" == "$n" && "$n2" == "$n" ]] || { echo "count mismatch $run $n1 $n2" >&2; exit 4; }
    echo "fastq_gate_ok $run"
}

verify_fastq SRR28349244 \
  e2510f83241b0aa5f497351955de852091feb792e71b75a4f70b867866841e6e \
  edb13a02069cec584d1f7d8ff7a60605f996bb492d599fdaf701cc844323b3ab \
  1600722 \
  df654a7884e4daad915afdd4042b3dab29090226acbabbceed0028c9b9c7890f \
  934e6b76160e60e145fed5ed7eb480f10c1857cdf4212d9a367f2574b120741d
verify_fastq SRR28349243 \
  effe0219b3274dcbd1cb9d227637bbf9b382eec75fb193334a6e43cd6cd31329 \
  a7914f0a3327b46574359e92cd87160b44c56f1775a0810834afc3605de5c681 \
  1514507 \
  ad8eade53611ed73137ba7c154c3b6ab53185eb0f5f380fe2f843111a256977d \
  296634a7eb3607e3d3e0dde128b3ae2e0da4c48a56fb45c40508ab54ca3bb3ce

echo "=== build release from HEAD ==="
cargo build --release -p prokadiff --bin prokadiff
PROKADIFF="$REPO_DIR/target/release/prokadiff"

run_one() {
    local run="$1" out="$2"
    local dir="$REPO_DIR/data/prjna1088182/fastq/$run"
    if [[ -e "$out/output.gd" ]]; then
        echo "refusing to overwrite $out" >&2
        exit 3
    fi
    mkdir -p "$out"
    echo "run $run out $out start $(date -Iseconds)"
    /usr/bin/time -v -o "$out/telemetry.txt" \
        "$PROKADIFF" evidence \
        --ref "$REF" \
        --fastq "$dir/${run}_1.fastq.gz" "$dir/${run}_2.fastq.gz" \
        --threads "${SLURM_CPUS_PER_TASK:-8}" \
        --outdir "$out"
    echo "run $run end $(date -Iseconds)"
}

run_one SRR28349244 "$BASE/SRR28349244"
run_one SRR28349243 "$BASE/SRR28349243"
run_one SRR28349244 "$BASE/SRR28349244_rerun"
echo "smoke jobs finished $(date -Iseconds)"
