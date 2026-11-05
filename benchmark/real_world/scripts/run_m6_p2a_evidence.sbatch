#!/usr/bin/env bash
#SBATCH --job-name=m6p2a_evidence
#SBATCH --partition=qcpu_18i
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=16:00:00
#SBATCH --array=0-22%4
#SBATCH --output=benchmark/real_world/results/prjna1088182/m6_p2a/slurm_%A_%a.log
#SBATCH --error=benchmark/real_world/results/prjna1088182/m6_p2a/slurm_%A_%a.err

# Reference-relative evidence for one eligible paired run.
# 8 CPU / 32G / 16h. P1b used about 1.5 h and 2.9 GB at 1.6 million spots.
# The largest new run has about 5.0 million spots. No --starter. No breseq.
# Uses the P1b release binary. Array tasks do not run cargo: a later
# cargo build --release relinked a different artifact and was reverted.
set -euo pipefail
export PATH="/hpcfs/fhome/caizhh/.conda/envs/prokadiff/bin:$PATH"
export RUST_BACKTRACE=1

REPO_DIR="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff"
cd "$REPO_DIR"
PROKADIFF="$REPO_DIR/target/release/prokadiff"
EXPECTED_BINARY_SHA256="02077c761b1c9bd83cf1b626a3e2e955d8d2a1c63abb3e35852421e0dbc20bfd"
BASE="$REPO_DIR/benchmark/real_world/results/prjna1088182/m6_p2a"

mapfile -t RUNS < <(python3 "$REPO_DIR/benchmark/real_world/scripts/build_m6_p2a_cohort.py" --list-new-runs)
if [[ "${#RUNS[@]}" -ne 23 ]]; then
    echo "expected 23 new paired runs, found ${#RUNS[@]}" >&2
    exit 2
fi
RUN="${RUNS[$SLURM_ARRAY_TASK_ID]}"
OUT="$BASE/$RUN"
mkdir -p "$BASE"
echo "job ${SLURM_JOB_ID:-manual} task ${SLURM_ARRAY_TASK_ID} host $(hostname) run $RUN start $(date -Iseconds)"
echo "commit $(git rev-parse HEAD)"
got_sha="$(sha256sum "$PROKADIFF" | awk '{print $1}')"
if [[ "$got_sha" != "$EXPECTED_BINARY_SHA256" ]]; then
    echo "refusing to run binary $got_sha; P1b binary is $EXPECTED_BINARY_SHA256" >&2
    exit 5
fi
echo "binary_sha256 $got_sha"
echo "bowtie2 $(bowtie2 --version | head -1)"
echo "cpus ${SLURM_CPUS_PER_TASK:-8} mem 32G"

if [[ -e "$OUT/output.gd" ]]; then
    echo "refusing to overwrite $OUT" >&2
    exit 3
fi

gate="$(python3 "$REPO_DIR/benchmark/real_world/scripts/m6_p2a_stage_gate.py" --run "$RUN")"
printf '%s\n' "$gate"
REF="$(printf '%s\n' "$gate" | awk -F '\t' '{print $6}')"
ACCESSION="$(printf '%s\n' "$gate" | awk -F '\t' '{print $3}')"
SAMPLE="$(printf '%s\n' "$gate" | awk -F '\t' '{print $2}')"

prov="$REPO_DIR/data/prjna1088182/fastq/$RUN/run_provenance.tsv"
if [[ ! -s "$prov" ]]; then
    echo "missing staging provenance for $RUN" >&2
    exit 4
fi
python3 - "$prov" "$RUN" "$ACCESSION" "$REPO_DIR" << 'PY'
import csv, hashlib, subprocess, sys
from pathlib import Path
prov, run, accession, repo_s = sys.argv[1:]
repo = Path(repo_s)
with open(prov, encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if len(rows) != 1 or rows[0]["run"] != run:
    raise SystemExit(f"provenance identity failed for {run}")
row = rows[0]
if row["staging_status"] != "READY_FOR_EXECUTION" or row["spot_consistency"] != "MATCH":
    raise SystemExit(f"{run} staging_status={row['staging_status']} consistency={row['spot_consistency']}")
if row["execution_reference"] != accession:
    raise SystemExit(f"{run} provenance reference {row['execution_reference']} != {accession}")
if row["r1_reads"] != row["r2_reads"] or int(row["r1_reads"]) <= 0:
    raise SystemExit(f"{run} read counts are not a nonempty pair")
for label, rel, gz_expect, raw_expect in (
    ("r1", row["r1_path"], row["r1_compressed_sha256"], row["r1_uncompressed_sha256"]),
    ("r2", row["r2_path"], row["r2_compressed_sha256"], row["r2_uncompressed_sha256"]),
):
    path = repo / rel
    gz_digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if gz_digest != gz_expect:
        raise SystemExit(f"{run} {label} gzip sha {gz_digest} != {gz_expect}")
    proc = subprocess.Popen(["gzip", "-dc", str(path)], stdout=subprocess.PIPE)
    assert proc.stdout is not None
    raw = hashlib.sha256()
    count = 0
    for line in proc.stdout:
        raw.update(line)
        count += 1
    if proc.wait() != 0:
        raise SystemExit(f"{run} {label} gzip decode failed")
    if raw.hexdigest() != raw_expect:
        raise SystemExit(f"{run} {label} content sha mismatch")
    if count % 4 != 0 or count // 4 != int(row["r1_reads"]):
        raise SystemExit(f"{run} {label} count {count // 4} != {row['r1_reads']}")
print(f"fastq_gate_ok {run} reads={row['r1_reads']}")
PY

mkdir -p "$OUT"
{
    echo "sample $SAMPLE"
    echo "run $RUN"
    echo "reference $ACCESSION"
    echo "reference_path $REF"
    echo "commit $(git rev-parse HEAD)"
    echo "binary_sha256 $(sha256sum "$PROKADIFF" | awk '{print $1}')"
    echo "bowtie2 $(bowtie2 --version | head -1)"
    echo "command $PROKADIFF evidence --ref $REF --fastq $REPO_DIR/data/prjna1088182/fastq/$RUN/${RUN}_1.fastq.gz $REPO_DIR/data/prjna1088182/fastq/$RUN/${RUN}_2.fastq.gz --threads ${SLURM_CPUS_PER_TASK:-8} --outdir $OUT"
    echo "slurm_job_id ${SLURM_JOB_ID:-manual}"
    echo "cpus ${SLURM_CPUS_PER_TASK:-8}"
    echo "requested_memory 32G"
} > "$OUT/command.txt"

echo "run $RUN start $(date -Iseconds)"
/usr/bin/time -v -o "$OUT/telemetry.txt" \
    "$PROKADIFF" evidence \
    --ref "$REF" \
    --fastq "$REPO_DIR/data/prjna1088182/fastq/$RUN/${RUN}_1.fastq.gz" \
            "$REPO_DIR/data/prjna1088182/fastq/$RUN/${RUN}_2.fastq.gz" \
    --threads "${SLURM_CPUS_PER_TASK:-8}" \
    --outdir "$OUT"
echo "run $RUN end $(date -Iseconds) exit 0"
