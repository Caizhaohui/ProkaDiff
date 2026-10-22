#!/usr/bin/env bash
#SBATCH --job-name=m6p3c_oracle
#SBATCH --partition=qcpu_18i
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-4%2
#SBATCH --output=benchmark/real_world/results/prjna1088182/m6_p3c/slurm_%A_%a.log
#SBATCH --error=benchmark/real_world/results/prjna1088182/m6_p3c/slurm_%A_%a.err

# ==============================================================================
# ProkaDiff M6-P3c Structural Comparator & Reconciliation Pipeline
#
# Rule 4 & 5: breseq 0.40.2 is invoked strictly as an offline computational
# comparator on Slurm partition qcpu_18i. It is NEVER used as an engine.
#
# Bounded Quintet Cohort:
#   Index 0: SRR28349173 (Cas9_110_111, MG1655 / NC_000913.3)
#   Index 1: SRR28349176 (Cas9_104_105, MG1655 / NC_000913.3)
#   Index 2: SRR28349244 (Cas9_1,       MG1655 / NC_000913.3) [Pilot]
#   Index 3: SRR28349231 (I-SceI_3..6,  BW25113 / CP009273.1)
#   Index 4: SRR28349226 (control_1..4, BW25113 / CP009273.1)
# ==============================================================================

set -euo pipefail
export PATH="/hpcfs/fhome/caizhh/.conda/envs/prokadiff/bin:$PATH"

REPO_DIR="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff"
cd "$REPO_DIR"

EXPECTED_BRESEQ_VERSION="0.40.2"
EXPECTED_BRESEQ_SHA256="e7572467a5f25f0e3a0fb8b423d0ec24ff1901237c275a414c9392c87f8ec202"
EXPECTED_BOWTIE2_VERSION="2.5.4"

# 1. Resolve Target Run Accession
RUNS=(
    "SRR28349173"
    "SRR28349176"
    "SRR28349244"
    "SRR28349231"
    "SRR28349226"
)

if [[ $# -ge 1 && -n "$1" ]]; then
    RUN="$1"
    echo "Running explicit single run target: $RUN"
else
    TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
    if [[ "$TASK_ID" -lt 0 || "$TASK_ID" -ge "${#RUNS[@]}" ]]; then
        echo "Error: SLURM_ARRAY_TASK_ID $TASK_ID out of bounds (0-${#RUNS[@]})" >&2
        exit 2
    fi
    RUN="${RUNS[$TASK_ID]}"
    echo "Running array task $TASK_ID target: $RUN"
fi

BASE_DIR="$REPO_DIR/benchmark/real_world/results/prjna1088182/m6_p3c"
OUT_RUN="$BASE_DIR/$RUN"
BRESEQ_DIR="$OUT_RUN/breseq"
mkdir -p "$BRESEQ_DIR"

echo "=== M6-P3c breseq Oracle Execution ==="
echo "job ${SLURM_JOB_ID:-manual} task ${SLURM_ARRAY_TASK_ID:-none} host $(hostname) run $RUN start $(date -Iseconds)"
echo "commit $(git rev-parse HEAD)"

# 2. Binary and Version Integrity Gate
BRESEQ_BIN="$(which breseq)"
BOWTIE2_BIN="$(which bowtie2)"

actual_breseq_ver="$(breseq --version 2>&1 | head -1 | awk '{print $2}')"
if [[ "$actual_breseq_ver" != "$EXPECTED_BRESEQ_VERSION" ]]; then
    echo "Error: breseq version mismatch: got $actual_breseq_ver, expected $EXPECTED_BRESEQ_VERSION" >&2
    exit 5
fi

actual_breseq_sha="$(sha256sum "$BRESEQ_BIN" | awk '{print $1}')"
if [[ "$actual_breseq_sha" != "$EXPECTED_BRESEQ_SHA256" ]]; then
    echo "Error: breseq binary SHA256 mismatch: got $actual_breseq_sha, expected $EXPECTED_BRESEQ_SHA256" >&2
    exit 5
fi

echo "breseq_version $actual_breseq_ver"
echo "breseq_sha256 $actual_breseq_sha"
echo "bowtie2 $(bowtie2 --version | head -1)"
echo "cpus ${SLURM_CPUS_PER_TASK:-8} mem 32G"

# 3. Storage Preflight Gate (>= 40 GB free)
avail_kb="$(df -k --output=avail "$BASE_DIR" | tail -1 | tr -d ' ')"
if (( avail_kb < 40 * 1024 * 1024 )); then
    echo "Error: BLOCKED_STORAGE: less than 40GB available ($avail_kb KB)" >&2
    exit 6
fi
echo "storage_gate_ok: free space $(awk "BEGIN {print int($avail_kb / 1024 / 1024)}") GB"

# 4. Resolve Reference and Staging Gates
if [[ "$RUN" == "SRR28349231" || "$RUN" == "SRR28349226" ]]; then
    ACCESSION="CP009273.1"
    REF_GBK="$REPO_DIR/data/references/m6_identity/CP009273.1.gbk"
    EXPECTED_REF_FILE_SHA256="3440cc8d28d643dda30e8ed9054656bc093abf61895077678ef4ddad4343df3a"
    EXPECTED_SEQ_SHA256="3761072f8884add5e9facdafc583de5d81a4cf5dbc7e51023ba966f3cf4887ce"
    SAMPLE="BW25113"
else
    ACCESSION="NC_000913.3"
    REF_GBK="$REPO_DIR/data/references/MG1655.gbk"
    EXPECTED_REF_FILE_SHA256="313752bfe10d4aaf04abdb541a1309bce9a6b66f059e75fe215898e736a5a8a5"
    EXPECTED_SEQ_SHA256="e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0"
    SAMPLE="MG1655"
fi

if [[ ! -s "$REF_GBK" ]]; then
    echo "Error: reference GenBank file not found: $REF_GBK" >&2
    exit 4
fi

actual_ref_sha="$(sha256sum "$REF_GBK" | awk '{print $1}')"
if [[ "$actual_ref_sha" != "$EXPECTED_REF_FILE_SHA256" ]]; then
    echo "Error: reference file hash mismatch for $REF_GBK: got $actual_ref_sha, expected $EXPECTED_REF_FILE_SHA256" >&2
    exit 4
fi
echo "reference_gate_ok: $ACCESSION ($REF_GBK)"

# 5. FASTQ Input Verification
PROV_FILE="$REPO_DIR/data/prjna1088182/fastq/$RUN/run_provenance.tsv"
if [[ ! -s "$PROV_FILE" ]]; then
    echo "Error: missing staging provenance for $RUN at $PROV_FILE" >&2
    exit 4
fi

python3 - "$PROV_FILE" "$RUN" "$ACCESSION" "$REPO_DIR" << 'PY'
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
for label, rel, gz_expect in (
    ("r1", row["r1_path"], row["r1_compressed_sha256"]),
    ("r2", row["r2_path"], row["r2_compressed_sha256"]),
):
    path = repo / rel
    if not path.is_file():
        raise SystemExit(f"{run} {label} file not found: {path}")
    gz_digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if gz_digest != gz_expect:
        raise SystemExit(f"{run} {label} gzip sha {gz_digest} != {gz_expect}")
print(f"fastq_gate_ok {run} reads={row['r1_reads']}")
PY

R1_FASTQ="$REPO_DIR/data/prjna1088182/fastq/$RUN/${RUN}_1.fastq.gz"
R2_FASTQ="$REPO_DIR/data/prjna1088182/fastq/$RUN/${RUN}_2.fastq.gz"

# 6. Idempotency Check
if [[ -e "$BRESEQ_DIR/output/output.gd" ]]; then
    echo "Refusing to overwrite existing breseq run at $BRESEQ_DIR/output/output.gd" >&2
    exit 3
fi

# 7. Record Execution Command Metadata
{
    echo "run $RUN"
    echo "sample $SAMPLE"
    echo "reference $ACCESSION"
    echo "reference_gbk $REF_GBK"
    echo "commit $(git rev-parse HEAD)"
    echo "breseq_version $EXPECTED_BRESEQ_VERSION"
    echo "breseq_binary_sha256 $EXPECTED_BRESEQ_SHA256"
    echo "bowtie2_version $(bowtie2 --version | head -1)"
    echo "command breseq -j ${SLURM_CPUS_PER_TASK:-8} -r $REF_GBK -o $BRESEQ_DIR $R1_FASTQ $R2_FASTQ"
    echo "slurm_job_id ${SLURM_JOB_ID:-manual}"
    echo "cpus ${SLURM_CPUS_PER_TASK:-8}"
    echo "requested_memory 32G"
    echo "time_start $(date -Iseconds)"
} > "$BRESEQ_DIR/command.txt"

# 8. breseq Oracle Invocation with Telemetry
echo "breseq execution start: $(date -Iseconds)"
/usr/bin/time -v -o "$BRESEQ_DIR/telemetry.txt" \
    breseq \
    -j "${SLURM_CPUS_PER_TASK:-8}" \
    -r "$REF_GBK" \
    -o "$BRESEQ_DIR" \
    "$R1_FASTQ" \
    "$R2_FASTQ"
echo "breseq execution end: $(date -Iseconds)"

if [[ ! -s "$BRESEQ_DIR/output/output.gd" ]]; then
    echo "Error: breseq completed but output.gd was not created at $BRESEQ_DIR/output/output.gd" >&2
    exit 7
fi

# 9. Automated M6-P3c Structural Reconciliation
PK_GD="$REPO_DIR/benchmark/real_world/results/prjna1088182/m6_p2a/$RUN/output.gd"
if [[ ! -s "$PK_GD" ]]; then
    PK_GD="$REPO_DIR/benchmark/real_world/results/prjna1088182/m6_p1b/$RUN/output.gd"
fi

if [[ -s "$PK_GD" ]]; then
    echo "Running M6-P3c structural reconciliation against ProkaDiff: $PK_GD"
    python3 "$REPO_DIR/benchmark/real_world/scripts/reconcile_m6_p3c.py" \
        --prokadiff-gd "$PK_GD" \
        --comparator-gd "$BRESEQ_DIR/output/output.gd" \
        --reference-gbk "$REF_GBK" \
        --out-tsv "$OUT_RUN/reconciled_structural_events.tsv" \
        --json > "$OUT_RUN/reconciliation_summary.json"
    echo "Reconciliation complete: $OUT_RUN/reconciled_structural_events.tsv"
else
    echo "Notice: ProkaDiff GD not found at $PK_GD; reconciliation skipped."
fi

echo "M6-P3c breseq oracle job for $RUN completed successfully at $(date -Iseconds)"
