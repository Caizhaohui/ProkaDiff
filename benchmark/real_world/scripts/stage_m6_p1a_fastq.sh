#!/usr/bin/env bash
# Stage one authorized M6-P1a smoke SRA run to paired FASTQ.
# Does not run ProkaDiff, Bowtie2, or breseq.
set -euo pipefail

RUN="${1:?SRA run accession required}"
REPO_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_DIR"

SRA_BIN="$REPO_DIR/data/tools/sra-tools/bin"
export PATH="$SRA_BIN:$PATH"

DEST="$REPO_DIR/data/prjna1088182/fastq/$RUN"
TMP="$REPO_DIR/data/prjna1088182/tmp/$RUN"
mkdir -p "$DEST" "$TMP/fqd"
PROV="$DEST/run_provenance.tsv"
LOG="$DEST/staging.log"
START_ISO="$(date -Iseconds)"

write_blocked() {
    local status="$1"
    local note="$2"
    {
        printf 'run\t%s\n' "$RUN"
        printf 'staging_status\t%s\n' "$status"
        printf 'notes\t%s\n' "$note"
        printf 'start_time\t%s\n' "$START_ISO"
        printf 'end_time\t%s\n' "$(date -Iseconds)"
        printf 'hostname\t%s\n' "$(hostname)"
        printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-manual}"
        printf 'exit_code\t%s\n' "1"
    } > "$PROV"
}

fail() {
    local status="$1"
    local note="$2"
    echo "BLOCKED $status $note" | tee -a "$LOG" >&2
    write_blocked "$status" "$note"
    exit 1
}

{
    echo "=== M6-P1a stage $RUN ==="
    echo "start $START_ISO"
    echo "host $(hostname)"
    echo "job ${SLURM_JOB_ID:-manual}"
    echo "prefetch $($SRA_BIN/prefetch --version 2>&1 | head -1)"
    echo "fasterq-dump $($SRA_BIN/fasterq-dump --version 2>&1 | head -1)"
    echo "vdb-validate $($SRA_BIN/vdb-validate --version 2>&1 | head -1)"
    echo "pigz $(pigz --version 2>&1 | head -1)"
} | tee "$LOG"

if ! python3 - "$RUN" "$REPO_DIR" << 'PY' | tee -a "$LOG"
import hashlib, sys
from pathlib import Path
sys.path.insert(0, sys.argv[2] + "/benchmark/real_world/scripts")
from m6_reference_policy import (
    MG1655_EXECUTION, RELATIONSHIPS, assert_authorized_smoke_run,
    read_tsv, require_execution_reference,
)
run = sys.argv[1]
assert_authorized_smoke_run(run)
rows = read_tsv(RELATIONSHIPS)
require_execution_reference(run, MG1655_EXECUTION, rows)
gbk = Path(sys.argv[2]) / "data/references/MG1655.gbk"
seq = []
on = False
for line in gbk.open():
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
expected = "e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0"
if len(bases) != 4641652 or digest != expected:
    raise SystemExit(f"reference hash mismatch len={len(bases)} sha={digest}")
print(f"reference_gate NC_000913.3 len=4641652 sha256={digest}")
PY
then
    fail BLOCKED_REFERENCE "reference policy or local NC_000913.3 sequence hash failed"
fi

avail="$(df -B1 --output=avail "$DEST" | tail -1 | tr -d ' ')"
need=$((40 * 1024 * 1024 * 1024))
if [[ "$avail" -lt "$need" ]]; then
    fail BLOCKED_STORAGE "free bytes $avail below 40GiB"
fi

# max-size is in KB. 5_000_000 KB is about 4.8 GiB, above the ~155 MB SRA objects.
if ! prefetch --max-size 5000000 --output-directory "$TMP" "$RUN" >> "$LOG" 2>&1; then
    fail BLOCKED_DOWNLOAD "prefetch failed"
fi

sra_file="$(find "$TMP" -type f -name "${RUN}.sra" -o -type f -name "${RUN}.sralite" | head -1)"
if [[ -z "$sra_file" ]]; then
    fail BLOCKED_DOWNLOAD "prefetch produced no ${RUN}.sra"
fi
if ! vdb-validate "$sra_file" >> "$LOG" 2>&1; then
    fail BLOCKED_SRA_VALIDATION "vdb-validate failed for $sra_file"
fi

threads="${SLURM_CPUS_PER_TASK:-4}"
if ! fasterq-dump "$sra_file" --split-files --threads "$threads" \
    --outdir "$DEST" --temp "$TMP/fqd" >> "$LOG" 2>&1; then
    fail BLOCKED_DOWNLOAD "fasterq-dump failed"
fi

shopt -s nullglob
extra=()
for path in "$DEST"/*.fastq "$DEST"/*.fastq.gz; do
    base="$(basename "$path")"
    case "$base" in
        "${RUN}_1.fastq"|"${RUN}_2.fastq") ;;
        *) extra+=("$base") ;;
    esac
done
if [[ ${#extra[@]} -gt 0 ]]; then
    fail BLOCKED_PAIRING "unexpected FASTQ retained and not discarded: ${extra[*]}"
fi
r1="$DEST/${RUN}_1.fastq"
r2="$DEST/${RUN}_2.fastq"
if [[ ! -s "$r1" || ! -s "$r2" ]]; then
    fail BLOCKED_PAIRING "missing or empty paired FASTQ"
fi

r1_lines="$(wc -l < "$r1")"
r2_lines="$(wc -l < "$r2")"
if (( r1_lines % 4 != 0 || r2_lines % 4 != 0 )); then
    fail BLOCKED_PAIRING "FASTQ line counts are not multiples of 4"
fi
r1_reads="$((r1_lines / 4))"
r2_reads="$((r2_lines / 4))"
if [[ "$r1_reads" != "$r2_reads" || "$r1_reads" == "0" ]]; then
    fail BLOCKED_PAIRING "R1=$r1_reads R2=$r2_reads"
fi

r1_content="$(sha256sum "$r1" | awk '{print $1}')"
r2_content="$(sha256sum "$r2" | awk '{print $1}')"
pigz -n -p "$threads" -f "$r1" "$r2"
r1_gz="$DEST/${RUN}_1.fastq.gz"
r2_gz="$DEST/${RUN}_2.fastq.gz"
if [[ ! -s "$r1_gz" || ! -s "$r2_gz" ]]; then
    fail BLOCKED_CHECKSUM "compression did not produce both gzip files"
fi
r1_gz_sha="$(sha256sum "$r1_gz" | awk '{print $1}')"
r2_gz_sha="$(sha256sum "$r2_gz" | awk '{print $1}')"
r1_bytes="$(stat -c %s "$r1_gz")"
r2_bytes="$(stat -c %s "$r2_gz")"

if ! python3 - "$RUN" "$REPO_DIR" "$r1_reads" "$r1_content" "$r2_content" "$r1_gz_sha" "$r2_gz_sha" "$r1_bytes" "$r2_bytes" << 'PY' > "$PROV"
import csv, sys
from pathlib import Path
sys.path.insert(0, sys.argv[2] + "/benchmark/real_world/scripts")
run, repo, reads, c1, c2, z1, z2, b1, b2 = sys.argv[1:]
runinfo = list(csv.DictReader(open(Path(repo)/"benchmark/real_world/manifests/prjna1088182_sra_runinfo.csv")))
meta = next(row for row in runinfo if row["Run"] == run)
spots = int(meta["spots"])
observed = int(float(reads))
if observed == spots:
    consistency = "MATCH"
    status = "READY_FOR_EXECUTION"
    notes = "paired FASTQ count equals SRA spots"
elif observed > 0 and observed < spots:
    consistency = "UNRESOLVED_DIFFERENCE"
    status = "BLOCKED_PAIRING"
    notes = f"paired reads {observed} are below SRA spots {spots}"
else:
    consistency = "UNRESOLVED_DIFFERENCE"
    status = "BLOCKED_PAIRING"
    notes = f"paired reads {observed} disagree with SRA spots {spots}"
rel = f"data/prjna1088182/fastq/{run}"
fields = [
    "run","sample_name","biosample","experiment","submission","strain","library_layout",
    "execution_reference","reference_length","reference_sequence_sha256",
    "r1_path","r2_path","r1_reads","r2_reads","r1_uncompressed_sha256","r2_uncompressed_sha256",
    "r1_compressed_sha256","r2_compressed_sha256","r1_bytes","r2_bytes",
    "sra_toolkit_version","compression_tool","compression_version","slurm_job_id",
    "staging_status","spot_consistency","notes",
]
row = {
    "run": run,
    "sample_name": meta["SampleName"],
    "biosample": meta["BioSample"],
    "experiment": meta["Experiment"],
    "submission": meta["Submission"],
    "strain": meta["ScientificName"],
    "library_layout": meta["LibraryLayout"],
    "execution_reference": "NC_000913.3",
    "reference_length": "4641652",
    "reference_sequence_sha256": "e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0",
    "r1_path": f"{rel}/{run}_1.fastq.gz",
    "r2_path": f"{rel}/{run}_2.fastq.gz",
    "r1_reads": str(observed),
    "r2_reads": str(observed),
    "r1_uncompressed_sha256": c1,
    "r2_uncompressed_sha256": c2,
    "r1_compressed_sha256": z1,
    "r2_compressed_sha256": z2,
    "r1_bytes": b1,
    "r2_bytes": b2,
    "sra_toolkit_version": "3.4.1",
    "compression_tool": "pigz",
    "compression_version": "2.8 -n",
    "slurm_job_id": __import__("os").environ.get("SLURM_JOB_ID", "manual"),
    "staging_status": status,
    "spot_consistency": consistency,
    "notes": notes,
}
writer = csv.DictWriter(sys.stdout, fieldnames=fields, delimiter="\t", lineterminator="\n")
writer.writeheader()
writer.writerow(row)
if status != "READY_FOR_EXECUTION":
    raise SystemExit(3)
PY
then
    echo "end $(date -Iseconds) status_not_ready" | tee -a "$LOG"
    exit 3
fi
echo "end $(date -Iseconds)" | tee -a "$LOG"
echo "READY_FOR_EXECUTION $RUN" | tee -a "$LOG"
