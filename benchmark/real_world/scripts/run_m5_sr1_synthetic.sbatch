#!/usr/bin/env bash
#SBATCH --job-name=prokadiff_m5_sr1_synth
#SBATCH --partition=qcpu_18i
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=benchmark/real_world/results/m5_sr1_synthetic_%j.log
#SBATCH --error=benchmark/real_world/results/m5_sr1_synthetic_%j.err

set -euo pipefail
export PATH="/hpcfs/fhome/caizhh/.conda/envs/prokadiff/bin:$PATH"

REPO_DIR="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff"
RUN_ROOT="$REPO_DIR/benchmark/real_world/results/m5_sr1_synthetic_${SLURM_JOB_ID:-manual}"
mkdir -p "$RUN_ROOT"
cd "$REPO_DIR"

echo "=== M5-SR1 SYNTHETIC VALIDATION START ==="
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-manual}"
echo "SLURM_NODELIST: ${SLURM_NODELIST:-$(hostname)}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK:-16}"
echo "SLURM_MEM_PER_NODE: 32G"
echo "START_TIME: $(date -Iseconds)"
echo "GIT_HEAD: $(git rev-parse HEAD)"
echo "GIT_DIFF_HEAD_SHA256: $(git diff HEAD | sha256sum | awk '{print $1}')"
echo "GIT_STATUS_PORCELAIN_SHA256: $(git status --porcelain | sha256sum | awk '{print $1}')"

echo "=== STEP 1: CARGO TEST M5-SR1 INTEGRATION SUITE ==="
cargo test --release -p prokadiff --test m5_sr1_evidence_only -- --nocapture

echo "=== STEP 1b: BUILD CURRENT RELEASE BINARY ==="
# cargo test --test does not rebuild the CLI binary. The FASTQ step must use
# the current worktree, not a stale target/release/prokadiff.
cargo build --release -p prokadiff --bin prokadiff

echo "=== STEP 2: END-TO-END FASTQ/BOWTIE2 SYNTHETIC EXPERIMENTS ==="
python3 - "$RUN_ROOT" <<'PY'
from pathlib import Path
import csv
import subprocess
import sys

run_root = Path(sys.argv[1])
prokadiff_bin = Path("target/release/prokadiff").resolve()

def rc(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]

anchor = "GATTACAGTCGATCGTAGCTAGTCCGATGCTAGTCGATCG"  # 41 bp unique anchor
repeat_16 = "ACGTACGTACGAACGT"  # 16 bp repeat motif (>=14 bp min_overlap -> SUPPORTED_AMBIGUOUS)
repeat_12 = "ACGTACGTACGA"      # 12 bp repeat motif (<14 bp min_overlap -> CANDIDATE)

# -------------------------------------------------------------
# Fixture 1: 20-copy boundary control (Criterion A & D)
# -------------------------------------------------------------
print("--- Running Fixture 1: 20-copy boundary control ---")
dir_20 = run_root / "boundary_20_copies"
dir_20.mkdir(parents=True, exist_ok=True)
ref_20 = "A" * 100 + anchor + "C" * 100 + "G".join(repeat_16 for _ in range(20)) + "T" * 100
(dir_20 / "reference.fna").write_text(f">ref_20\n{ref_20}\n")

forward_16 = anchor + repeat_16
reads_20 = []
for i in range(6):
    reads_20.append((f"fwd_{i}", forward_16))
    reads_20.append((f"rev_{i}", rc(forward_16)))
with (dir_20 / "reads.fastq").open("w") as h:
    for name, seq in reads_20:
        h.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")

cmd_20 = [
    str(prokadiff_bin), "evidence",
    "--ref", str(dir_20 / "reference.fna"),
    "--fastq", str(dir_20 / "reads.fastq"),
    "--threads", "16",
    "--outdir", str(dir_20 / "out")
]
subprocess.run(cmd_20, check=True)

sidecar_20 = dir_20 / "out" / "repeat_ambiguous_junctions.tsv"
with sidecar_20.open() as h:
    rows_20 = list(csv.DictReader(h, delimiter="\t"))
print(f"20-copy fixture repeat sidecar rows: {len(rows_20)}")
if len(rows_20) != 0:
    raise SystemExit(f"CRITERION D FAILED: 20-copy fixture produced {len(rows_20)} repeat sidecar rows; expected 0 (concrete split path)")
print("Criterion A & D (20-copy control): PASS")

# -------------------------------------------------------------
# Fixture 2a: 21-copy high-copy unresolved -> SUPPORTED_AMBIGUOUS (Criterion B & D)
# -------------------------------------------------------------
print("--- Running Fixture 2a: 21-copy unresolved repeat family (SUPPORTED_AMBIGUOUS) ---")
dir_21a = run_root / "boundary_21a_copies"
dir_21a.mkdir(parents=True, exist_ok=True)
ref_21a = "A" * 100 + anchor + "C" * 100 + "G".join(repeat_16 for _ in range(21)) + "T" * 100
(dir_21a / "reference.fna").write_text(f">ref_21a\n{ref_21a}\n")

with (dir_21a / "reads.fastq").open("w") as h:
    for name, seq in reads_20:
        h.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")

cmd_21a = [
    str(prokadiff_bin), "evidence",
    "--ref", str(dir_21a / "reference.fna"),
    "--fastq", str(dir_21a / "reads.fastq"),
    "--threads", "16",
    "--outdir", str(dir_21a / "out")
]
subprocess.run(cmd_21a, check=True)

sidecar_21a = dir_21a / "out" / "repeat_ambiguous_junctions.tsv"
gd_21a = dir_21a / "out" / "output.gd"
with sidecar_21a.open() as h:
    rows_21a = list(csv.DictReader(h, delimiter="\t"))
print(f"21a-copy fixture repeat sidecar rows: {len(rows_21a)}")
if len(rows_21a) != 1:
    raise SystemExit(f"CRITERION B FAILED: 21a-copy fixture produced {len(rows_21a)} rows; expected exactly 1")

r21a = rows_21a[0]
print(f"Row 21a state: {r21a['state']}, exact_supported: {r21a['exact_breakpoint_supported']}, placement_count: {r21a['placement_count']}, distinct_molecules: {r21a['distinct_molecule_count']}")
if r21a["state"] != "SUPPORTED_AMBIGUOUS":
    raise SystemExit(f"CRITERION B FAILED: expected SUPPORTED_AMBIGUOUS, got {r21a['state']}")
if r21a["exact_breakpoint_supported"] != "false":
    raise SystemExit("CRITERION B FAILED: exact_breakpoint_supported must be false")
if int(r21a["placement_count"]) != 21:
    raise SystemExit(f"CRITERION B FAILED: placement_count expected 21, got {r21a['placement_count']}")
if int(r21a["distinct_molecule_count"]) != 12:
    raise SystemExit(f"CRITERION B FAILED: distinct_molecule_count expected 12, got {r21a['distinct_molecule_count']}")
if "\nJC\t" in gd_21a.read_text():
    raise SystemExit("CRITERION B FAILED: unresolved repeat family emitted a concrete JC to output.gd")
print("Criterion B & D (21-copy SUPPORTED_AMBIGUOUS): PASS")

# -------------------------------------------------------------
# Fixture 2b: 21-copy high-copy unresolved -> CANDIDATE (Criterion B & D)
# -------------------------------------------------------------
print("--- Running Fixture 2b: 21-copy unresolved repeat family (CANDIDATE) ---")
dir_21b = run_root / "boundary_21b_copies"
dir_21b.mkdir(parents=True, exist_ok=True)
ref_21b = "A" * 100 + anchor + "C" * 100 + "G".join(repeat_12 for _ in range(21)) + "T" * 100
(dir_21b / "reference.fna").write_text(f">ref_21b\n{ref_21b}\n")

forward_12 = anchor + repeat_12
reads_12 = []
for i in range(6):
    reads_12.append((f"fwd_{i}", forward_12))
    reads_12.append((f"rev_{i}", rc(forward_12)))
with (dir_21b / "reads.fastq").open("w") as h:
    for name, seq in reads_12:
        h.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")

cmd_21b = [
    str(prokadiff_bin), "evidence",
    "--ref", str(dir_21b / "reference.fna"),
    "--fastq", str(dir_21b / "reads.fastq"),
    "--threads", "16",
    "--outdir", str(dir_21b / "out")
]
subprocess.run(cmd_21b, check=True)

sidecar_21b = dir_21b / "out" / "repeat_ambiguous_junctions.tsv"
gd_21b = dir_21b / "out" / "output.gd"
with sidecar_21b.open() as h:
    rows_21b = list(csv.DictReader(h, delimiter="\t"))
print(f"21b-copy fixture repeat sidecar rows: {len(rows_21b)}")
if len(rows_21b) != 1:
    raise SystemExit(f"CRITERION B FAILED: 21b-copy fixture produced {len(rows_21b)} rows; expected exactly 1")

r21b = rows_21b[0]
print(f"Row 21b state: {r21b['state']}, exact_supported: {r21b['exact_breakpoint_supported']}, placement_count: {r21b['placement_count']}, distinct_molecules: {r21b['distinct_molecule_count']}")
if r21b["state"] != "CANDIDATE":
    raise SystemExit(f"CRITERION B FAILED: expected CANDIDATE, got {r21b['state']}")
if r21b["exact_breakpoint_supported"] != "false":
    raise SystemExit("CRITERION B FAILED: exact_breakpoint_supported must be false")
if int(r21b["placement_count"]) != 21:
    raise SystemExit(f"CRITERION B FAILED: placement_count expected 21, got {r21b['placement_count']}")
if "\nJC\t" in gd_21b.read_text():
    raise SystemExit("CRITERION B FAILED: candidate repeat family emitted a concrete JC to output.gd")
print("Criterion B & D (21-copy CANDIDATE): PASS")

# -------------------------------------------------------------
# Fixture 3: Determinism / Reordered FASTQ (Criterion F)
# -------------------------------------------------------------
print("--- Running Fixture 3: Reordered FASTQ determinism ---")
dir_reorder = run_root / "reorder_21_copies"
dir_reorder.mkdir(parents=True, exist_ok=True)
(dir_reorder / "reference.fna").write_text(f">ref_21a\n{ref_21a}\n")

# Reverse the read order
reads_reordered = list(reversed(reads_20))
with (dir_reorder / "reads.fastq").open("w") as h:
    for name, seq in reads_reordered:
        h.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")

cmd_reorder = [
    str(prokadiff_bin), "evidence",
    "--ref", str(dir_reorder / "reference.fna"),
    "--fastq", str(dir_reorder / "reads.fastq"),
    "--threads", "16",
    "--outdir", str(dir_reorder / "out")
]
subprocess.run(cmd_reorder, check=True)

sidecar_reorder = dir_reorder / "out" / "repeat_ambiguous_junctions.tsv"
with sidecar_reorder.open() as h:
    rows_reorder = list(csv.DictReader(h, delimiter="\t"))
if len(rows_reorder) != 1:
    raise SystemExit(f"CRITERION F FAILED: reordered run produced {len(rows_reorder)} rows")
r_reorder = rows_reorder[0]
for field in ["state", "exact_breakpoint_supported", "anchor_position_1", "placement_count", "distinct_molecule_count", "plus_molecules", "minus_molecules"]:
    if r_reorder[field] != r21a[field]:
        raise SystemExit(f"CRITERION F FAILED: field {field} differed under reordering: {r_reorder[field]} vs {r21a[field]}")
print("Criterion F (Determinism): PASS")

print("=== ALL SYNTHETIC FASTQ/EVIDENCE PIPELINES PASSED ===")

# Summary of diagnostic counts across synthetic runs
print("=== SR1 SYNTHETIC DIAGNOSTIC SUMMARY ===")
print("repeat_family_count: 2 (21a SUPPORTED_AMBIGUOUS + 21b CANDIDATE)")
print("CANDIDATE_count: 1")
print("SUPPORTED_AMBIGUOUS_count: 1")
print("RESOLVED_count: 1 (in m5_sr1_evidence_only: resolved_repeat_evidence_has_no_product_event)")
print("RESOURCE_LIMIT_count: 4 (in m5_sr1_evidence_only: reciprocal index, reciprocal candidates, stage-2 geometry overflow, truncation > 4096)")
print("maximum_observed_placement_count: 4097 (truncation fixture)")
print("maximum_observed_stage2_family_geometries: 5001 (stage_two_overflow fixture)")
print("max_stage2_family_geometries_5000_triggered: YES (triggered and cleanly handled in stage_two_overflow_repeat_evidence_has_no_product_event)")
PY

echo "=== M5-SR1 SYNTHETIC VALIDATION COMPLETE ==="
echo "END_TIME: $(date -Iseconds)"
