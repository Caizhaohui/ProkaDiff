#!/usr/bin/env bash
#SBATCH --job-name=m6p3a_locus
#SBATCH --partition=qcpu_18i
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-2%2
#SBATCH --output=benchmark/real_world/results/prjna1088182/m6_p3a/slurm_%A_%a.log
#SBATCH --error=benchmark/real_world/results/prjna1088182/m6_p3a/slurm_%A_%a.err

# Stage-1 Bowtie2 only. Not prokadiff evidence and not a junction pass.
# BAM sort uses samtools, not the production in-memory noodles sort.
set -euo pipefail
export PATH="/hpcfs/fhome/caizhh/.conda/envs/prokadiff/bin:${PATH}"
REPO="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff"
SAMTOOLS="$REPO/data/tools/samtools-env/bin/samtools"
cd "$REPO"
RUNS=(SRR28349204 SRR28349176 SRR28349230)
RUN="${RUNS[$SLURM_ARRAY_TASK_ID]}"
OUT="$REPO/benchmark/real_world/results/prjna1088182/m6_p3a/$RUN"
WORK="$OUT/work"
mkdir -p "$OUT" "$WORK"
echo "job ${SLURM_JOB_ID:-manual} task ${SLURM_ARRAY_TASK_ID} run $RUN host $(hostname) start $(date -Iseconds)"
echo "bowtie2 $(bowtie2 --version | head -1)"
echo "samtools $($SAMTOOLS --version | head -1)"
avail="$(df -B1 --output=avail "$OUT" | tail -1 | tr -d ' ')"
need=$((40 * 1024 * 1024 * 1024))
if [[ "$avail" -lt "$need" ]]; then
    echo "free bytes $avail below 40GiB" >&2
    exit 2
fi

python3 - "$RUN" "$REPO" "$WORK" << 'PY'
import csv, hashlib, sys
from pathlib import Path
run, repo_s, work_s = sys.argv[1:]
repo = Path(repo_s)
sys.path.insert(0, str(repo / "benchmark/real_world/scripts"))
from plan_m6_p3a import build_targets
from m6_p2a_stage_gate import chromosome_bases, verify_local_reference
targets = [row for row in build_targets() if row["run"] == run and row["diagnostic_priority"] == "1"]
if not targets:
    raise SystemExit(f"{run} has no priority-1 target")
refs = {row["reference"] for row in targets}
if len(refs) != 1:
    raise SystemExit(f"{run} mixed references {refs}")
accession = refs.pop()
digest = verify_local_reference(accession)
prov_path = repo / "data/prjna1088182/fastq" / run / "run_provenance.tsv"
prov = list(csv.DictReader(prov_path.open(encoding="utf-8"), delimiter="\t"))
if len(prov) != 1 or prov[0]["run"] != run:
    raise SystemExit("provenance identity failed")
row = prov[0]
if row["sample_name"] != targets[0]["sample"] or row["execution_reference"] != accession:
    raise SystemExit("sample or reference mismatch")
if row["staging_status"] != "READY_FOR_EXECUTION" or row["r1_reads"] != row["r2_reads"]:
    raise SystemExit("staging is not a matched pair")
for label in ("r1", "r2"):
    path = repo / row[f"{label}_path"]
    got = hashlib.sha256(path.read_bytes()).hexdigest()
    if got != row[f"{label}_compressed_sha256"]:
        raise SystemExit(f"{label} checksum mismatch")
print(f"gate_ok {run} {targets[0]['sample']} {accession} {digest} reads={row['r1_reads']}")
fasta = Path(work_s) / "reference.fa"
if accession == "NC_000913.3":
    bases = chromosome_bases(repo / "data/references/MG1655.gbk")
    fasta.write_text(f">NC_000913.3\n" + "\n".join(bases[i:i+80] for i in range(0, len(bases), 80)) + "\n", encoding="utf-8")
else:
    fasta.write_text((repo / "data/references/m6_identity/CP009273.1.fasta").read_text(encoding="utf-8"), encoding="utf-8")
Path(work_s, "targets.tsv").write_text(
    "coordinate\tref\talt\tsample\treference\n" + "".join(
        f"{item['coordinate']}\t{item['ref']}\t{item['alt']}\t{item['sample']}\t{item['reference']}\n" for item in targets
    ),
    encoding="utf-8",
)
PY

FASTA="$WORK/reference.fa"
bowtie2-build --threads "${SLURM_CPUS_PER_TASK:-8}" "$FASTA" "$WORK/bt2ref"
bowtie2 --local --no-unal -p "${SLURM_CPUS_PER_TASK:-8}" -x "$WORK/bt2ref" \
  -1 "$REPO/data/prjna1088182/fastq/$RUN/${RUN}_1.fastq.gz" \
  -2 "$REPO/data/prjna1088182/fastq/$RUN/${RUN}_2.fastq.gz" \
  --un-conc "$WORK/unconc.fastq" --un "$WORK/unpaired.fastq" \
  -S "$WORK/aligned.sam"
"$SAMTOOLS" sort -@ "${SLURM_CPUS_PER_TASK:-8}" -o "$OUT/aligned.bam" "$WORK/aligned.sam"
"$SAMTOOLS" index "$OUT/aligned.bam"
{
  echo "run $RUN"
  echo "bowtie2 $(bowtie2 --version | head -1)"
  echo "samtools $($SAMTOOLS --version | head -1)"
  echo "command bowtie2 --local --no-unal -p ${SLURM_CPUS_PER_TASK:-8} -x $WORK/bt2ref -1 $REPO/data/prjna1088182/fastq/$RUN/${RUN}_1.fastq.gz -2 $REPO/data/prjna1088182/fastq/$RUN/${RUN}_2.fastq.gz --un-conc $WORK/unconc.fastq --un $WORK/unpaired.fastq -S $WORK/aligned.sam"
  echo "bam_sort samtools sort then samtools index; production uses noodles in-memory sort and does not write a bai"
  echo "slurm_job_id ${SLURM_JOB_ID:-manual}"
} > "$OUT/command.txt"

python3 - "$RUN" "$OUT" "$REPO" << 'PY'
import csv, subprocess, sys
from pathlib import Path
run, out_s, repo_s = sys.argv[1:]
out = Path(out_s)
repo = Path(repo_s)
samtools = repo / "data/tools/samtools-env/bin/samtools"
targets = list(csv.DictReader((out / "work/targets.tsv").open(), delimiter="\t"))
rows = []
for target in targets:
    view = out / "work" / f"locus_{target['coordinate']}.sam"
    name = target["reference"]
    region = f"{name}:{target['coordinate']}-{target['coordinate']}"
    with view.open("w", encoding="utf-8") as handle:
        subprocess.check_call([str(samtools), "view", "-h", str(out / "aligned.bam"), region], stdout=handle)
    dest = out / "work" / f"locus_{target['coordinate']}.tsv"
    subprocess.check_call([
        "python3", str(repo / "benchmark/real_world/scripts/measure_m6_p3a_locus.py"),
        "--sam", str(view), "--position", target["coordinate"], "--ref", target["ref"], "--alt", target["alt"],
        "--run", run, "--sample", target["sample"], "--reference", target["reference"], "--out", str(dest),
    ])
    rows.extend(csv.DictReader(dest.open(), delimiter="\t"))
fields = list(rows[0])
with (out / "locus_pileup.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
print(f"loci {len(rows)} written")
PY
echo "run $RUN end $(date -Iseconds)"
