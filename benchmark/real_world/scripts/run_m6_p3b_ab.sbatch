#!/usr/bin/env bash
#SBATCH --job-name=m6p3b_ab
#SBATCH --partition=qcpu_18i
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-2%3
#SBATCH --output=benchmark/real_world/results/prjna1088182/m6_p3b/slurm_%A_%a.log
#SBATCH --error=benchmark/real_world/results/prjna1088182/m6_p3b/slurm_%A_%a.err

set -euo pipefail
export PATH="/hpcfs/fhome/caizhh/.conda/envs/prokadiff/bin:$PATH"
export RUST_BACKTRACE=1

REPO_DIR="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff"
cd "$REPO_DIR"

PROKADIFF="$REPO_DIR/target/release/prokadiff"
EXPECTED_BINARY_SHA256="02077c761b1c9bd83cf1b626a3e2e955d8d2a1c63abb3e35852421e0dbc20bfd"
BASE="$REPO_DIR/benchmark/real_world/results/prjna1088182/m6_p3b"
mkdir -p "$BASE"

RUNS=(SRR28349226 SRR28349231 SRR28349230)
RUN="${RUNS[$SLURM_ARRAY_TASK_ID]}"
OUT="$BASE/$RUN"
OUT_A="$OUT/cond_A"
OUT_B="$OUT/cond_B"
mkdir -p "$OUT" "$OUT_A" "$OUT_B"

echo "job ${SLURM_JOB_ID:-manual} task ${SLURM_ARRAY_TASK_ID} host $(hostname) run $RUN start $(date -Iseconds)"
echo "commit $(git rev-parse HEAD)"

# Binary verification
got_sha="$(sha256sum "$PROKADIFF" | awk '{print $1}')"
if [[ "$got_sha" != "$EXPECTED_BINARY_SHA256" ]]; then
    echo "refusing to run binary $got_sha; expected $EXPECTED_BINARY_SHA256" >&2
    exit 5
fi
echo "binary_sha256 $got_sha"
echo "bowtie2 $(bowtie2 --version | head -1)"
echo "cpus ${SLURM_CPUS_PER_TASK:-8} mem 32G"

# Free space verification
avail="$(df -B1 --output=avail "$OUT" | tail -1 | tr -d ' ')"
need=$((40 * 1024 * 1024 * 1024))
if [[ "$avail" -lt "$need" ]]; then
    echo "free bytes $avail below 40GiB" >&2
    exit 2
fi

# References
REF_A="$REPO_DIR/data/references/m6_identity/CP009273.1.fasta"
REF_B="$REPO_DIR/data/references/m6_identity/CP009273.1.gbk"
if [[ ! -s "$REF_A" ]]; then
    echo "missing reference A: $REF_A" >&2
    exit 3
fi
if [[ ! -s "$REF_B" ]]; then
    echo "missing reference B: $REF_B" >&2
    exit 3
fi

# FASTQs
R1="$REPO_DIR/data/prjna1088182/fastq/$RUN/${RUN}_1.fastq.gz"
R2="$REPO_DIR/data/prjna1088182/fastq/$RUN/${RUN}_2.fastq.gz"
if [[ ! -s "$R1" || ! -s "$R2" ]]; then
    echo "missing fastqs for $RUN: $R1 or $R2" >&2
    exit 4
fi

# === Condition A: FASTA Baseline (opts.repeats == 0) ===
if [[ -s "$OUT_A/output.gd" ]]; then
    echo "cond_A already complete for $RUN"
elif [[ -s "$REPO_DIR/benchmark/real_world/results/prjna1088182/m6_p2a/$RUN/output.gd" && "${P3B_REUSE_P2A:-1}" == "1" ]]; then
    echo "seeding cond_A from verified P2a baseline for $RUN"
    mkdir -p "$OUT_A/work"
    cp -r "$REPO_DIR/benchmark/real_world/results/prjna1088182/m6_p2a/$RUN/"* "$OUT_A/"
    cp "$REF_A" "$OUT_A/work/reference.fa"
    {
        echo "condition Condition_A_FASTA_baseline"
        echo "reference $REF_A"
        echo "run $RUN"
        echo "source seeded_from_p2a_verified_baseline"
        echo "binary_sha256 $got_sha"
    } > "$OUT_A/command.txt"
else
    echo "running cond_A for $RUN start $(date -Iseconds)"
    mkdir -p "$OUT_A/work"
    {
        echo "condition Condition_A_FASTA_baseline"
        echo "reference $REF_A"
        echo "run $RUN"
        echo "binary_sha256 $got_sha"
        echo "command $PROKADIFF evidence --ref $REF_A --fastq $R1 $R2 --threads ${SLURM_CPUS_PER_TASK:-8} --keep-bam --outdir $OUT_A"
    } > "$OUT_A/command.txt"
    /usr/bin/time -v -o "$OUT_A/telemetry.txt" \
        "$PROKADIFF" evidence \
        --ref "$REF_A" \
        --fastq "$R1" "$R2" \
        --threads "${SLURM_CPUS_PER_TASK:-8}" \
        --keep-bam \
        --outdir "$OUT_A"
    cp "$REF_A" "$OUT_A/work/reference.fa"
    echo "cond_A complete for $RUN at $(date -Iseconds)"
fi

# === Condition B: GenBank Annotation (opts.repeats == 401) ===
if [[ -s "$OUT_B/output.gd" ]]; then
    echo "cond_B already complete for $RUN"
else
    echo "running cond_B for $RUN start $(date -Iseconds)"
    mkdir -p "$OUT_B"
    {
        echo "condition Condition_B_GenBank_annotation"
        echo "reference $REF_B"
        echo "run $RUN"
        echo "binary_sha256 $got_sha"
        echo "command $PROKADIFF evidence --ref $REF_B --fastq $R1 $R2 --threads ${SLURM_CPUS_PER_TASK:-8} --keep-bam --outdir $OUT_B"
    } > "$OUT_B/command.txt"
    /usr/bin/time -v -o "$OUT_B/telemetry.txt" \
        "$PROKADIFF" evidence \
        --ref "$REF_B" \
        --fastq "$R1" "$R2" \
        --threads "${SLURM_CPUS_PER_TASK:-8}" \
        --keep-bam \
        --outdir "$OUT_B"
    echo "cond_B complete for $RUN at $(date -Iseconds)"
fi

echo "task ${SLURM_ARRAY_TASK_ID} run $RUN finished successfully at $(date -Iseconds)"
