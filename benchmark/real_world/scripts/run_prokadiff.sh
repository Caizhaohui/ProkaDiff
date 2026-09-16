#!/usr/bin/env bash
# benchmark/real_world/scripts/run_prokadiff.sh
# Runs ProkaDiff differential audit with runtime and peak RSS telemetry.
# Strict Rule 11: Execute on Slurm compute node (qcpu_18i).

set -euo pipefail

PROKADIFF_BIN="${PROKADIFF_BIN:-./target/release/prokadiff}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

STARTER_R1=""
STARTER_R2=""
EDITED_R1=""
EDITED_R2=""
REF=""
OUTDIR=""
INTENDED=""
EDITOR="cas9"
SPACER=""
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --starter-r1) STARTER_R1="$2"; shift 2 ;;
        --starter-r2) STARTER_R2="$2"; shift 2 ;;
        --edited-r1) EDITED_R1="$2"; shift 2 ;;
        --edited-r2) EDITED_R2="$2"; shift 2 ;;
        --ref) REF="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --intended) INTENDED="$2"; shift 2 ;;
        --editor) EDITOR="$2"; shift 2 ;;
        --spacer) SPACER="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        *) EXTRA_ARGS+=("$1"); shift ;;
    esac
done

if [[ -z "$STARTER_R1" || -z "$EDITED_R1" || -z "$REF" || -z "$OUTDIR" ]]; then
    echo "Usage: $0 --starter-r1 <R1> [--starter-r2 <R2>] --edited-r1 <R1> [--edited-r2 <R2>] --ref <REF> --outdir <DIR> [--intended <TSV>] [--editor <NAME>] [--spacer <SEQ>] [--threads <N>]"
    exit 1
fi

mkdir -p "$OUTDIR"

CMD=(
    "$PROKADIFF_BIN"
    --starter "$STARTER_R1"
)
if [[ -n "$STARTER_R2" ]]; then
    CMD+=( "$STARTER_R2" )
fi

CMD+=( --edited "$EDITED_R1" )
if [[ -n "$EDITED_R2" ]]; then
    CMD+=( "$EDITED_R2" )
fi

CMD+=(
    --ref "$REF"
    --outdir "$OUTDIR"
    --threads "$THREADS"
    --editor "$EDITOR"
)

if [[ -n "$INTENDED" ]]; then
    CMD+=( --intended "$INTENDED" )
fi

if [[ -n "$SPACER" ]]; then
    CMD+=( --spacer "$SPACER" )
fi

if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
    CMD+=( "${EXTRA_ARGS[@]}" )
fi

echo "=== Running ProkaDiff ==="
echo "Command: ${CMD[*]}"
echo "Start time: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"

/usr/bin/time -v -o "${OUTDIR}/time.log" "${CMD[@]}" > "${OUTDIR}/prokadiff.stdout" 2> "${OUTDIR}/prokadiff.stderr"

echo "End time: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
echo "ProkaDiff finished successfully. Output in ${OUTDIR}"
