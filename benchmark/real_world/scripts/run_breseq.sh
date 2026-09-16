#!/usr/bin/env bash
# benchmark/real_world/scripts/run_breseq.sh
# Runs breseq oracle with runtime and peak RSS telemetry.
# breseq is ORACLE ONLY (Rule 4 & 5). Never use as runtime engine.

set -euo pipefail

BRESEQ_BIN="${BRESEQ_BIN:-/hpcfs/fhome/caizhh/.conda/envs/prokadiff/bin/breseq}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

R1=""
R2=""
REF=""
OUTDIR=""
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --r1) R1="$2"; shift 2 ;;
        --r2) R2="$2"; shift 2 ;;
        --ref) REF="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        *) EXTRA_ARGS+=("$1"); shift ;;
    esac
done

if [[ -z "$R1" || -z "$REF" || -z "$OUTDIR" ]]; then
    echo "Usage: $0 --r1 <R1> [--r2 <R2>] --ref <REF> --outdir <DIR> [--threads <N>]"
    exit 1
fi

mkdir -p "$OUTDIR"

CMD=(
    "$BRESEQ_BIN"
    -j "$THREADS"
    -r "$REF"
    -o "$OUTDIR"
    "$R1"
)

if [[ -n "$R2" ]]; then
    CMD+=( "$R2" )
fi

if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
    CMD+=( "${EXTRA_ARGS[@]}" )
fi

echo "=== Running breseq Oracle ==="
echo "Command: ${CMD[*]}"
echo "Start time: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"

/usr/bin/time -v -o "${OUTDIR}/time.log" "${CMD[@]}" > "${OUTDIR}/breseq.stdout" 2> "${OUTDIR}/breseq.stderr"

echo "End time: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
echo "breseq oracle finished. Output in ${OUTDIR}"
