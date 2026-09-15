#!/usr/bin/env bash
# run_crispritz.sh — Run CRISPRitz 2.8.1 against a reference and guide with bulge support.
# STATUS: BLOCKED_EXTERNAL_DEPENDENCY — CRISPRitz 2.8.1 is not installed.
#
# Usage:
#   bash run_crispritz.sh <reference.fa> <guide_file.txt> <outdir> [max_dna_bulge] [max_rna_bulge]
#
# CRISPRitz is available from:
#   https://github.com/pinellolab/CRISPRitz
set -euo pipefail

CRISPRITZ="${CRISPRITZ:-crispritz}"
REFERENCE="${1:-}"
GUIDE_FILE="${2:-}"
OUTDIR="${3:-$(pwd)/crispritz_output}"
MAX_DNA_BULGE="${4:-2}"
MAX_RNA_BULGE="${5:-2}"

# ── dependency check ────────────────────────────────────────────────────────
if ! command -v "$CRISPRITZ" &>/dev/null; then
    echo "BLOCKED_EXTERNAL_DEPENDENCY: crispritz binary not found in PATH." >&2
    echo "Install CRISPRitz 2.8.1: https://github.com/pinellolab/CRISPRitz" >&2
    exit 2
fi

if [[ -z "$REFERENCE" ]] || [[ ! -f "$REFERENCE" ]]; then
    echo "Error: reference FASTA not found: $REFERENCE" >&2
    exit 1
fi

if [[ -z "$GUIDE_FILE" ]] || [[ ! -f "$GUIDE_FILE" ]]; then
    echo "Error: guide file not found: $GUIDE_FILE" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

# CRISPRitz search with bulge support
"$CRISPRITZ" search \
    "$REFERENCE" \
    "$GUIDE_FILE" \
    "$OUTDIR/crispritz_result" \
    -cas_alg SpCas9 \
    -mm 4 \
    -bDNA "$MAX_DNA_BULGE" \
    -bRNA "$MAX_RNA_BULGE" \
    -r

RAW_OUT="$OUTDIR/crispritz_result.targets.txt"
mkdir -p "$(dirname "$0")/../golden"

if [[ -f "$RAW_OUT" ]]; then
    cp "$RAW_OUT" "$(dirname "$0")/../golden/raw_crispritz.txt"
    echo "[run_crispritz] Raw output archived to golden/raw_crispritz.txt"
else
    echo "Warning: expected output $RAW_OUT not found" >&2
fi
