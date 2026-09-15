#!/usr/bin/env bash
# run_flashfry.sh — Run FlashFry 1.15 against a reference FASTA and spacer guide.
# STATUS: BLOCKED_EXTERNAL_DEPENDENCY — FlashFry 1.15 is not installed.
#
# Usage:
#   bash run_flashfry.sh <reference.fa> <guide_sequence> <outdir>
#
# FlashFry produces a TSV with columns including:
#   contig  start  stop  target  orientation  Doench2016CFDScore  Hsu2013
#
# Expected FlashFry version: 1.15
# Download: https://github.com/mckennalab/FlashFry/releases/tag/1.15
set -euo pipefail

FLASHFRY_JAR="${FLASHFRY_JAR:-}"
REFERENCE="${1:-}"
GUIDE="${2:-}"
OUTDIR="${3:-$(pwd)/flashfry_output}"

# ── dependency check ────────────────────────────────────────────────────────
if [[ -z "$FLASHFRY_JAR" ]] || [[ ! -f "$FLASHFRY_JAR" ]]; then
    echo "BLOCKED_EXTERNAL_DEPENDENCY: FlashFry JAR not found." >&2
    echo "Set FLASHFRY_JAR=/path/to/flashfry-1.15.jar before running." >&2
    exit 2
fi

if ! command -v java &>/dev/null; then
    echo "BLOCKED_EXTERNAL_DEPENDENCY: java not found in PATH." >&2
    exit 2
fi

if [[ -z "$REFERENCE" ]] || [[ ! -f "$REFERENCE" ]]; then
    echo "Error: reference FASTA not found: $REFERENCE" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

DB_PREFIX="$OUTDIR/flashfry_db"
OUTPUT_TSV="$OUTDIR/raw_flashfry.tsv"

# Step 1: Build FlashFry database
java -jar "$FLASHFRY_JAR" index \
    --tmpLocation "$OUTDIR/tmp" \
    --database "$DB_PREFIX" \
    --reference "$REFERENCE" \
    --enzyme spcas9ngg

# Step 2: Discover off-target sites
java -jar "$FLASHFRY_JAR" discover \
    --database "$DB_PREFIX" \
    --fasta "$REFERENCE" \
    --output "$OUTPUT_TSV" \
    --guide "$GUIDE" \
    --maxMismatch 4 \
    --flankingSequence 10

echo "[run_flashfry] Raw output written to: $OUTPUT_TSV"

# Step 3: Score sites for CFD and Hsu2013
SCORED_TSV="$OUTDIR/raw_flashfry_scored.tsv"
java -jar "$FLASHFRY_JAR" score \
    --input "$OUTPUT_TSV" \
    --output "$SCORED_TSV" \
    --scoringMetrics doench2016cfd,hsu2013

echo "[run_flashfry] Scored output written to: $SCORED_TSV"

# Copy to golden/
mkdir -p "$(dirname "$0")/../golden"
cp "$SCORED_TSV" "$(dirname "$0")/../golden/raw_flashfry.tsv"
echo "[run_flashfry] Archived to golden/raw_flashfry.tsv"
