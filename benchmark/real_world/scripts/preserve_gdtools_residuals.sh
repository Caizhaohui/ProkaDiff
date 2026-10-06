#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "usage: $0 LEFT.gd RIGHT.gd OUTPUT_DIR LABEL" >&2
    exit 2
fi
left="$1"
right="$2"
out_dir="$3"
label="$4"
mkdir -p "$out_dir"
left="$(realpath "$left")"
right="$(realpath "$right")"
(
    cd "$out_dir"
    gdtools SUBTRACT "$left" "$right"
    mv -f output.gd "${label}_A_minus_B.gd"
    gdtools SUBTRACT "$right" "$left"
    mv -f output.gd "${label}_B_minus_A.gd"
)
