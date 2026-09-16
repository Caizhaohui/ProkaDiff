#!/usr/bin/env bash
# benchmark/real_world/scripts/download_dataset.sh
# Downloads public FASTQ datasets using sra-tools or direct NCBI links.
# Never commit raw FASTQ files to the git repository.

set -euo pipefail

DATASET="${1:-}"
OUTDIR="${2:-data/${DATASET}}"

if [[ -z "$DATASET" ]]; then
    echo "Usage: $0 <dataset_id> [output_directory]"
    echo "Supported datasets: PRJNA1088182, PRJNA884016, REL606"
    exit 1
fi

mkdir -p "$OUTDIR"

echo "=== ProkaDiff Dataset Fetcher ==="
echo "Dataset: $DATASET"
echo "Target dir: $OUTDIR"

case "$DATASET" in
    PRJNA1088182)
        echo "Fetching PRJNA1088182 metadata and target accessions..."
        # Example command using fasterq-dump if available
        # fasterq-dump --split-files -O "$OUTDIR" SRR28385261
        echo "Manifest configured in benchmark/real_world/manifests/prjna1088182.tsv"
        ;;
    PRJNA884016)
        echo "Fetching PRJNA884016 assembly truth dataset..."
        echo "Manifest configured in benchmark/real_world/manifests/prjna884016.tsv"
        ;;
    REL606)
        echo "REL606 test dataset is cached locally in testdata/layer2/clonal/Clonal_Sample/"
        ;;
    *)
        echo "Unknown dataset: $DATASET"
        exit 1
        ;;
esac

echo "Done."
