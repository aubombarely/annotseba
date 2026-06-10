#!/usr/bin/env bash
# Quick launcher — run from inside the genome_pipeline/ directory.
# Usage: bash run_pipeline.sh [extra snakemake args]
#   e.g. bash run_pipeline.sh --dryrun
#   e.g. bash run_pipeline.sh --cores 16

set -euo pipefail

CORES="${SNAKEMAKE_CORES:-8}"

snakemake \
    --cores "$CORES" \
    --rerun-incomplete \
    --keep-going \
    "$@"
