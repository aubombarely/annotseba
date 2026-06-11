#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

VERSION=$(cat "$SCRIPT_DIR/VERSION" 2>/dev/null || echo "unknown")

usage() {
    cat <<EOF
annotseba v${VERSION}

A Snakemake pipeline to download eukaryotic genomes from NCBI and assess
their quality and completeness, including assembly statistics and annotation
quality evaluation.

Usage:
    bash run_annotseba.sh [options] [snakemake options]

Options:
    -h, --help            Show this help message and exit
    -v, --version         Show version and exit
    -n, --dryrun          Dry-run: show what would be executed without running
    -c, --cores INT       Number of CPU cores to use (default: \$SNAKEMAKE_CORES or 8)
    -f, --config_file FILE Path to a config YAML file (merged on top of default config/config.yaml)
    -a, --accessions FILE Path to accessions TSV file (overrides config.yaml)
    --keep_source         Keep raw NCBI files ({acc}.fna and {acc}.gff3) after renaming

Snakemake pass-through:
    Any additional arguments are passed directly to Snakemake.
    e.g. --forceall, --until <rule>, --rerun-triggers mtime

Examples:
    bash run_annotseba.sh --cores 16
    bash run_annotseba.sh --dryrun
    bash run_annotseba.sh --accessions my_species.tsv --cores 16
    bash run_annotseba.sh --cores 16 --forceall
    bash run_annotseba.sh --until run_busco

Input:
    Edit accessions.txt and config/config.yaml before running.
    See README.md for full documentation.
EOF
}

# Convert a path to absolute using the caller's working directory.
# Must be called before cd-ing to SCRIPT_DIR.
abspath() { [[ "$1" = /* ]] && echo "$1" || echo "$PWD/$1"; }

CORES="${SNAKEMAKE_CORES:-8}"
SNAKEMAKE_ARGS=()

# Parse arguments — resolve file paths to absolute here, before cd
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        -v|--version)
            echo "annotseba v${VERSION}"
            exit 0
            ;;
        -n|--dryrun)
            SNAKEMAKE_ARGS+=("--dryrun")
            shift
            ;;
        -c|--cores)
            CORES="$2"
            shift 2
            ;;
        -f|--config_file)
            SNAKEMAKE_ARGS+=("--configfile" "$(abspath "$2")")
            shift 2
            ;;
        -a|--accessions)
            SNAKEMAKE_ARGS+=("--config" "accessions_file=$(abspath "$2")")
            shift 2
            ;;
        --keep_source)
            SNAKEMAKE_ARGS+=("--config" "keep_source=true")
            shift
            ;;
        *)
            SNAKEMAKE_ARGS+=("$1")
            shift
            ;;
    esac
done

# cd to the repo root so Snakemake finds the Snakefile and all internal
# relative paths (config/config.yaml, accessions.txt, VERSION) resolve correctly.
cd "$SCRIPT_DIR"

snakemake \
    --cores "$CORES" \
    --rerun-incomplete \
    --keep-going \
    "${SNAKEMAKE_ARGS[@]}"
