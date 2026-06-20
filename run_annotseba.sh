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
    -h, --help              Show this help message and exit
    -v, --version           Show version and exit
    -n, --dryrun            Dry-run: show what would be executed without running
    -c, --cores INT         Number of CPU cores to use (default: \$SNAKEMAKE_CORES or 8)
    -f, --config_file FILE  Path to a config YAML (merged on top of default config/config.yaml)
    -a, --accessions FILE   Path to accessions TSV file (overrides config.yaml)
    --keep_source           Keep raw NCBI files ({acc}.fna and {acc}.gff3) after renaming
    -r, --rundir DIR        Top-level run directory (default: annotseba_run); v1 only
                            results/, logs/ and benchmarks/ are created inside

Phytolaeno v2 integration:
    --byspecies_root DIR    Path to a Phytolaeno BySpecies/ root. Outputs land inside
                            the pre-existing species tree (v2 layout).
    --accession ACC         Single NCBI accession (e.g. GCF_000001405.40).
                            Use with --species instead of an accessions TSV file.
    --species SPECIES       Species name matching generate_species_dir.py (e.g. Homo_sapiens).
                            Required when --accession is used.
    --taxid TAXID           NCBI taxid (optional; auto-read from taxonomy.json if omitted).
    --run_gaqet             Enable GAQET annotation quality evaluation
    --no_gaqet              Disable GAQET annotation quality evaluation

Snakemake pass-through:
    Any additional arguments are forwarded directly to Snakemake.
    e.g. --forceall, --until <rule>, --rerun-triggers mtime

Examples:
    # Standard v1 run (accessions.txt in current directory)
    bash run_annotseba.sh --cores 16
    bash run_annotseba.sh --dryrun
    bash run_annotseba.sh --accessions my_species.tsv --cores 16

    # Phytolaeno v2 integration — single accession
    bash run_annotseba.sh --byspecies_root /data/BySpecies \\
        --accession GCF_000001405.40 --species Homo_sapiens --cores 16

    # Phytolaeno v2 integration — batch file
    bash run_annotseba.sh --byspecies_root /data/BySpecies \\
        --accessions batch.tsv --cores 32

Input:
    Edit accessions.txt and config/config.yaml before running.
    See README.md for full documentation.
EOF
}

# Convert a path to absolute using the caller's working directory.
# Must be called before any cd; keeps relative paths usable.
abspath() { [[ "$1" = /* ]] && echo "$1" || echo "$PWD/$1"; }

CORES="${SNAKEMAKE_CORES:-8}"
SNAKEMAKE_ARGS=()
_SINGLE_ACCESSION=""
_SINGLE_SPECIES=""
_SINGLE_TAXID="NA"
_TMP_ACC_FILE=""

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
        -r|--rundir)
            SNAKEMAKE_ARGS+=("--config" "rundir=$2")
            shift 2
            ;;
        --byspecies_root)
            SNAKEMAKE_ARGS+=("--config" "byspecies_root=$(abspath "$2")")
            shift 2
            ;;
        --accession)
            _SINGLE_ACCESSION="$2"
            shift 2
            ;;
        --species)
            _SINGLE_SPECIES="$2"
            shift 2
            ;;
        --taxid)
            _SINGLE_TAXID="$2"
            shift 2
            ;;
        --run_gaqet)
            SNAKEMAKE_ARGS+=("--config" "run_gaqet=true")
            shift
            ;;
        --no_gaqet)
            SNAKEMAKE_ARGS+=("--config" "run_gaqet=false")
            shift
            ;;
        *)
            SNAKEMAKE_ARGS+=("$1")
            shift
            ;;
    esac
done

# Single-accession mode: build a temp accessions file so the Snakefile
# sees the standard format without requiring a pre-existing TSV.
if [[ -n "$_SINGLE_ACCESSION" ]]; then
    if [[ -z "$_SINGLE_SPECIES" ]]; then
        echo "ERROR: --species is required when --accession is used" >&2
        exit 1
    fi
    _TMP_ACC_FILE=$(mktemp /tmp/annotseba_acc_XXXXXX.tsv)
    printf '%s\t%s\t%s\n' "$_SINGLE_SPECIES" "$_SINGLE_ACCESSION" "$_SINGLE_TAXID" \
        > "$_TMP_ACC_FILE"
    SNAKEMAKE_ARGS+=("--config" "accessions_file=$_TMP_ACC_FILE")
    echo "annotseba: single-accession mode — ${_SINGLE_SPECIES} / ${_SINGLE_ACCESSION}" \
         "(taxid=${_SINGLE_TAXID})"
fi

# Remove temp accessions file on exit
cleanup() {
    [[ -n "$_TMP_ACC_FILE" && -f "$_TMP_ACC_FILE" ]] && rm -f "$_TMP_ACC_FILE"
}
trap cleanup EXIT

snakemake \
    --snakefile "$SCRIPT_DIR/Snakefile" \
    --cores "$CORES" \
    --rerun-incomplete \
    --keep-going \
    "${SNAKEMAKE_ARGS[@]}"
