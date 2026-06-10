# annotseba

A Snakemake pipeline to download eukaryotic genomes from NCBI and assess their quality and completeness.

## Pipeline overview

```
accessions.txt
      │
      ▼
download_genome        (NCBI datasets CLI)
  ├── {acc}.fna
  └── {acc}.gff3
      │
      ├──▶ run_quast          → results/{acc}/quast/
      ├──▶ run_assembly_stats → results/{acc}/assembly_stats/
      └──▶ run_busco          → results/{acc}/busco/   (one run per lineage)
                                        │
                                        ▼
                                    multiqc  → results/multiqc/multiqc_report.html
```

## Tools used

| Tool | Purpose |
|------|---------|
| [NCBI datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) | Download genome FASTA and GFF3 annotation |
| [QUAST](https://quast.sourceforge.net/) | Assembly statistics (N50, contigs, gene features) |
| [assembly-stats](https://github.com/sanger-pathogens/assembly-stats) | Lightweight assembly statistics |
| [BUSCO](https://busco.ezlab.org/) | Genome completeness against conserved gene sets |
| [MultiQC](https://multiqc.info/) | Aggregate all results into a single HTML report |

## Installation

```bash
conda env create -f envs/pipeline.yaml
conda activate genome_pipeline
```

## Usage

**1. Add your accession IDs to `accessions.txt`** (one GCA_* or GCF_* per line):
```
GCA_000001405.29
GCA_000001735.4
```

**2. Edit `config/config.yaml`** — at minimum set the BUSCO lineage(s) for your organisms:
```yaml
busco_lineages:
  - embryophyta_odb10
  - viridiplantae_odb10
```

Common lineages:

| Organism group | Lineage |
|----------------|---------|
| Plants | `embryophyta_odb10`, `viridiplantae_odb10` |
| Animals | `metazoa_odb10`, `vertebrata_odb10`, `mammalia_odb10`, `insecta_odb10` |
| Fungi | `fungi_odb10`, `ascomycota_odb10`, `basidiomycota_odb10` |

Full list: https://busco.ezlab.org/list_of_lineages.html

**3. Dry-run to check the execution plan:**
```bash
snakemake --dryrun --cores 8
```

**4. Run the pipeline:**
```bash
bash run_pipeline.sh --cores 8
```

## Output structure

```
results/
└── {accession}/
    ├── genome/
    │   ├── {accession}.fna        # genome FASTA
    │   └── {accession}.gff3       # genome annotation
    ├── quast/                     # QUAST assembly report
    ├── assembly_stats/            # assembly-stats output
    └── busco/                     # one subdirectory per lineage
        └── {accession}/
            └── short_summary.specific.{lineage}.{accession}.txt
results/multiqc/
    └── multiqc_report.html        # aggregated report
```

## Configuration

All settings are in `config/config.yaml`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `accessions_file` | `accessions.txt` | Path to accession list |
| `busco_lineages` | `[embryophyta_odb10]` | List of BUSCO lineages to run |
| `busco_downloads_path` | `busco_downloads/` | Local cache for BUSCO databases |
| `quast_min_contig` | `500` | Minimum contig length for QUAST |
| `threads` | `8` | Threads per job |
