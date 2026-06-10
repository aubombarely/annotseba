# annotseba

A Snakemake pipeline to download eukaryotic genomes from NCBI and assess their quality and completeness, including genome assembly statistics and annotation quality evaluation.

## Pipeline overview

```
accessions.txt  (species<TAB>accession)
      │
      ▼
download_genome                        (NCBI datasets CLI)
  ├── results/{species}/{acc}/genome/{acc}.fna
  └── results/{species}/{acc}/genome/{acc}.gff3
      │
      ▼
rename_fasta                           (NCBI_FastaRename)
  ├── results/{species}/{acc}/genome/{acc}_renamed.fasta
  └── results/{species}/{acc}/genome/{acc}_renamed.equiv_seqID.txt
      │
      ├──▶ AssemblyQC/
      │      ├── run_quast          → results/{species}/{acc}/AssemblyQC/quast/
      │      ├── run_assembly_stats → results/{species}/{acc}/AssemblyQC/assembly_stats/
      │      └── run_busco          → results/{species}/{acc}/AssemblyQC/busco/
      │                                        (one run per lineage)
      └──▶ AnnotationQC/
             └── write_gaqet_yaml + run_gaqet → results/{species}/{acc}/AnnotationQC/gaqet/
                                    │
                                    ▼
                                multiqc  → results/multiqc/multiqc_report.html
```

## Tools used

| Tool | Purpose |
|------|---------|
| [NCBI datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) | Download genome FASTA and GFF3 annotation |
| [NCBI_FastaRename](https://github.com/aubombarely/GenoToolBox/blob/master/AnnotThis/NCBI_FastaRename) | Rename sequence IDs with a consistent prefix |
| [QUAST](https://quast.sourceforge.net/) | Assembly statistics (N50, contigs, gene features via GFF3) |
| [assembly-stats](https://github.com/sanger-pathogens/assembly-stats) | Lightweight assembly statistics |
| [BUSCO](https://busco.ezlab.org/) | Genome completeness against conserved gene sets (multiple lineages supported) |
| [GAQET2](https://github.com/victorgcb1987/GAQET2) | Genome annotation quality evaluation |
| [MultiQC](https://multiqc.info/) | Aggregate all results into a single HTML report |

## Installation

```bash
conda env create -f envs/pipeline.yaml
conda activate genome_pipeline
```

> **Note:** GAQET2 has additional dependencies (InterproScan, OMAmer, etc.) that may require manual installation. See the [GAQET2 install docs](https://github.com/victorgcb1987/GAQET2) for details.

Ensure `NCBI_FastaRename` is on your `$PATH` or set its full path in `config/config.yaml`.

## Usage

**1. Add your entries to `accessions.txt`** (tab-separated `species<TAB>accession<TAB>taxa_id`, one per line):
```
Homo_sapiens	GCA_000001405.29	9606
Arabidopsis_thaliana	GCA_000001735.4	3702
```
Species names should not contain spaces — use underscores. The `taxa_id` is the NCBI Taxonomy ID, required when using the OMARK analysis in GAQET2.

**2. Edit `config/config.yaml`** — key settings to review before running:

- BUSCO lineages for your organisms:
```yaml
busco_lineages:
  - embryophyta_odb10
  - viridiplantae_odb10
```

- Sequence ID rename prefix:
```yaml
rename_prefix: "Sp"
```

- GAQET2 analyses to run:
```yaml
gaqet_analyses:
  - AGAT
  - BUSCO
```

**3. Dry-run to check the execution plan:**
```bash
snakemake --dryrun --cores 8
```

**4. Run the pipeline:**
```bash
bash run_annotseba.sh --cores 8
```

## Output structure

```
results/
└── {species}/
    └── {accession}/
        ├── genome/
        │   ├── {accession}.fna                       # raw genome FASTA
        │   ├── {accession}.gff3                      # genome annotation
        │   ├── {accession}_renamed.fasta             # renamed sequence IDs
        │   └── {accession}_renamed.equiv_seqID.txt   # old → new ID mapping
        ├── AssemblyQC/
        │   ├── quast/                                # QUAST assembly report
        │   ├── assembly_stats/                       # assembly-stats output
        │   └── busco/                                # one subdirectory per lineage
        │       └── {accession}/
        │           └── short_summary.specific.{lineage}.{accession}.txt
        └── AnnotationQC/
            └── gaqet/                                # GAQET2 annotation QC
                ├── gaqet_config.yaml
                └── {accession}_GAQET.stats.tsv
results/multiqc/
    └── multiqc_report.html                           # aggregated report
logs/
└── {rule}/
    └── {species}/
        └── {accession}.log
```

## Configuration

All settings are in `config/config.yaml`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `accessions_file` | `accessions.txt` | Path to accession list |
| `outdir` | `results` | Root directory for all output files |
| `rename_prefix` | `Sp` | Prefix for renamed sequence IDs |
| `ncbi_fasta_rename_script` | `NCBI_FastaRename` | Path to the rename script |
| `busco_lineages` | `[embryophyta_odb10]` | List of BUSCO lineages to run |
| `busco_downloads_path` | `busco_downloads/` | Local cache for BUSCO databases |
| `quast_min_contig` | `500` | Minimum contig length for QUAST |
| `gaqet_analyses` | `[AGAT, BUSCO]` | GAQET2 analyses to run |
| `omark_db` | `""` | OMAmer database path (required for OMARK analysis) |
| `detenga_db` | `""` | DeTEnGA database path (required for DETENGA analysis) |
| `threads` | `8` | Threads per job |

### BUSCO lineages reference

| Organism group | Lineage |
|----------------|---------|
| Plants | `embryophyta_odb10`, `viridiplantae_odb10` |
| Animals | `metazoa_odb10`, `vertebrata_odb10`, `mammalia_odb10`, `insecta_odb10` |
| Fungi | `fungi_odb10`, `ascomycota_odb10`, `basidiomycota_odb10` |

Full list: https://busco.ezlab.org/list_of_lineages.html
