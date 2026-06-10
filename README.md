# annotseba v0.1.0

A Snakemake pipeline to download eukaryotic genomes from NCBI and assess their quality and completeness, including genome assembly statistics and annotation quality evaluation.

## Pipeline overview

```
accessions.txt  (species<TAB>accession<TAB>taxa_id)
      │
      ▼
download_genome                        (NCBI datasets CLI)
  ├── {outdir}/{species}/{acc}/genome/{acc}.fna   [temp — deleted after renaming]
  └── {outdir}/{species}/{acc}/genome/{acc}.gff3
      │
      ▼
rename_fasta                           (NCBI_FastaRename)
  ├── {outdir}/{species}/{acc}/genome/{acc}_renamed.fasta
  └── {outdir}/{species}/{acc}/genome/{acc}_renamed.equiv_seqID.txt
      │
      ├──▶ AssemblyQC/
      │      ├── run_quast          → {outdir}/{species}/{acc}/AssemblyQC/quast/
      │      ├── run_assembly_stats → {outdir}/{species}/{acc}/AssemblyQC/assembly_stats/
      │      └── run_busco          → {outdir}/{species}/{acc}/AssemblyQC/busco/
      │                                        (one run per lineage)
      ├──▶ AnnotationQC/
      │      └── write_gaqet_yaml + run_gaqet → {outdir}/{species}/{acc}/AnnotationQC/gaqet/
      │                                    │
      │                                    ▼
      │                                multiqc  → {outdir}/multiqc/multiqc_report.html
      │
      └──▶ compress  (after all QC)
             ├── {outdir}/{species}/{acc}/genome/{acc}_renamed.fasta.gz
             └── {outdir}/{species}/{acc}/genome/{acc}.gff3.gz
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
conda activate annotseba
```

> **GAQET2 must be installed separately** before running the pipeline, as its heavy dependencies (InterproScan, OMAmer, TEsorter, etc.) can conflict with other tools. Follow the [GAQET2 install docs](https://github.com/victorgcb1987/GAQET2) and ensure the `GAQET` command is available on your `$PATH` before running annotseba.

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
{outdir}/
└── {species}/
    └── {accession}/
        ├── genome/
        │   ├── {accession}.gff3                      # genome annotation
        │   ├── {accession}.gff3.gz                   # compressed after QC
        │   ├── {accession}_renamed.fasta             # renamed sequence IDs
        │   ├── {accession}_renamed.fasta.gz          # compressed after QC
        │   └── {accession}_renamed.equiv_seqID.txt   # old → new ID mapping
        │   (note: raw {accession}.fna is deleted automatically after renaming)
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
{outdir}/multiqc/
    └── multiqc_report.html                           # aggregated report
logs/
└── {rule}/
    └── {species}/
        └── {accession}.log                           # compress, busco also include lineage
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

## Third-party tools

annotseba would not be possible without the following tools. Please cite them appropriately in your work.

### Direct dependencies

| Tool | Repository | Reference |
|------|-----------|-----------|
| NCBI Datasets CLI | https://github.com/ncbi/datasets | Sayers et al. 2022, Nucleic Acids Res. |
| NCBI_FastaRename | https://github.com/aubombarely/GenoToolBox | — |
| QUAST | https://github.com/ablab/quast | Gurevich et al. 2013, Bioinformatics |
| assembly-stats | https://github.com/sanger-pathogens/assembly-stats | — |
| BUSCO | https://gitlab.com/ezlab/busco | Manni et al. 2021, Mol. Biol. Evol. |
| GAQET2 | https://github.com/victorgcb1987/GAQET2 | — |
| MultiQC | https://github.com/MultiQC/MultiQC | Ewels et al. 2016, Bioinformatics |

### GAQET2 dependencies

These tools are used internally by GAQET2 depending on which analyses are enabled:

| Tool | Repository |
|------|-----------|
| AGAT | https://github.com/NBISweden/AGAT |
| GFFread | https://github.com/gpertea/gffread |
| OMAmer | https://github.com/DessimozLab/omamer |
| OMArk | https://github.com/DessimozLab/OMArk |
| BUSCO | https://gitlab.com/ezlab/busco |
| Diamond | https://github.com/bbuchfink/diamond |
| TEsorter | https://github.com/zhangrengang/TEsorter |
| InterProScan | https://github.com/ebi-pf-team/interproscan |
| PSAURON | https://github.com/salzberg-lab/PSAURON |
| DeTEnGA | https://github.com/victorgcb1987/DeTEnGA |
