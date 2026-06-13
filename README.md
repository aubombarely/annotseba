# annotseba v0.2.0

A Snakemake pipeline to download eukaryotic genomes from NCBI and assess their quality and completeness, including genome assembly statistics and annotation quality evaluation.

## Pipeline overview

```
accessions.txt  (species<TAB>accession<TAB>taxa_id[<TAB>prefix])
      │
      ▼
download_genome                        (NCBI datasets CLI)
  ├── {outdir}/{species}/{acc}/genome/{acc}.fna   [temp unless --keep_source]
  └── {outdir}/{species}/{acc}/genome/{acc}.gff3  [temp unless --keep_source]
      │
      │  If NCBI returns no FASTA or GFF3, empty placeholder files are
      │  created and the accession is skipped in all downstream rules.
      │
      ▼
rename_fasta                           (NCBI_FastaRename)
  ├── {outdir}/{species}/{acc}/genome/{species}_{acc}_asb.fasta
  └── {outdir}/{species}/{acc}/genome/{species}_{acc}_asb.equiv_seqID.txt
      │
      ▼
rename_gff3                            (AGAT agat_sq_rename_seqid.pl)
  └── {outdir}/{species}/{acc}/genome/{species}_{acc}_asb.gff3
      │
      ├──▶ AssemblyQC/
      │      ├── run_quast          → {outdir}/{species}/{acc}/AssemblyQC/quast/
      │      ├── run_assembly_stats → {outdir}/{species}/{acc}/AssemblyQC/assembly_stats/
      │      └── run_busco          → {outdir}/{species}/{acc}/AssemblyQC/busco/
      │                                        (one run per lineage in busco_lineages)
      ├──▶ AnnotationQC/
      │      └── write_gaqet_yaml + run_gaqet → {outdir}/{species}/{acc}/AnnotationQC/gaqet/
      │
      └──▶ compress  (after all QC is complete)
             ├── {outdir}/{species}/{acc}/genome/{species}_{acc}_asb.fasta.gz
             └── {outdir}/{species}/{acc}/genome/{species}_{acc}_asb.gff3.gz
                    │
                    ▼
             merge_gaqet_stats  → {outdir}/report/all_species_GAQET.stats.tsv
                    │
                    ▼
             run_gaqet_plot     → {outdir}/report/all_species_GAQET.plot.{fmt}
                    │
                    ▼
             generate_report    → {outdir}/report/annotseba_AssemblyQC.tsv
                                  {outdir}/report/annotseba_AnnotationQC.tsv
                                  {outdir}/report/annotseba_report.html
                    │
                    ▼
             compute_usage      → {outdir}/report/computer_usage.log
```

## Tools used

| Tool | Purpose |
|------|---------|
| [NCBI datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) | Download genome FASTA and GFF3 annotation |
| [NCBI_FastaRename](https://github.com/aubombarely/GenoToolBox/blob/master/AnnotThis/NCBI_FastaRename) | Rename sequence IDs with a consistent prefix |
| [AGAT](https://github.com/NBISweden/AGAT) | Rename GFF3 sequence IDs to match the renamed FASTA |
| [QUAST](https://quast.sourceforge.net/) | Assembly statistics (N50, contigs, gene features via GFF3) |
| [assembly-stats](https://github.com/sanger-pathogens/assembly-stats) | Lightweight assembly statistics |
| [BUSCO](https://busco.ezlab.org/) | Genome completeness against conserved gene sets (multiple lineages) |
| [GAQET2](https://github.com/victorgcb1987/GAQET2) | Genome annotation quality evaluation |
| matplotlib / pandas | Custom HTML+TSV report generation |

## Installation

```bash
conda env create -f envs/pipeline.yaml
conda activate annotseba
```

> **GAQET2 must be installed separately** before running the pipeline, as its heavy dependencies (InterproScan, OMAmer, TEsorter, etc.) can conflict with other tools. Follow the [GAQET2 install docs](https://github.com/victorgcb1987/GAQET2) and ensure the `GAQET` and `GAQET_PLOT` commands are available on your `$PATH`.

Ensure `NCBI_FastaRename` is on your `$PATH` or set its full path in `config/config.yaml`.

## Usage

**1. Add your entries to `accessions.txt`** (tab-separated, one per line):
```
Homo_sapiens          GCA_000001405.29  9606  Hs
Arabidopsis_thaliana  GCA_000001735.4   3702  At
Oryza_sativa          GCA_001433935.1   4530  NA
```

| Column | Required | Description |
|--------|----------|-------------|
| `species` | Yes | Species name — no spaces, use underscores |
| `accession` | Yes | NCBI accession (GCA_* or GCF_*) |
| `taxa_id` | Yes | NCBI Taxonomy ID — required for OMARK analysis in GAQET2 |
| `prefix` | No | Prefix for renamed sequence IDs (e.g. `Hs`, `At`) — use `NA` or omit to fall back to `rename_prefix` in `config.yaml` |

**2. Edit `config/config.yaml`** — key settings to review before running:

- BUSCO lineages for your organisms:
```yaml
busco_lineages:
  - embryophyta_odb10
  - viridiplantae_odb10
```

- GAQET2 analyses to run:
```yaml
gaqet_analyses:
  - AGAT
  - BUSCO
  - PSAURON
  - DETENGA
  - OMARK
  - PROTHOMOLOGY
```

- Protein homology databases (required only for PROTHOMOLOGY):
```yaml
prothomology_dbs:
  TREMBL:     "/path/to/uniprot_trembl.dmnd"
  SWISSPROT:  "/path/to/uniprot_sprot.dmnd"
```

**3. Dry-run to check the execution plan:**
```bash
bash run_annotseba.sh --dryrun
```

**4. Run the pipeline:**
```bash
bash run_annotseba.sh --cores 8
```

You can also supply a separate config file or accessions file without editing the defaults:
```bash
bash run_annotseba.sh --config_file my_project.yaml --accessions my_species.tsv --cores 16
```

## run_annotseba.sh options

| Option | Description |
|--------|-------------|
| `-h, --help` | Show help and exit |
| `-v, --version` | Show version and exit |
| `-n, --dryrun` | Show execution plan without running |
| `-c, --cores INT` | Number of CPU cores (default: `$SNAKEMAKE_CORES` or 8) |
| `-f, --config_file FILE` | Config YAML merged on top of `config/config.yaml` (later values win) |
| `-a, --accessions FILE` | Path to accessions TSV (overrides `config.yaml`) |
| `--keep_source` | Keep raw NCBI files (`{acc}.fna`, `{acc}.gff3`) after renaming |

Any unrecognised options are passed directly to Snakemake (e.g. `--forceall`, `--until <rule>`).

## Output structure

```
{outdir}/
└── {species}/
    └── {accession}/
        ├── genome/
        │   ├── {species}_{accession}_asb.fasta             # assembly with renamed seq IDs
        │   ├── {species}_{accession}_asb.fasta.gz          # compressed after QC
        │   ├── {species}_{accession}_asb.gff3              # annotation with renamed seq IDs
        │   ├── {species}_{accession}_asb.gff3.gz           # compressed after QC
        │   ├── {species}_{accession}_asb.equiv_seqID.txt   # old → new seq ID mapping
        │   ├── {accession}.fna                             # raw NCBI FASTA (--keep_source only)
        │   └── {accession}.gff3                            # raw NCBI GFF3  (--keep_source only)
        ├── AssemblyQC/
        │   ├── quast/                                # QUAST assembly report
        │   ├── assembly_stats/                       # assembly-stats output
        │   └── busco/
        │       └── {accession}/
        │           └── short_summary.specific.{lineage}.{accession}.txt
        └── AnnotationQC/
            └── gaqet/
                ├── gaqet_config.yaml
                └── {accession}_GAQET.stats.tsv
{outdir}/report/
    ├── all_species_GAQET.stats.tsv       # merged GAQET stats (all accessions)
    ├── all_species_GAQET.plot.{fmt}      # GAQET_PLOT output
    ├── annotseba_AssemblyQC.tsv          # assembly QC summary table
    ├── annotseba_AnnotationQC.tsv        # annotation QC summary table
    ├── annotseba_report.html             # self-contained HTML report
    └── computer_usage.log               # disk usage + carbon footprint estimate
benchmarks/
└── {rule}/{species}/{accession}.tsv      # Snakemake benchmark (CPU/wall time, memory)
logs/
└── {rule}/{species}/{accession}.log
```

## HTML report

The pipeline produces a self-contained HTML report (`annotseba_report.html`) with no external dependencies. It contains:

- **Assembly QC table** — sortable table with QUAST metrics (N50, L50, contig counts, GC%, etc.) and BUSCO completeness (Complete / Single / Duplicated / Fragmented / Missing) for each configured lineage.
- **Assembly size barplot** — total assembly size (Mb) per species.
- **Cumulative contig length plot** — QUAST-style cumulative plot for all species.
- **Annotation QC table** — sortable table with all GAQET metrics, grouped by analysis type with colour-coded column headers:

  | Group | Colour | Columns |
  |-------|--------|---------|
  | General | grey | Species, Accession, NCBI_TaxID, Assembly_Version, Annotation_Version |
  | AGAT | blue | Gene model counts, lengths, UTR stats, model completeness |
  | BUSCO | green | `Annotation_BUSCO_*` |
  | PSAURON | orange | `PSAURON SCORE` |
  | DETENGA | red | `DETENGA_FPV`, `DETENGA_FP%` |
  | OMARK | purple | `OMArk Consistency/Completeness/Species Composition` |
  | PROTHOMOLOGY | brown | `ProteinsWithTREMBLHits (%)`, `ProteinsWithSWISSPROTHits (%)` |

- **GAQET plot** — image produced by `GAQET_PLOT` across all species, embedded directly in the HTML.

## Resource & carbon usage log

After the report is generated, the `compute_usage` rule writes `{outdir}/report/computer_usage.log` with:

- **Disk usage** — total space used by the output directory (`du`).
- **Compute time** — total CPU time and wall time aggregated from Snakemake benchmark files.
- **Carbon footprint estimate** — energy (kWh) and CO₂ equivalent (g) based on:

  ```
  energy (kWh) = total_CPU_time (h) × cpu_tdp_per_core (W) / 1000
  CO2 (g)      = energy × carbon_intensity (gCO2/kWh)
  ```

  Default values use the IEA 2022 world average (475 gCO₂/kWh, 10 W/core). Tune both in `config.yaml` for your hardware and location (see [electricitymaps.com](https://app.electricitymaps.com)).

## Graceful degradation

If NCBI returns no FASTA or GFF3 for an accession (e.g. the record exists but has no genomic sequence available), the pipeline:

1. Logs a warning to `logs/download/{species}/{accession}.log`.
2. Creates empty placeholder files so Snakemake considers the rule successful.
3. Skips all downstream rules for that accession (rename, QC, compression).
4. Excludes the accession from the final report silently.

All other accessions continue to completion and are included in the report normally.

## Configuration reference

All settings live in `config/config.yaml`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `accessions_file` | `accessions.txt` | Path to accession list |
| `outdir` | `results` | Root directory for all output files |
| `rename_prefix` | `Sp` | Default prefix for renamed sequence IDs |
| `ncbi_fasta_rename_script` | `NCBI_FastaRename` | Path to the rename script |
| `busco_lineages` | `[embryophyta_odb10]` | List of BUSCO lineages to run |
| `busco_downloads_path` | `busco_downloads/` | Local cache for BUSCO lineage databases |
| `quast_min_contig` | `500` | Minimum contig length (bp) for QUAST |
| `gaqet_analyses` | `[AGAT, BUSCO]` | GAQET2 analyses to run |
| `omark_db` | `""` | OMAmer database path (required for OMARK) |
| `detenga_db` | `""` | DeTEnGA database path (required for DETENGA) |
| `prothomology_dbs` | `{}` | Protein databases as `TAG: path` pairs (required for PROTHOMOLOGY) |
| `gaqet_plot_format` | `png` | Output format for GAQET_PLOT (`png`, `jpeg`, `svg`, `pdf`) |
| `keep_source` | `false` | Keep raw NCBI `{acc}.fna` / `{acc}.gff3` after renaming |
| `carbon_intensity` | `475` | Carbon intensity in gCO₂/kWh for footprint estimate |
| `cpu_tdp_per_core` | `10` | CPU TDP per core in watts for footprint estimate |
| `threads` | `8` | Threads per job |

### BUSCO lineages reference

| Organism group | Lineages |
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
| AGAT | https://github.com/NBISweden/AGAT | Dainat et al. 2021, JOSS |
| QUAST | https://github.com/ablab/quast | Gurevich et al. 2013, Bioinformatics |
| assembly-stats | https://github.com/sanger-pathogens/assembly-stats | — |
| BUSCO | https://gitlab.com/ezlab/busco | Manni et al. 2021, Mol. Biol. Evol. |
| GAQET2 | https://github.com/victorgcb1987/GAQET2 | — |

### GAQET2 dependencies

These tools are used internally by GAQET2 depending on which analyses are enabled:

| Tool | Repository |
|------|-----------|
| GFFread | https://github.com/gpertea/gffread |
| OMAmer | https://github.com/DessimozLab/omamer |
| OMArk | https://github.com/DessimozLab/OMArk |
| BUSCO | https://gitlab.com/ezlab/busco |
| Diamond | https://github.com/bbuchfink/diamond |
| TEsorter | https://github.com/zhangrengang/TEsorter |
| InterProScan | https://github.com/ebi-pf-team/interproscan |
| PSAURON | https://github.com/salzberg-lab/PSAURON |
| DeTEnGA | https://github.com/victorgcb1987/DeTEnGA |
