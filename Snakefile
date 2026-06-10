import os

configfile: "config/config.yaml"

# ── Load accession list ────────────────────────────────────────────────────────
with open(config["accessions_file"]) as fh:
    ACCESSIONS = [l.strip() for l in fh if l.strip() and not l.startswith("#")]

LINEAGES = config["busco_lineages"]

# ── Target rule ───────────────────────────────────────────────────────────────
rule all:
    input:
        expand("results/{acc}/genome/{acc}.gff3", acc=ACCESSIONS),
        expand("results/{acc}/quast/report.tsv", acc=ACCESSIONS),
        expand("results/{acc}/assembly_stats/stats.txt", acc=ACCESSIONS),
        expand(
            "results/{acc}/busco/{acc}/short_summary.specific.{lineage}.{acc}.txt",
            acc=ACCESSIONS,
            lineage=LINEAGES,
        ),
        "results/multiqc/multiqc_report.html",

# ── Download genome from NCBI ─────────────────────────────────────────────────
rule download_genome:
    output:
        fasta="results/{acc}/genome/{acc}.fna",
        gff3="results/{acc}/genome/{acc}.gff3",
    params:
        workdir="results/{acc}/genome",
        zipfile="results/{acc}/genome/ncbi_dataset.zip",
    log:
        "logs/download/{acc}.log",
    retries: 3
    shell:
        """
        mkdir -p {params.workdir}

        datasets download genome accession {wildcards.acc} \
            --include genome,gff3 \
            --filename {params.zipfile} \
            2>{log}

        unzip -o {params.zipfile} -d {params.workdir}/tmp >>{log} 2>&1

        # Move the first FASTA found into the expected output path
        fasta=$(find {params.workdir}/tmp -name "*.fna" | head -1)
        if [ -z "$fasta" ]; then
            echo "ERROR: no .fna file found in downloaded archive" >>{log}
            exit 1
        fi
        mv "$fasta" {output.fasta}

        # Move the GFF3 annotation file
        gff=$(find {params.workdir}/tmp -name "*.gff" | head -1)
        if [ -z "$gff" ]; then
            echo "ERROR: no .gff file found in downloaded archive" >>{log}
            exit 1
        fi
        mv "$gff" {output.gff3}

        rm -rf {params.workdir}/tmp {params.zipfile}
        """

# ── QUAST assembly statistics ─────────────────────────────────────────────────
rule run_quast:
    input:
        fasta="results/{acc}/genome/{acc}.fna",
        gff3="results/{acc}/genome/{acc}.gff3",
    output:
        report="results/{acc}/quast/report.tsv",
    params:
        outdir="results/{acc}/quast",
        threads=config.get("threads", 4),
        min_contig=config.get("quast_min_contig", 500),
    log:
        "logs/quast/{acc}.log",
    shell:
        """
        quast.py {input.fasta} \
            --features {input.gff3} \
            --output-dir {params.outdir} \
            --threads {params.threads} \
            --min-contig {params.min_contig} \
            --no-gzip \
            >{log} 2>&1
        """

# ── assembly-stats ────────────────────────────────────────────────────────────
rule run_assembly_stats:
    input:
        fasta="results/{acc}/genome/{acc}.fna",
    output:
        stats="results/{acc}/assembly_stats/stats.txt",
    log:
        "logs/assembly_stats/{acc}.log",
    shell:
        """
        mkdir -p $(dirname {output.stats})
        assembly-stats {input.fasta} >{output.stats} 2>{log}
        """

# ── BUSCO completeness assessment ─────────────────────────────────────────────
rule run_busco:
    input:
        fasta="results/{acc}/genome/{acc}.fna",
    output:
        summary="results/{acc}/busco/{acc}/short_summary.specific.{lineage}.{acc}.txt",
    params:
        outdir="results/{acc}/busco",
        threads=config.get("threads", 4),
        mode="genome",
        downloads_path=config.get("busco_downloads_path", "busco_downloads"),
    log:
        "logs/busco/{acc}.log",
    shell:
        """
        busco \
            -i {input.fasta} \
            -o {wildcards.acc} \
            --out_path {params.outdir} \
            -l {wildcards.lineage} \
            -m {params.mode} \
            -c {params.threads} \
            --download_path {params.downloads_path} \
            -f \
            >{log} 2>&1
        """

# ── MultiQC aggregate report ───────────────────────────────────────────────────
rule multiqc:
    input:
        quast=expand("results/{acc}/quast/report.tsv", acc=ACCESSIONS),
        busco=expand(
            "results/{acc}/busco/{acc}/short_summary.specific.{lineage}.{acc}.txt",
            acc=ACCESSIONS,
            lineage=LINEAGES,
        ),
    output:
        report="results/multiqc/multiqc_report.html",
    params:
        outdir="results/multiqc",
    log:
        "logs/multiqc.log",
    shell:
        """
        multiqc results/ \
            --outdir {params.outdir} \
            --force \
            >{log} 2>&1
        """
