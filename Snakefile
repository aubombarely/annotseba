import os

configfile: "config/config.yaml"

# ── Load accession list (species\taccession) ───────────────────────────────────
SAMPLES = []
ACC_TO_SPECIES = {}
with open(config["accessions_file"]) as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        species, acc = line.split("\t")
        SAMPLES.append((species, acc))
        ACC_TO_SPECIES[acc] = species

SPECIES  = [s for s, a in SAMPLES]
ACCESSIONS = [a for s, a in SAMPLES]
LINEAGES = config["busco_lineages"]

# ── Target rule ───────────────────────────────────────────────────────────────
rule all:
    input:
        expand("results/{species}/{acc}/genome/{acc}.gff3",
               zip, species=SPECIES, acc=ACCESSIONS),
        expand("results/{species}/{acc}/genome/{acc}_renamed.fasta",
               zip, species=SPECIES, acc=ACCESSIONS),
        expand("results/{species}/{acc}/genome/{acc}_renamed.equiv_seqID.txt",
               zip, species=SPECIES, acc=ACCESSIONS),
        expand("results/{species}/{acc}/AssemblyQC/quast/report.tsv",
               zip, species=SPECIES, acc=ACCESSIONS),
        expand("results/{species}/{acc}/AssemblyQC/assembly_stats/stats.txt",
               zip, species=SPECIES, acc=ACCESSIONS),
        [
            f"results/{s}/{a}/AssemblyQC/busco/{a}/short_summary.specific.{l}.{a}.txt"
            for s, a in SAMPLES
            for l in LINEAGES
        ],
        expand("results/{species}/{acc}/AnnotationQC/gaqet/{acc}_GAQET.stats.tsv",
               zip, species=SPECIES, acc=ACCESSIONS),
        "results/multiqc/multiqc_report.html",

# ── Download genome from NCBI ─────────────────────────────────────────────────
rule download_genome:
    output:
        fasta="results/{species}/{acc}/genome/{acc}.fna",
        gff3="results/{species}/{acc}/genome/{acc}.gff3",
    params:
        workdir="results/{species}/{acc}/genome",
        zipfile="results/{species}/{acc}/genome/ncbi_dataset.zip",
    log:
        "logs/download/{species}/{acc}.log",
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

# ── Rename FASTA sequence IDs ─────────────────────────────────────────────────
rule rename_fasta:
    input:
        fasta="results/{species}/{acc}/genome/{acc}.fna",
    output:
        fasta="results/{species}/{acc}/genome/{acc}_renamed.fasta",
        equiv="results/{species}/{acc}/genome/{acc}_renamed.equiv_seqID.txt",
    params:
        prefix=config.get("rename_prefix", "seq"),
        out_basename="results/{species}/{acc}/genome/{acc}_renamed",
        script=config.get("ncbi_fasta_rename_script", "NCBI_FastaRename"),
    log:
        "logs/rename/{species}/{acc}.log",
    shell:
        """
        {params.script} \
            -f {input.fasta} \
            -p {params.prefix} \
            -o {params.out_basename} \
            >{log} 2>&1
        """

# ── QUAST assembly statistics ─────────────────────────────────────────────────
rule run_quast:
    input:
        fasta="results/{species}/{acc}/genome/{acc}_renamed.fasta",
        gff3="results/{species}/{acc}/genome/{acc}.gff3",
    output:
        report="results/{species}/{acc}/AssemblyQC/quast/report.tsv",
    params:
        outdir="results/{species}/{acc}/AssemblyQC/quast",
        threads=config.get("threads", 4),
        min_contig=config.get("quast_min_contig", 500),
    log:
        "logs/quast/{species}/{acc}.log",
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
        fasta="results/{species}/{acc}/genome/{acc}_renamed.fasta",
    output:
        stats="results/{species}/{acc}/AssemblyQC/assembly_stats/stats.txt",
    log:
        "logs/assembly_stats/{species}/{acc}.log",
    shell:
        """
        mkdir -p $(dirname {output.stats})
        assembly-stats {input.fasta} >{output.stats} 2>{log}
        """

# ── BUSCO completeness assessment ─────────────────────────────────────────────
rule run_busco:
    input:
        fasta="results/{species}/{acc}/genome/{acc}_renamed.fasta",
    output:
        summary="results/{species}/{acc}/AssemblyQC/busco/{acc}/short_summary.specific.{lineage}.{acc}.txt",
    params:
        outdir="results/{species}/{acc}/AssemblyQC/busco",
        threads=config.get("threads", 4),
        mode="genome",
        downloads_path=config.get("busco_downloads_path", "busco_downloads"),
    log:
        "logs/busco/{species}/{acc}/{lineage}.log",
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

# ── GAQET2 annotation quality ─────────────────────────────────────────────────
rule write_gaqet_yaml:
    input:
        fasta="results/{species}/{acc}/genome/{acc}_renamed.fasta",
        gff3="results/{species}/{acc}/genome/{acc}.gff3",
    output:
        yaml="results/{species}/{acc}/AnnotationQC/gaqet/gaqet_config.yaml",
    params:
        outdir="results/{species}/{acc}/AnnotationQC/gaqet",
        threads=config.get("threads", 4),
        analyses=config.get("gaqet_analyses", ["AGAT", "BUSCO"]),
        busco_downloads=config.get("busco_downloads_path", "busco_downloads"),
        omark_db=config.get("omark_db", ""),
        detenga_db=config.get("detenga_db", ""),
    run:
        import yaml, os
        os.makedirs(params.outdir, exist_ok=True)
        cfg = {
            "ID":         wildcards.acc,
            "Assembly":   input.fasta,
            "Annotation": input.gff3,
            "Basedir":    params.outdir,
            "Threads":    params.threads,
            "Analysis":   list(params.analyses),
        }
        if "BUSCO" in params.analyses:
            cfg["BUSCO_lineages"] = params.busco_downloads
        if "OMARK" in params.analyses:
            cfg["OMARK_db"] = params.omark_db
        if "DETENGA" in params.analyses:
            cfg["DETENGA_db"] = params.detenga_db
        with open(output.yaml, "w") as fh:
            yaml.dump(cfg, fh, default_flow_style=False)

rule run_gaqet:
    input:
        yaml="results/{species}/{acc}/AnnotationQC/gaqet/gaqet_config.yaml",
    output:
        stats="results/{species}/{acc}/AnnotationQC/gaqet/{acc}_GAQET.stats.tsv",
    log:
        "logs/gaqet/{species}/{acc}.log",
    shell:
        """
        GAQET --YAML {input.yaml} >{log} 2>&1
        """

# ── MultiQC aggregate report ───────────────────────────────────────────────────
rule multiqc:
    input:
        quast=expand("results/{species}/{acc}/AssemblyQC/quast/report.tsv",
                     zip, species=SPECIES, acc=ACCESSIONS),
        busco=[
            f"results/{s}/{a}/AssemblyQC/busco/{a}/short_summary.specific.{l}.{a}.txt"
            for s, a in SAMPLES
            for l in LINEAGES
        ],
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
