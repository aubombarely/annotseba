import os

with open("VERSION") as _vf:
    VERSION = _vf.read().strip()

print(f"annotseba v{VERSION}")

configfile: "config/config.yaml"

# ── Load accession list (species\taccession\ttaxa_id[\tprefix]) ───────────────
SAMPLES = []
ACC_TO_SPECIES = {}
ACC_TO_TAXID   = {}
ACC_TO_PREFIX  = {}
with open(config["accessions_file"]) as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        species, acc, taxa_id = fields[0], fields[1], fields[2]
        prefix = fields[3] if len(fields) > 3 and fields[3].strip().upper() != "NA" \
                 else config.get("rename_prefix", "seq")
        SAMPLES.append((species, acc))
        ACC_TO_SPECIES[acc] = species
        ACC_TO_TAXID[acc]   = taxa_id
        ACC_TO_PREFIX[acc]  = prefix

SPECIES    = [s for s, a in SAMPLES]
ACCESSIONS = [a for s, a in SAMPLES]
LINEAGES   = config["busco_lineages"]
OUTDIR     = config.get("outdir", "results")
_ks = config.get("keep_source", False)
KEEP_SOURCE = _ks if isinstance(_ks, bool) else str(_ks).lower() in ("true", "1", "yes")

# ── Target rule ───────────────────────────────────────────────────────────────
rule all:
    input:
        expand(f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.fasta",
               zip, species=SPECIES, acc=ACCESSIONS),
        expand(f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_ann.gff3",
               zip, species=SPECIES, acc=ACCESSIONS),
        expand(f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.equiv_seqID.txt",
               zip, species=SPECIES, acc=ACCESSIONS),
        expand(f"{OUTDIR}/{{species}}/{{acc}}/AssemblyQC/quast/report.tsv",
               zip, species=SPECIES, acc=ACCESSIONS),
        expand(f"{OUTDIR}/{{species}}/{{acc}}/AssemblyQC/assembly_stats/stats.txt",
               zip, species=SPECIES, acc=ACCESSIONS),
        [
            f"{OUTDIR}/{s}/{a}/AssemblyQC/busco/{a}/short_summary.specific.{l}.{a}.txt"
            for s, a in SAMPLES
            for l in LINEAGES
        ],
        expand(f"{OUTDIR}/{{species}}/{{acc}}/AnnotationQC/gaqet/{{acc}}_GAQET.stats.tsv",
               zip, species=SPECIES, acc=ACCESSIONS),
        f"{OUTDIR}/multiqc/multiqc_report.html",
        expand(f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.fasta.gz",
               zip, species=SPECIES, acc=ACCESSIONS),
        expand(f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_ann.gff3.gz",
               zip, species=SPECIES, acc=ACCESSIONS),

# ── Download genome from NCBI ─────────────────────────────────────────────────
rule download_genome:
    output:
        fasta=(f"{OUTDIR}/{{species}}/{{acc}}/genome/{{acc}}.fna"
               if KEEP_SOURCE else
               temp(f"{OUTDIR}/{{species}}/{{acc}}/genome/{{acc}}.fna")),
        gff3=(f"{OUTDIR}/{{species}}/{{acc}}/genome/{{acc}}.gff3"
              if KEEP_SOURCE else
              temp(f"{OUTDIR}/{{species}}/{{acc}}/genome/{{acc}}.gff3")),
    params:
        workdir=f"{OUTDIR}/{{species}}/{{acc}}/genome",
        zipfile=f"{OUTDIR}/{{species}}/{{acc}}/genome/ncbi_dataset.zip",
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
        fasta=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{acc}}.fna",
    output:
        fasta=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.fasta",
        equiv=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.equiv_seqID.txt",
    params:
        prefix=lambda wildcards: ACC_TO_PREFIX[wildcards.acc],
        out_basename=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb",
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

# ── Rename GFF3 sequence IDs to match renamed FASTA ──────────────────────────
rule rename_gff3:
    input:
        gff3=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{acc}}.gff3",
        equiv=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.equiv_seqID.txt",
    output:
        gff3=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_ann.gff3",
    log:
        "logs/rename_gff3/{species}/{acc}.log",
    shell:
        """
        agat_sq_rename_seqid.pl \
            --gff {input.gff3} \
            --tsv {input.equiv} \
            --output {output.gff3} \
            >{log} 2>&1
        """

# ── QUAST assembly statistics ─────────────────────────────────────────────────
rule run_quast:
    input:
        fasta=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.fasta",
        gff3=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_ann.gff3",
    output:
        report=f"{OUTDIR}/{{species}}/{{acc}}/AssemblyQC/quast/report.tsv",
    params:
        outdir=f"{OUTDIR}/{{species}}/{{acc}}/AssemblyQC/quast",
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
        fasta=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.fasta",
    output:
        stats=f"{OUTDIR}/{{species}}/{{acc}}/AssemblyQC/assembly_stats/stats.txt",
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
        fasta=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.fasta",
    output:
        summary=f"{OUTDIR}/{{species}}/{{acc}}/AssemblyQC/busco/{{acc}}/short_summary.specific.{{lineage}}.{{acc}}.txt",
    params:
        outdir=f"{OUTDIR}/{{species}}/{{acc}}/AssemblyQC/busco",
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
        fasta=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.fasta",
        gff3=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_ann.gff3",
    output:
        yaml=f"{OUTDIR}/{{species}}/{{acc}}/AnnotationQC/gaqet/gaqet_config.yaml",
    params:
        outdir=f"{OUTDIR}/{{species}}/{{acc}}/AnnotationQC/gaqet",
        threads=config.get("threads", 4),
        analyses=config.get("gaqet_analyses", ["AGAT", "BUSCO"]),
        busco_downloads=config.get("busco_downloads_path", "busco_downloads"),
        omark_db=config.get("omark_db", ""),
        detenga_db=config.get("detenga_db", ""),
        taxa_id=lambda wildcards: ACC_TO_TAXID[wildcards.acc],
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
            cfg["OMARK_db"]  = params.omark_db
            cfg["taxid"]     = params.taxa_id
        if "DETENGA" in params.analyses:
            cfg["DETENGA_db"] = params.detenga_db
        with open(output.yaml, "w") as fh:
            yaml.dump(cfg, fh, default_flow_style=False)

rule run_gaqet:
    input:
        yaml=f"{OUTDIR}/{{species}}/{{acc}}/AnnotationQC/gaqet/gaqet_config.yaml",
    output:
        stats=f"{OUTDIR}/{{species}}/{{acc}}/AnnotationQC/gaqet/{{acc}}_GAQET.stats.tsv",
    log:
        "logs/gaqet/{species}/{acc}.log",
    shell:
        """
        GAQET --YAML {input.yaml} >{log} 2>&1
        """

# ── Compress genome files after all QC is done ────────────────────────────────
rule compress:
    input:
        fasta=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.fasta",
        gff3=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_ann.gff3",
        # Wait for all QC outputs before compressing
        quast=f"{OUTDIR}/{{species}}/{{acc}}/AssemblyQC/quast/report.tsv",
        assembly_stats=f"{OUTDIR}/{{species}}/{{acc}}/AssemblyQC/assembly_stats/stats.txt",
        busco=[
            f"{OUTDIR}/{{species}}/{{acc}}/AssemblyQC/busco/{{acc}}/short_summary.specific.{l}.{{acc}}.txt"
            for l in LINEAGES
        ],
        gaqet=f"{OUTDIR}/{{species}}/{{acc}}/AnnotationQC/gaqet/{{acc}}_GAQET.stats.tsv",
    output:
        fasta_gz=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.fasta.gz",
        gff3_gz=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_ann.gff3.gz",
    log:
        "logs/compress/{species}/{acc}.log",
    shell:
        """
        gzip -k {input.fasta} >{log} 2>&1
        gzip -k {input.gff3}  >>{log} 2>&1
        """

# ── MultiQC aggregate report ───────────────────────────────────────────────────
rule multiqc:
    input:
        quast=expand(f"{OUTDIR}/{{species}}/{{acc}}/AssemblyQC/quast/report.tsv",
                     zip, species=SPECIES, acc=ACCESSIONS),
        busco=[
            f"{OUTDIR}/{s}/{a}/AssemblyQC/busco/{a}/short_summary.specific.{l}.{a}.txt"
            for s, a in SAMPLES
            for l in LINEAGES
        ],
    output:
        report=f"{OUTDIR}/multiqc/multiqc_report.html",
    params:
        outdir=f"{OUTDIR}/multiqc",
    log:
        "logs/multiqc.log",
    shell:
        """
        multiqc {OUTDIR}/ \
            --outdir {params.outdir} \
            --force \
            >{log} 2>&1
        """
