import os

# workflow.basedir is the directory containing this Snakefile — always
# absolute regardless of where the user runs the pipeline from.
_basedir = workflow.basedir

with open(os.path.join(_basedir, "VERSION")) as _vf:
    VERSION = _vf.read().strip()

print(f"annotseba v{VERSION}")

configfile: os.path.join(_basedir, "config", "config.yaml")

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
KEEP_SOURCE     = _ks if isinstance(_ks, bool) else str(_ks).lower() in ("true", "1", "yes")
GAQET_PLOT_FMT  = config.get("gaqet_plot_format", "png")

# All benchmark TSVs produced by Snakemake — aggregated for compute_usage rule
BENCHMARK_FILES = (
    expand("benchmarks/download/{species}/{acc}.tsv",       zip, species=SPECIES, acc=ACCESSIONS) +
    expand("benchmarks/rename/{species}/{acc}.tsv",         zip, species=SPECIES, acc=ACCESSIONS) +
    expand("benchmarks/rename_gff3/{species}/{acc}.tsv",    zip, species=SPECIES, acc=ACCESSIONS) +
    expand("benchmarks/quast/{species}/{acc}.tsv",          zip, species=SPECIES, acc=ACCESSIONS) +
    expand("benchmarks/assembly_stats/{species}/{acc}.tsv", zip, species=SPECIES, acc=ACCESSIONS) +
    [f"benchmarks/busco/{s}/{a}/{l}.tsv" for s, a in SAMPLES for l in LINEAGES] +
    expand("benchmarks/gaqet/{species}/{acc}.tsv",          zip, species=SPECIES, acc=ACCESSIONS) +
    expand("benchmarks/compress/{species}/{acc}.tsv",       zip, species=SPECIES, acc=ACCESSIONS)
)

# ── Target rule ───────────────────────────────────────────────────────────────
rule all:
    input:
        expand(f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.fasta",
               zip, species=SPECIES, acc=ACCESSIONS),
        expand(f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.gff3",
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
        expand(f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.fasta.gz",
               zip, species=SPECIES, acc=ACCESSIONS),
        expand(f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.gff3.gz",
               zip, species=SPECIES, acc=ACCESSIONS),
        f"{OUTDIR}/report/annotseba_AssemblyQC.tsv",
        f"{OUTDIR}/report/annotseba_AnnotationQC.tsv",
        f"{OUTDIR}/report/annotseba_report.html",
        f"{OUTDIR}/report/computer_usage.log",

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
    benchmark:
        "benchmarks/download/{species}/{acc}.tsv"
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
            echo "WARNING: no .fna file found for {wildcards.acc} — creating empty placeholder" >>{log}
            touch {output.fasta} {output.gff3}
            rm -rf {params.workdir}/tmp {params.zipfile}
            exit 0
        fi
        mv "$fasta" {output.fasta}

        # Move the GFF3 annotation file
        gff=$(find {params.workdir}/tmp -name "*.gff" | head -1)
        if [ -z "$gff" ]; then
            echo "WARNING: no .gff file found for {wildcards.acc} — creating empty placeholder" >>{log}
            touch {output.gff3}
        else
            mv "$gff" {output.gff3}
        fi

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
    benchmark:
        "benchmarks/rename/{species}/{acc}.tsv"
    log:
        "logs/rename/{species}/{acc}.log",
    shell:
        """
        if [ ! -s {input.fasta} ]; then
            echo "SKIP ({wildcards.acc}): download failed, skipping rename" >{log}
            touch {output.fasta} {output.equiv}
            exit 0
        fi
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
        gff3=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.gff3",
    benchmark:
        "benchmarks/rename_gff3/{species}/{acc}.tsv"
    log:
        "logs/rename_gff3/{species}/{acc}.log",
    shell:
        """
        if [ ! -s {input.gff3} ]; then
            echo "SKIP ({wildcards.acc}): download failed, skipping GFF3 rename" >{log}
            touch {output.gff3}
            exit 0
        fi
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
        gff3=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.gff3",
    output:
        report=f"{OUTDIR}/{{species}}/{{acc}}/AssemblyQC/quast/report.tsv",
    params:
        outdir=f"{OUTDIR}/{{species}}/{{acc}}/AssemblyQC/quast",
        threads=config.get("threads", 4),
        min_contig=config.get("quast_min_contig", 500),
    benchmark:
        "benchmarks/quast/{species}/{acc}.tsv"
    log:
        "logs/quast/{species}/{acc}.log",
    shell:
        """
        if [ ! -s {input.fasta} ]; then
            echo "SKIP ({wildcards.acc}): download failed, skipping QUAST" >{log}
            mkdir -p {params.outdir}
            touch {output.report}
            exit 0
        fi
        quast.py {input.fasta} \
            --features {input.gff3} \
            --output-dir {params.outdir} \
            --threads {params.threads} \
            --min-contig {params.min_contig} \
            >{log} 2>&1
        """

# ── assembly-stats ────────────────────────────────────────────────────────────
rule run_assembly_stats:
    input:
        fasta=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.fasta",
    output:
        stats=f"{OUTDIR}/{{species}}/{{acc}}/AssemblyQC/assembly_stats/stats.txt",
    benchmark:
        "benchmarks/assembly_stats/{species}/{acc}.tsv"
    log:
        "logs/assembly_stats/{species}/{acc}.log",
    shell:
        """
        mkdir -p $(dirname {output.stats})
        if [ ! -s {input.fasta} ]; then
            echo "SKIP ({wildcards.acc}): download failed, skipping assembly-stats" >{log}
            touch {output.stats}
            exit 0
        fi
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
    benchmark:
        "benchmarks/busco/{species}/{acc}/{lineage}.tsv"
    log:
        "logs/busco/{species}/{acc}/{lineage}.log",
    shell:
        """
        if [ ! -s {input.fasta} ]; then
            echo "SKIP ({wildcards.acc}): download failed, skipping BUSCO" >{log}
            mkdir -p $(dirname {output.summary})
            touch {output.summary}
            exit 0
        fi
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
        gff3=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.gff3",
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
        if os.path.getsize(str(input.fasta)) == 0:
            open(str(output.yaml), "w").close()
        else:
            cfg = {
                "ID":         wildcards.acc,
                "Assembly":   input.fasta,
                "Annotation": input.gff3,
                "Basedir":    params.outdir,
                "Threads":    params.threads,
                "Analysis":   list(params.analyses),
            }
            if "BUSCO" in params.analyses:
                cfg["BUSCO_lineages"] = config["busco_lineages"]
            if "OMARK" in params.analyses:
                cfg["OMARK_db"]    = params.omark_db
                cfg["OMARK_taxid"] = params.taxa_id
            if "DETENGA" in params.analyses:
                cfg["DETENGA_db"] = params.detenga_db
            if "PROTHOMOLOGY" in params.analyses:
                cfg["PROTHOMOLOGY_tags"] = [
                    {k: v} for k, v in config.get("prothomology_dbs", {}).items()
                ]
            with open(output.yaml, "w") as fh:
                yaml.dump(cfg, fh, default_flow_style=False)

rule run_gaqet:
    input:
        yaml=f"{OUTDIR}/{{species}}/{{acc}}/AnnotationQC/gaqet/gaqet_config.yaml",
        fasta=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.fasta",
        gff3=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.gff3",
    output:
        stats=f"{OUTDIR}/{{species}}/{{acc}}/AnnotationQC/gaqet/{{acc}}_GAQET.stats.tsv",
    params:
        outbase=f"{OUTDIR}/{{species}}/{{acc}}/AnnotationQC/gaqet",
        taxa_id=lambda wildcards: ACC_TO_TAXID[wildcards.acc],
    benchmark:
        "benchmarks/gaqet/{species}/{acc}.tsv"
    log:
        "logs/gaqet/{species}/{acc}.log",
    shell:
        """
        if [ ! -s {input.fasta} ]; then
            echo "SKIP ({wildcards.acc}): download failed, skipping GAQET" >{log}
            touch {output.stats}
            exit 0
        fi
        GAQET \
            --yaml {input.yaml} \
            --species {wildcards.species} \
            --genome {input.fasta} \
            --annotation {input.gff3} \
            --taxid {params.taxa_id} \
            --outbase {params.outbase} \
            >{log} 2>&1
        """

# ── Compress genome files after all QC is done ────────────────────────────────
rule compress:
    input:
        fasta=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.fasta",
        gff3=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.gff3",
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
        gff3_gz=f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.gff3.gz",
    benchmark:
        "benchmarks/compress/{species}/{acc}.tsv"
    log:
        "logs/compress/{species}/{acc}.log",
    shell:
        """
        if [ ! -s {input.fasta} ]; then
            echo "SKIP ({wildcards.acc}): download failed, skipping compression" >{log}
            touch {output.fasta_gz} {output.gff3_gz}
            exit 0
        fi
        gzip -k {input.fasta} >{log} 2>&1
        gzip -k {input.gff3}  >>{log} 2>&1
        """

# ── Merge per-species GAQET stats into one TSV ────────────────────────────────
rule merge_gaqet_stats:
    input:
        expand(f"{OUTDIR}/{{species}}/{{acc}}/AnnotationQC/gaqet/{{acc}}_GAQET.stats.tsv",
               zip, species=SPECIES, acc=ACCESSIONS),
    output:
        merged=f"{OUTDIR}/report/all_species_GAQET.stats.tsv",
    run:
        import pandas as pd, os
        os.makedirs(os.path.dirname(output.merged), exist_ok=True)
        dfs = [pd.read_csv(f, sep="\t") for f in input if os.path.getsize(f) > 0]
        if dfs:
            pd.concat(dfs, ignore_index=True).to_csv(output.merged, sep="\t", index=False)
        else:
            open(output.merged, "w").close()

# ── GAQET_PLOT ─────────────────────────────────────────────────────────────────
rule run_gaqet_plot:
    input:
        merged=f"{OUTDIR}/report/all_species_GAQET.stats.tsv",
    output:
        plot=f"{OUTDIR}/report/all_species_GAQET.plot.{GAQET_PLOT_FMT}",
    log:
        "logs/gaqet_plot.log",
    shell:
        """
        if [ ! -s {input.merged} ]; then
            echo "SKIP: no GAQET stats available (all downloads may have failed)" >{log}
            touch {output.plot}
            exit 0
        fi
        GAQET_PLOT --input {input.merged} --output {output.plot} >{log} 2>&1
        """

# ── Custom HTML + TSV report ───────────────────────────────────────────────────
rule generate_report:
    input:
        fasta_gz=expand(
            f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb.fasta.gz",
            zip, species=SPECIES, acc=ACCESSIONS),
        quast=expand(
            f"{OUTDIR}/{{species}}/{{acc}}/AssemblyQC/quast/report.tsv",
            zip, species=SPECIES, acc=ACCESSIONS),
        busco=[
            f"{OUTDIR}/{s}/{a}/AssemblyQC/busco/{a}/short_summary.specific.{l}.{a}.txt"
            for s, a in SAMPLES for l in LINEAGES
        ],
        gaqet=expand(
            f"{OUTDIR}/{{species}}/{{acc}}/AnnotationQC/gaqet/{{acc}}_GAQET.stats.tsv",
            zip, species=SPECIES, acc=ACCESSIONS),
        gaqet_plot=f"{OUTDIR}/report/all_species_GAQET.plot.{GAQET_PLOT_FMT}",
    output:
        assembly_tsv=f"{OUTDIR}/report/annotseba_AssemblyQC.tsv",
        annotation_tsv=f"{OUTDIR}/report/annotseba_AnnotationQC.tsv",
        html=f"{OUTDIR}/report/annotseba_report.html",
    params:
        script=os.path.join(_basedir, "scripts", "generate_report.py"),
        outdir=OUTDIR,
        accessions=lambda _: config["accessions_file"],
        lineages=" ".join(LINEAGES),
        version=VERSION,
    log:
        "logs/generate_report.log",
    shell:
        """
        python {params.script} \
            --outdir {params.outdir} \
            --accessions {params.accessions} \
            --lineages {params.lineages} \
            --gaqet-plot {input.gaqet_plot} \
            --version {params.version} \
            >{log} 2>&1
        """

# ── Resource & carbon usage log ───────────────────────────────────────────────
rule compute_usage:
    input:
        report=f"{OUTDIR}/report/annotseba_report.html",
        benchmarks=BENCHMARK_FILES,
    output:
        log=f"{OUTDIR}/report/computer_usage.log",
    params:
        script=os.path.join(_basedir, "scripts", "compute_usage.py"),
        outdir=OUTDIR,
        cores=config.get("threads", 4),
        carbon=config.get("carbon_intensity", 475),
        tdp=config.get("cpu_tdp_per_core", 10),
        version=VERSION,
    log:
        "logs/compute_usage.log",
    shell:
        """
        python {params.script} \
            --outdir {params.outdir} \
            --benchmarks {input.benchmarks} \
            --cores {params.cores} \
            --carbon-intensity {params.carbon} \
            --tdp {params.tdp} \
            --version {params.version} \
            >{log} 2>&1
        """
