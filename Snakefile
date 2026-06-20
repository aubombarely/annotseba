import json
import os

# workflow.basedir is the directory containing this Snakefile — always
# absolute regardless of where the user runs the pipeline from.
_basedir = workflow.basedir

with open(os.path.join(_basedir, "VERSION")) as _vf:
    VERSION = _vf.read().strip()

print(f"annotseba v{VERSION}")

configfile: os.path.join(_basedir, "config", "config.yaml")

# ── Layout / feature flags ────────────────────────────────────────────────────
# byspecies_root: path to a Phytolaeno BySpecies/ root (enables v2 layout).
# When set, outputs land inside the pre-existing species tree created by
# generate_species_dir.py instead of a flat annotseba_run/ directory.
BYSPECIES_ROOT = config.get("byspecies_root", "").rstrip("/")
RUN_GAQET      = config.get("run_gaqet", True)
USE_V2         = bool(BYSPECIES_ROOT)

_ks = config.get("keep_source", False)
KEEP_SOURCE    = _ks if isinstance(_ks, bool) else str(_ks).lower() in ("true", "1", "yes")
GAQET_PLOT_FMT = config.get("gaqet_plot_format", "png")

# ── Load accession list ───────────────────────────────────────────────────────
# Format: species<TAB>accession[<TAB>taxa_id[<TAB>prefix]]
# In v2 mode taxa_id may be omitted — it is read from 00_Taxonomy/taxonomy.json.

def _taxid_from_taxonomy(species):
    """Read taxid from Phytolaeno 00_Taxonomy/taxonomy.json (v2 only)."""
    if not BYSPECIES_ROOT:
        return ""
    path = os.path.join(BYSPECIES_ROOT, species, "00_Taxonomy", "taxonomy.json")
    try:
        with open(path) as fh:
            return str(json.load(fh).get("taxid", ""))
    except (FileNotFoundError, KeyError, ValueError):
        return ""

SAMPLES        = []
ACC_TO_SPECIES = {}
ACC_TO_TAXID   = {}
ACC_TO_PREFIX  = {}
with open(config["accessions_file"]) as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        species, acc = fields[0], fields[1]
        # taxa_id: explicit in file > taxonomy.json (v2) > empty
        if len(fields) > 2 and fields[2].strip() and fields[2].strip().upper() != "NA":
            taxa_id = fields[2].strip()
        elif USE_V2:
            taxa_id = _taxid_from_taxonomy(species)
        else:
            taxa_id = ""
        prefix = (fields[3] if len(fields) > 3 and fields[3].strip().upper() != "NA"
                  else config.get("rename_prefix", "seq"))
        SAMPLES.append((species, acc))
        ACC_TO_SPECIES[acc] = species
        ACC_TO_TAXID[acc]   = taxa_id
        ACC_TO_PREFIX[acc]  = prefix

SPECIES    = [s for s, a in SAMPLES]
ACCESSIONS = [a for s, a in SAMPLES]
LINEAGES   = config["busco_lineages"]

# ── Directory layout ──────────────────────────────────────────────────────────
#
# v2 (USE_V2 = True):
#   {byspecies_root}/{species}/01_Genomes/{acc}/
#       00_SourceFiles/          — genome FASTA, equiv_seqID.txt, source files
#       01_Assembly_QC/          — QUAST, BUSCO, assembly-stats
#       02_WorkingFiles/         — logs, benchmarks, report (annotseba runtime)
#       03_Annotations/
#           01_Gene_model_annotations/
#               00_NCBI_RefSeq/  — renamed GFF3, GAQET output
#
# v1 (USE_V2 = False):
#   {RUNDIR}/results/{species}/{acc}/  (original flat layout)

if USE_V2:
    _V2B = f"{BYSPECIES_ROOT}/{{species}}/01_Genomes/{{acc}}"
    _SRC     = f"{_V2B}/00_SourceFiles"
    _AQC     = f"{_V2B}/01_Assembly_QC"
    _WK      = f"{_V2B}/02_WorkingFiles"
    _ANN     = f"{_V2B}/03_Annotations/01_Gene_model_annotations/00_NCBI_RefSeq"
    _LDIR    = f"{_WK}/logs"
    _BDIR    = f"{_WK}/benchmarks"

    # Report dir: per-assembly for single run, configurable for batch
    if len(SAMPLES) == 1:
        _RDIR = (f"{BYSPECIES_ROOT}/{SAMPLES[0][0]}/01_Genomes/"
                 f"{SAMPLES[0][1]}/02_WorkingFiles/report")
    else:
        _RDIR = config.get("report_dir", f"{BYSPECIES_ROOT}/_annotseba_reports")

    # File path templates (wildcards {species} and {acc} as literal strings)
    _SRC_FNA   = f"{_SRC}/{{acc}}.fna"
    _SRC_GFF   = f"{_SRC}/{{acc}}.gff3"
    _ASB_FASTA = f"{_SRC}/{{species}}_{{acc}}_asb.fasta"
    _EQUIV     = f"{_SRC}/{{species}}_{{acc}}_asb.equiv_seqID.txt"
    _ASB_GFF3  = f"{_ANN}/{{species}}_{{acc}}_asb.gff3"
    _FASTA_GZ  = f"{_SRC}/{{species}}_{{acc}}_asb.fasta.gz"
    _GFF3_GZ   = f"{_ANN}/{{species}}_{{acc}}_asb.gff3.gz"
    _QUAST_TSV = f"{_AQC}/quast/report.tsv"
    _ASTATS    = f"{_AQC}/assembly_stats/stats.txt"
    _BUSCO     = f"{_AQC}/busco/{{lineage}}/{{acc}}/short_summary.specific.{{lineage}}.{{acc}}.txt"
    _GAQET_YML = f"{_ANN}/gaqet_config.yaml"
    _GAQET_TSV = f"{_ANN}/{{species}}_{{acc}}_GAQET.stats.tsv"
    _MGAQET    = f"{_RDIR}/all_species_GAQET.stats.tsv"
    _GPLOT     = f"{_RDIR}/all_species_GAQET.plot.{GAQET_PLOT_FMT}"

    def _log(rule):   return f"{_LDIR}/{rule}.log"
    def _bench(rule): return f"{_BDIR}/{rule}.tsv"
    _BUSCO_LOG   = f"{_LDIR}/busco_{{lineage}}.log"
    _BUSCO_BENCH = f"{_BDIR}/busco_{{lineage}}.tsv"

    # Concrete benchmark paths for compute_usage
    _BFILES = (
        [f"{BYSPECIES_ROOT}/{s}/01_Genomes/{a}/02_WorkingFiles/benchmarks/download.tsv"
         for s, a in SAMPLES] +
        [f"{BYSPECIES_ROOT}/{s}/01_Genomes/{a}/02_WorkingFiles/benchmarks/rename.tsv"
         for s, a in SAMPLES] +
        [f"{BYSPECIES_ROOT}/{s}/01_Genomes/{a}/02_WorkingFiles/benchmarks/rename_gff3.tsv"
         for s, a in SAMPLES] +
        [f"{BYSPECIES_ROOT}/{s}/01_Genomes/{a}/02_WorkingFiles/benchmarks/quast.tsv"
         for s, a in SAMPLES] +
        [f"{BYSPECIES_ROOT}/{s}/01_Genomes/{a}/02_WorkingFiles/benchmarks/assembly_stats.tsv"
         for s, a in SAMPLES] +
        [f"{BYSPECIES_ROOT}/{s}/01_Genomes/{a}/02_WorkingFiles/benchmarks/busco_{l}.tsv"
         for s, a in SAMPLES for l in LINEAGES] +
        ([f"{BYSPECIES_ROOT}/{s}/01_Genomes/{a}/02_WorkingFiles/benchmarks/gaqet.tsv"
          for s, a in SAMPLES] if RUN_GAQET else []) +
        [f"{BYSPECIES_ROOT}/{s}/01_Genomes/{a}/02_WorkingFiles/benchmarks/compress.tsv"
         for s, a in SAMPLES]
    )

else:  # v1 legacy layout
    RUNDIR   = config.get("rundir",   "annotseba_run")
    OUTDIR   = config.get("outdir",   f"{RUNDIR}/results")
    _LOGDIR  = config.get("logdir",   f"{RUNDIR}/logs")
    _BENCHD  = config.get("benchdir", f"{RUNDIR}/benchmarks")
    _RDIR    = f"{OUTDIR}/report"

    _BASE      = f"{OUTDIR}/{{species}}/{{acc}}"
    _SRC_FNA   = f"{_BASE}/genome/{{acc}}.fna"
    _SRC_GFF   = f"{_BASE}/genome/{{acc}}.gff3"
    _ASB_FASTA = f"{_BASE}/genome/{{species}}_{{acc}}_asb.fasta"
    _EQUIV     = f"{_BASE}/genome/{{species}}_{{acc}}_asb.equiv_seqID.txt"
    _ASB_GFF3  = f"{_BASE}/genome/{{species}}_{{acc}}_asb.gff3"
    _FASTA_GZ  = f"{_BASE}/genome/{{species}}_{{acc}}_asb.fasta.gz"
    _GFF3_GZ   = f"{_BASE}/genome/{{species}}_{{acc}}_asb.gff3.gz"
    _QUAST_TSV = f"{_BASE}/AssemblyQC/quast/report.tsv"
    _ASTATS    = f"{_BASE}/AssemblyQC/assembly_stats/stats.txt"
    _BUSCO     = f"{_BASE}/AssemblyQC/busco/{{lineage}}/{{acc}}/short_summary.specific.{{lineage}}.{{acc}}.txt"
    _GAQET_YML = f"{_BASE}/AnnotationQC/gaqet/gaqet_config.yaml"
    _GAQET_TSV = f"{_BASE}/AnnotationQC/gaqet/{{species}}_{{acc}}_GAQET.stats.tsv"
    _MGAQET    = f"{_RDIR}/all_species_GAQET.stats.tsv"
    _GPLOT     = f"{_RDIR}/all_species_GAQET.plot.{GAQET_PLOT_FMT}"

    def _log(rule):   return f"{_LOGDIR}/{rule}/{{species}}/{{acc}}.log"
    def _bench(rule): return f"{_BENCHD}/{rule}/{{species}}/{{acc}}.tsv"
    _BUSCO_LOG   = f"{_LOGDIR}/busco/{{species}}/{{acc}}/{{lineage}}.log"
    _BUSCO_BENCH = f"{_BENCHD}/busco/{{species}}/{{acc}}/{{lineage}}.tsv"

    _BFILES = (
        expand(f"{_BENCHD}/download/{{species}}/{{acc}}.tsv",       zip, species=SPECIES, acc=ACCESSIONS) +
        expand(f"{_BENCHD}/rename/{{species}}/{{acc}}.tsv",         zip, species=SPECIES, acc=ACCESSIONS) +
        expand(f"{_BENCHD}/rename_gff3/{{species}}/{{acc}}.tsv",    zip, species=SPECIES, acc=ACCESSIONS) +
        expand(f"{_BENCHD}/quast/{{species}}/{{acc}}.tsv",          zip, species=SPECIES, acc=ACCESSIONS) +
        expand(f"{_BENCHD}/assembly_stats/{{species}}/{{acc}}.tsv", zip, species=SPECIES, acc=ACCESSIONS) +
        [f"{_BENCHD}/busco/{s}/{a}/{l}.tsv" for s, a in SAMPLES for l in LINEAGES] +
        (expand(f"{_BENCHD}/gaqet/{{species}}/{{acc}}.tsv",         zip, species=SPECIES, acc=ACCESSIONS)
         if RUN_GAQET else []) +
        expand(f"{_BENCHD}/compress/{{species}}/{{acc}}.tsv",       zip, species=SPECIES, acc=ACCESSIONS)
    )

# ── Target rule ───────────────────────────────────────────────────────────────
rule all:
    input:
        expand(_ASB_FASTA,  zip, species=SPECIES, acc=ACCESSIONS),
        expand(_EQUIV,      zip, species=SPECIES, acc=ACCESSIONS),
        expand(_QUAST_TSV,  zip, species=SPECIES, acc=ACCESSIONS),
        expand(_ASTATS,     zip, species=SPECIES, acc=ACCESSIONS),
        [_BUSCO.format(species=s, acc=a, lineage=l)
         for s, a in SAMPLES for l in LINEAGES],
        *(expand(_GAQET_TSV, zip, species=SPECIES, acc=ACCESSIONS) if RUN_GAQET else []),
        expand(_FASTA_GZ,   zip, species=SPECIES, acc=ACCESSIONS),
        *(expand(_GFF3_GZ,  zip, species=SPECIES, acc=ACCESSIONS) if RUN_GAQET else []),
        f"{_RDIR}/annotseba_AssemblyQC.tsv",
        f"{_RDIR}/annotseba_AnnotationQC.tsv",
        f"{_RDIR}/annotseba_report.html",
        f"{_RDIR}/computer_usage.log",

# ── Download genome from NCBI ─────────────────────────────────────────────────
rule download_genome:
    output:
        fasta=(_SRC_FNA if KEEP_SOURCE else temp(_SRC_FNA)),
        gff3= (_SRC_GFF if KEEP_SOURCE else temp(_SRC_GFF)),
    params:
        workdir=_SRC if USE_V2 else f"{OUTDIR}/{{species}}/{{acc}}/genome",
        zipfile=(f"{_SRC}/ncbi_dataset.zip" if USE_V2
                 else f"{OUTDIR}/{{species}}/{{acc}}/genome/ncbi_dataset.zip"),
    benchmark: _bench("download")
    log:        _log("download")
    retries: 3
    shell:
        """
        mkdir -p {params.workdir}

        datasets download genome accession {wildcards.acc} \
            --include genome,gff3 \
            --filename {params.zipfile} \
            2>{log}

        unzip -o {params.zipfile} -d {params.workdir}/tmp >>{log} 2>&1

        fasta=$(find {params.workdir}/tmp -name "*.fna" | head -1)
        if [ -z "$fasta" ]; then
            echo "WARNING: no .fna file found for {wildcards.acc} — creating empty placeholder" >>{log}
            touch {output.fasta} {output.gff3}
            rm -rf {params.workdir}/tmp {params.zipfile}
            exit 0
        fi
        mv "$fasta" {output.fasta}

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
        fasta=_SRC_FNA,
    output:
        fasta=_ASB_FASTA,
        equiv=_EQUIV,
    params:
        prefix=lambda wildcards: ACC_TO_PREFIX[wildcards.acc],
        out_basename=(_SRC + "/{species}_{acc}_asb" if USE_V2
                      else f"{OUTDIR}/{{species}}/{{acc}}/genome/{{species}}_{{acc}}_asb"),
        script=config.get("ncbi_fasta_rename_script", "NCBI_FastaRename"),
    benchmark: _bench("rename")
    log:        _log("rename")
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
        gff3=_SRC_GFF,
        equiv=_EQUIV,
    output:
        gff3=_ASB_GFF3,
    benchmark: _bench("rename_gff3")
    log:        _log("rename_gff3")
    shell:
        """
        if [ ! -s {input.gff3} ]; then
            echo "SKIP ({wildcards.acc}): download failed, skipping GFF3 rename" >{log}
            touch {output.gff3}
            exit 0
        fi
        _gff_in=$(realpath {input.gff3})
        _tsv_in=$(realpath {input.equiv})
        _gff_out=$(realpath -m {output.gff3})
        _log=$(realpath -m {log})
        mkdir -p $(dirname $_gff_out)
        cd $(dirname $_gff_in)
        agat_sq_rename_seqid.pl \
            --gff "$_gff_in" \
            --tsv "$_tsv_in" \
            --output "$_gff_out" \
            >"$_log" 2>&1
        """

# ── QUAST assembly statistics ─────────────────────────────────────────────────
rule run_quast:
    input:
        fasta=_ASB_FASTA,
        gff3=_ASB_GFF3,
    output:
        report=_QUAST_TSV,
    params:
        outdir=(_AQC + "/quast" if USE_V2
                else f"{OUTDIR}/{{species}}/{{acc}}/AssemblyQC/quast"),
        threads=config.get("threads", 4),
        min_contig=config.get("quast_min_contig", 500),
    benchmark: _bench("quast")
    log:        _log("quast")
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
        fasta=_ASB_FASTA,
    output:
        stats=_ASTATS,
    benchmark: _bench("assembly_stats")
    log:        _log("assembly_stats")
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
        fasta=_ASB_FASTA,
    output:
        summary=_BUSCO,
    params:
        outdir=(_AQC + "/busco/{lineage}" if USE_V2
                else f"{OUTDIR}/{{species}}/{{acc}}/AssemblyQC/busco/{{lineage}}"),
        threads=config.get("threads", 4),
        mode="genome",
        downloads_path=config.get("busco_downloads_path", "busco_downloads"),
    benchmark: _BUSCO_BENCH
    log:        _BUSCO_LOG
    shell:
        """
        if [ ! -s {input.fasta} ]; then
            echo "SKIP ({wildcards.acc}): download failed, skipping BUSCO" >{log}
            mkdir -p $(dirname {output.summary})
            touch {output.summary}
            exit 0
        fi
        rm -rf {params.outdir}/{wildcards.acc}
        busco \
            -i {input.fasta} \
            -o {wildcards.acc} \
            --out_path {params.outdir} \
            -l {wildcards.lineage} \
            -m {params.mode} \
            -c {params.threads} \
            --download_path {params.downloads_path} \
            >{log} 2>&1
        """

# ── GAQET2 annotation quality (optional: --run_gaqet / run_gaqet: true) ───────
if RUN_GAQET:
    rule write_gaqet_yaml:
        input:
            fasta=_ASB_FASTA,
            gff3=_ASB_GFF3,
        output:
            yaml=_GAQET_YML,
        params:
            outdir=(_ANN if USE_V2
                    else f"{OUTDIR}/{{species}}/{{acc}}/AnnotationQC/gaqet"),
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
            yaml=_GAQET_YML,
            fasta=_ASB_FASTA,
            gff3=_ASB_GFF3,
        output:
            stats=_GAQET_TSV,
        params:
            outbase=(_ANN if USE_V2
                     else f"{OUTDIR}/{{species}}/{{acc}}/AnnotationQC/gaqet"),
            taxa_id=lambda wildcards: ACC_TO_TAXID[wildcards.acc],
        benchmark: _bench("gaqet")
        log:        _log("gaqet")
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
            mv {params.outbase}/{wildcards.species}_GAQET.stats.tsv {output.stats}
            """

# ── Compress genome files after all QC is done ────────────────────────────────
rule compress:
    input:
        fasta=_ASB_FASTA,
        gff3=_ASB_GFF3,
        quast=_QUAST_TSV,
        assembly_stats=_ASTATS,
        busco=lambda wildcards: [
            _BUSCO.format(species=wildcards.species, acc=wildcards.acc, lineage=l)
            for l in LINEAGES
        ],
        gaqet=lambda wildcards: (
            [_GAQET_TSV.format(species=wildcards.species, acc=wildcards.acc)]
            if RUN_GAQET else []
        ),
    output:
        fasta_gz=_FASTA_GZ,
        gff3_gz= _GFF3_GZ,
    benchmark: _bench("compress")
    log:        _log("compress")
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
        expand(_GAQET_TSV, zip, species=SPECIES, acc=ACCESSIONS) if RUN_GAQET else [],
    output:
        merged=_MGAQET,
    run:
        import pandas as pd, os
        os.makedirs(os.path.dirname(output.merged), exist_ok=True)
        if not RUN_GAQET or not input:
            open(output.merged, "w").close()
        else:
            dfs = [pd.read_csv(f, sep="\t") for f in input if os.path.getsize(f) > 0]
            if dfs:
                pd.concat(dfs, ignore_index=True).to_csv(output.merged, sep="\t", index=False)
            else:
                open(output.merged, "w").close()

# ── GAQET plot ────────────────────────────────────────────────────────────────
rule run_gaqet_plot:
    input:
        merged=_MGAQET,
    output:
        plot=_GPLOT,
    log:
        f"{_RDIR}/gaqet_plot.log"
    shell:
        """
        if [ ! -s {input.merged} ]; then
            echo "SKIP: no GAQET stats available" >{log}
            touch {output.plot}
            exit 0
        fi
        GAQET_PLOT --input {input.merged} --output {output.plot} >{log} 2>&1
        """

# ── Custom HTML + TSV report ───────────────────────────────────────────────────
rule generate_report:
    input:
        fasta_gz=expand(_FASTA_GZ, zip, species=SPECIES, acc=ACCESSIONS),
        quast=expand(_QUAST_TSV,   zip, species=SPECIES, acc=ACCESSIONS),
        busco=[_BUSCO.format(species=s, acc=a, lineage=l)
               for s, a in SAMPLES for l in LINEAGES],
        gaqet=(expand(_GAQET_TSV, zip, species=SPECIES, acc=ACCESSIONS)
               if RUN_GAQET else [_MGAQET]),
        gaqet_plot=_GPLOT,
    output:
        assembly_tsv=f"{_RDIR}/annotseba_AssemblyQC.tsv",
        annotation_tsv=f"{_RDIR}/annotseba_AnnotationQC.tsv",
        html=f"{_RDIR}/annotseba_report.html",
    params:
        script=os.path.join(_basedir, "scripts", "generate_report.py"),
        report_dir=_RDIR,
        accessions=lambda _: config["accessions_file"],
        lineages=" ".join(LINEAGES),
        version=VERSION,
        layout_flag=(f"--byspecies_root {BYSPECIES_ROOT}" if USE_V2
                     else f"--outdir {OUTDIR}"),
    log:
        f"{_RDIR}/generate_report.log"
    shell:
        """
        python {params.script} \
            {params.layout_flag} \
            --report_dir {params.report_dir} \
            --accessions {params.accessions} \
            --lineages {params.lineages} \
            --gaqet-plot {input.gaqet_plot} \
            --version {params.version} \
            >{log} 2>&1
        """

# ── Resource & carbon usage log ───────────────────────────────────────────────
rule compute_usage:
    input:
        report=f"{_RDIR}/annotseba_report.html",
        benchmarks=_BFILES,
    output:
        log=f"{_RDIR}/computer_usage.log",
    params:
        script=os.path.join(_basedir, "scripts", "compute_usage.py"),
        report_dir=_RDIR,
        cores=config.get("threads", 4),
        carbon=config.get("carbon_intensity", 475),
        tdp=config.get("cpu_tdp_per_core", 10),
        version=VERSION,
    log:
        f"{_RDIR}/compute_usage.log"
    shell:
        """
        python {params.script} \
            --outdir {params.report_dir} \
            --benchmarks {input.benchmarks} \
            --cores {params.cores} \
            --carbon-intensity {params.carbon} \
            --tdp {params.tdp} \
            --version {params.version} \
            >{log} 2>&1
        """
