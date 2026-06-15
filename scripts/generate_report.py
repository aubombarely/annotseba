#!/usr/bin/env python3
"""
generate_report.py
Generate a self-contained HTML report and two TSV summaries
(AssemblyQC, AnnotationQC) from annotseba pipeline output.
"""

import argparse
import base64
import gzip
import io
import os
import re
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_accessions(path):
    """Return list of (species, acc) from accessions.txt."""
    samples = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) >= 2:
                samples.append((fields[0], fields[1]))
    return samples


def parse_quast(path):
    """Parse QUAST report.tsv -> dict of metric: value."""
    d = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] != "Assembly":
                d[parts[0]] = parts[1]
    return d


def parse_busco(path):
    """Parse BUSCO short_summary -> dict with Complete/Single/Duplicated/Fragmented/Missing."""
    with open(path) as fh:
        for line in fh:
            m = re.search(
                r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)",
                line,
            )
            if m:
                return {
                    "Complete":   float(m.group(1)),
                    "Single":     float(m.group(2)),
                    "Duplicated": float(m.group(3)),
                    "Fragmented": float(m.group(4)),
                    "Missing":    float(m.group(5)),
                    "Total":      int(m.group(6)),
                }
    return {}


def get_contig_lengths(fasta_path):
    """Return sorted-descending list of contig lengths from a (gzipped) FASTA."""
    lengths, cur = [], 0
    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    with opener(fasta_path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if cur:
                    lengths.append(cur)
                cur = 0
            else:
                cur += len(line)
    if cur:
        lengths.append(cur)
    return sorted(lengths, reverse=True)


def lineage_label(lin):
    """embryophyta_odb10 -> embryophyta"""
    return re.sub(r"_odb\d+$", "", lin)


# ---------------------------------------------------------------------------
# Plots  (all return base64 PNG strings)
# ---------------------------------------------------------------------------

PALETTE = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
    "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD",
]

# ── GAQET column grouping ─────────────────────────────────────────────────────
_GAQET_GROUP_ORDER = ["General", "AGAT", "BUSCO", "PSAURON", "DETENGA", "OMARK", "PROTHOMOLOGY"]

GAQET_GROUP_COLORS = {
    "General":      "#6c757d",
    "AGAT":         "#4C72B0",
    "BUSCO":        "#55A868",
    "PSAURON":      "#DD8452",
    "DETENGA":      "#C44E52",
    "OMARK":        "#8172B3",
    "PROTHOMOLOGY": "#937860",
}

# Rules are checked in order; first match wins. AGAT is the catch-all.
_GAQET_GROUP_RULES = [
    ("General",      lambda c: c in {"Species", "Accession", "NCBI_TaxID",
                                     "Assembly_Version", "Annotation_Version"}),
    ("BUSCO",        lambda c: "BUSCO" in c.upper()),
    ("PSAURON",      lambda c: "PSAURON" in c.upper()),
    ("DETENGA",      lambda c: "DETENGA" in c.upper()),
    ("OMARK",        lambda c: "OMARK" in c.upper()),
    ("PROTHOMOLOGY", lambda c: "PROTEIN" in c.upper()
                              or "TREMBL" in c.upper()
                              or "SWISSPROT" in c.upper()),
    ("AGAT",         lambda c: True),   # everything else is an AGAT gene-model stat
]


def group_gaqet_columns(columns):
    """Return [(group_name, [col, ...]), ...] preserving _GAQET_GROUP_ORDER."""
    buckets = {g: [] for g in _GAQET_GROUP_ORDER}
    for col in columns:
        for grp, test in _GAQET_GROUP_RULES:
            if test(col):
                buckets[grp].append(col)
                break
    return [(g, buckets[g]) for g in _GAQET_GROUP_ORDER if buckets[g]]


def _fig_to_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=130)
    buf.seek(0)
    data = base64.b64encode(buf.read()).decode()
    plt.close(fig)
    return data


def plot_assembly_sizes(df):
    labels = [s.replace("_", " ") for s in df["Species"]]
    sizes  = pd.to_numeric(df["Total_length_bp"], errors="coerce") / 1e6
    n = len(labels)
    fig, ax = plt.subplots(figsize=(max(5, n * 0.9 + 1), 4))
    bars = ax.bar(range(n), sizes, color=PALETTE[:n], edgecolor="white")
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=9)
    ax.set_ylabel("Assembly size (Mb)")
    ax.set_title("Total assembly size per species")
    ax.bar_label(bars, fmt="%.0f Mb", fontsize=8, padding=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    return _fig_to_b64(fig)


def plot_cumulative(contig_data):
    """contig_data: list of (label, lengths_sorted_desc)"""
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, (label, lengths) in enumerate(contig_data):
        if not lengths:
            continue
        cumsum = pd.Series(lengths).cumsum() / 1e6
        ax.plot(range(1, len(cumsum) + 1), cumsum,
                color=PALETTE[i % len(PALETTE)],
                label=label.replace("_", " "), linewidth=1.5)
    ax.set_xlabel("Contig index (sorted by length, descending)")
    ax.set_ylabel("Cumulative length (Mb)")
    ax.set_title("Cumulative contig length")
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    return _fig_to_b64(fig)


def embed_image(path):
    """Embed any image file as a data-URI src string."""
    ext  = Path(path).suffix.lstrip(".").lower()
    mime = {"png": "image/png", "jpg": "image/jpeg",
            "jpeg": "image/jpeg", "svg": "image/svg+xml"}.get(ext, "image/png")
    with open(path, "rb") as fh:
        data = base64.b64encode(fh.read()).decode()
    return f"data:{mime};base64,{data}"


# ---------------------------------------------------------------------------
# HTML helpers
# ---------------------------------------------------------------------------

_CSS = """
<style>
body{font-family:Arial,sans-serif;margin:30px;color:#222;}
h1{color:#2c5f8a;}
h2{color:#3a7ab5;border-bottom:2px solid #3a7ab5;padding-bottom:4px;margin-top:40px;}
h3{color:#555;margin-top:20px;}
p.meta{color:#666;font-size:.9em;}
.plots{display:flex;flex-wrap:wrap;gap:24px;margin-top:16px;}
.plot-box{text-align:center;}
.plot-box img{max-width:700px;width:100%;border:1px solid #ddd;border-radius:4px;}
.tbl-wrap{overflow-x:auto;margin-top:12px;}
table{border-collapse:collapse;font-size:.82em;min-width:600px;}
th,td{border:1px solid #ccc;padding:5px 9px;text-align:right;white-space:nowrap;}
th{background:#3a7ab5;color:#fff;cursor:pointer;user-select:none;}
th:first-child,td:first-child,th:nth-child(2),td:nth-child(2){text-align:left;}
tr:nth-child(even){background:#f4f8fc;}
tr:hover{background:#dce8f4;}
th[data-sorted=asc]::after{content:' ▲';font-size:.75em;}
th[data-sorted=desc]::after{content:' ▼';font-size:.75em;}
tr.grp-hdr th{cursor:default;text-align:center;font-size:.78em;
  letter-spacing:.06em;text-transform:uppercase;border-bottom:2px solid #fff;}
tr.grp-hdr th::after{content:'' !important;}
</style>
"""

_JS = """
<script>
function sortTable(th){
  var tb=th.closest('table').tBodies[0];
  var rows=Array.from(tb.rows);
  var col=Array.from(th.parentNode.children).indexOf(th);
  var asc=th.dataset.asc!=='true';
  th.dataset.asc=asc;
  th.parentNode.querySelectorAll('th').forEach(function(t){delete t.dataset.sorted;});
  th.dataset.sorted=asc?'asc':'desc';
  rows.sort(function(a,b){
    var av=a.cells[col].textContent.trim();
    var bv=b.cells[col].textContent.trim();
    var an=parseFloat(av),bn=parseFloat(bv);
    if(!isNaN(an)&&!isNaN(bn)) return asc?an-bn:bn-an;
    return asc?av.localeCompare(bv):bv.localeCompare(av);
  });
  rows.forEach(function(r){tb.appendChild(r);});
}
</script>
"""


def _df_to_html_table(df, table_id):
    thead = "".join(
        f'<th onclick="sortTable(this)">{c}</th>' for c in df.columns
    )
    rows = "".join(
        "<tr>" + "".join(f"<td>{v}</td>" for v in row) + "</tr>"
        for row in df.itertuples(index=False, name=None)
    )
    return (f'<div class="tbl-wrap">'
            f'<table id="{table_id}">'
            f"<thead><tr>{thead}</tr></thead>"
            f"<tbody>{rows}</tbody>"
            f"</table></div>")


def _df_to_grouped_html_table(df, table_id, col_groups):
    """Two-row-header sortable table with colour-coded group spans."""
    # Header row 1 — group name cells with colspan
    hrow1 = ""
    for grp_name, cols in col_groups:
        color = GAQET_GROUP_COLORS.get(grp_name, "#3a7ab5")
        hrow1 += (f'<th colspan="{len(cols)}" '
                  f'style="background:{color}">'
                  f'{grp_name}</th>')

    # Header row 2 — individual sortable column headers
    all_cols = [col for _, cols in col_groups for col in cols]
    hrow2 = "".join(
        f'<th onclick="sortTable(this)">{c}</th>' for c in all_cols
    )

    # Body — reorder df to match grouped column order
    existing = [c for c in all_cols if c in df.columns]
    rows = "".join(
        "<tr>" + "".join(f"<td>{v}</td>" for v in row) + "</tr>"
        for row in df.reindex(columns=existing).itertuples(index=False, name=None)
    )

    return (f'<div class="tbl-wrap">'
            f'<table id="{table_id}">'
            f'<thead>'
            f'<tr class="grp-hdr">{hrow1}</tr>'
            f'<tr>{hrow2}</tr>'
            f'</thead>'
            f'<tbody>{rows}</tbody>'
            f'</table></div>')


def _png_img(b64, alt):
    if not b64:
        return ""
    return (f'<div class="plot-box">'
            f'<img src="data:image/png;base64,{b64}" alt="{alt}">'
            f'</div>')


def build_html(asm_df, ann_df, asm_size_b64, cumlen_b64,
               gaqet_img_src, version, date_str):
    ann_col_groups = group_gaqet_columns(list(ann_df.columns)) if not ann_df.empty else []
    ann_table = (_df_to_grouped_html_table(ann_df, "tbl_ann", ann_col_groups)
                 if ann_col_groups else "<p>No annotation data available.</p>")

    gaqet_html = ""
    if gaqet_img_src:
        gaqet_html = (f'<div class="plot-box">'
                      f'<h3>GAQET annotation quality</h3>'
                      f'<img src="{gaqet_img_src}" alt="GAQET plot">'
                      f'</div>')

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>annotseba QC Report</title>
{_CSS}
{_JS}
</head>
<body>
<h1>annotseba QC Report</h1>
<p class="meta">Version: {version}&nbsp;|&nbsp;Generated: {date_str}</p>

<h2>Assembly QC</h2>
{_df_to_html_table(asm_df, "tbl_asm")}
<div class="plots">
  {_png_img(asm_size_b64, "Assembly sizes")}
  {_png_img(cumlen_b64,   "Cumulative contig length")}
</div>

<h2>Annotation QC</h2>
{ann_table}
<div class="plots">
  {gaqet_html}
</div>

</body>
</html>"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Generate annotseba HTML+TSV reports")
    ap.add_argument("--outdir",     required=True,  help="Pipeline OUTDIR")
    ap.add_argument("--accessions", required=True,  help="Path to accessions.txt")
    ap.add_argument("--lineages",   required=True,  nargs="+", help="BUSCO lineages")
    ap.add_argument("--gaqet-plot", required=True,  help="Path to GAQET_PLOT output image")
    ap.add_argument("--version",    default="unknown")
    args = ap.parse_args()

    outdir   = Path(args.outdir)
    lineages = args.lineages
    samples  = parse_accessions(args.accessions)

    report_dir = outdir / "report"
    report_dir.mkdir(parents=True, exist_ok=True)

    # ── AssemblyQC data ────────────────────────────────────────────────────
    print("Collecting AssemblyQC data ...", flush=True)
    asm_rows    = []
    contig_data = []

    for species, acc in samples:
        base = outdir / species / acc
        row  = {"Species": species, "Accession": acc}

        qpath = base / "AssemblyQC" / "quast" / "report.tsv"
        if qpath.exists() and qpath.stat().st_size > 0:
            q = parse_quast(qpath)
            row["Total_length_bp"]   = q.get("Total length",   "NA")
            row["Num_contigs"]       = q.get("# contigs",       "NA")
            row["Largest_contig_bp"] = q.get("Largest contig",  "NA")
            row["N50_bp"]            = q.get("N50",             "NA")
            row["N75_bp"]            = q.get("N75",             "NA")
            row["N90_bp"]            = q.get("N90",             "NA")
            row["L50"]               = q.get("L50",             "NA")
            row["L75"]               = q.get("L75",             "NA")
            row["L90"]               = q.get("L90",             "NA")
            row["GC_pct"]            = q.get("GC (%)",          "NA")
        else:
            print(f"  Warning: QUAST report missing for {species}/{acc}", file=sys.stderr)

        fasta_gz = base / "genome" / f"{species}_{acc}_asb.fasta.gz"
        if fasta_gz.exists() and fasta_gz.stat().st_size > 0:
            try:
                lengths = get_contig_lengths(fasta_gz)
                contig_data.append((species, lengths))
            except Exception as e:
                print(f"  Warning: contig lengths failed for {species}/{acc}: {e}",
                      file=sys.stderr)

        for lin in lineages:
            lbl   = lineage_label(lin)
            bpath = (base / "AssemblyQC" / "busco" / lin / acc /
                     f"short_summary.specific.{lin}.{acc}.txt")
            if bpath.exists() and bpath.stat().st_size > 0:
                b = parse_busco(bpath)
                row[f"BUSCO_{lbl}_C%"]  = b.get("Complete",   "NA")
                row[f"BUSCO_{lbl}_S%"]  = b.get("Single",     "NA")
                row[f"BUSCO_{lbl}_D%"]  = b.get("Duplicated", "NA")
                row[f"BUSCO_{lbl}_F%"]  = b.get("Fragmented", "NA")
                row[f"BUSCO_{lbl}_M%"]  = b.get("Missing",    "NA")
                row[f"BUSCO_{lbl}_n"]   = b.get("Total",      "NA")
            else:
                for s in ["C%", "S%", "D%", "F%", "M%", "n"]:
                    row[f"BUSCO_{lbl}_{s}"] = "NA"

        asm_rows.append(row)

    asm_df = pd.DataFrame(asm_rows)

    # ── AnnotationQC data ──────────────────────────────────────────────────
    print("Collecting AnnotationQC data ...", flush=True)
    ann_dfs = []
    for species, acc in samples:
        gpath = (outdir / species / acc /
                 "AnnotationQC" / "gaqet" / f"{species}_{acc}_GAQET.stats.tsv")
        if gpath.exists() and gpath.stat().st_size > 0:
            try:
                df = pd.read_csv(gpath, sep="\t")
                if "Accession" not in df.columns:
                    df.insert(0, "Accession", acc)
                if "Species" not in df.columns:
                    df.insert(0, "Species", species)
                ann_dfs.append(df)
            except Exception as e:
                print(f"  Warning: could not parse {gpath}: {e}", file=sys.stderr)
        else:
            print(f"  Warning: GAQET stats missing or empty for {species}/{acc}: {gpath}",
                  file=sys.stderr)

    ann_df = (pd.concat(ann_dfs, ignore_index=True)
              if ann_dfs else pd.DataFrame(columns=["Species", "Accession"]))

    # ── Write TSVs ─────────────────────────────────────────────────────────
    asm_tsv = report_dir / "annotseba_AssemblyQC.tsv"
    ann_tsv = report_dir / "annotseba_AnnotationQC.tsv"
    asm_df.to_csv(asm_tsv, sep="\t", index=False)
    ann_df.to_csv(ann_tsv, sep="\t", index=False)
    print(f"AssemblyQC TSV   -> {asm_tsv}")
    print(f"AnnotationQC TSV -> {ann_tsv}")

    # ── Plots ──────────────────────────────────────────────────────────────
    print("Generating plots ...", flush=True)
    asm_size_b64 = cumlen_b64 = None

    if not asm_df.empty and "Total_length_bp" in asm_df.columns:
        try:
            asm_size_b64 = plot_assembly_sizes(asm_df)
        except Exception as e:
            print(f"  Warning: assembly size plot failed: {e}", file=sys.stderr)

    if contig_data:
        try:
            cumlen_b64 = plot_cumulative(contig_data)
        except Exception as e:
            print(f"  Warning: cumulative length plot failed: {e}", file=sys.stderr)

    gaqet_img_src = None
    if os.path.exists(args.gaqet_plot) and os.path.getsize(args.gaqet_plot) > 0:
        try:
            gaqet_img_src = embed_image(args.gaqet_plot)
        except Exception as e:
            print(f"  Warning: could not embed GAQET plot: {e}", file=sys.stderr)

    # ── Write HTML ─────────────────────────────────────────────────────────
    html_path = report_dir / "annotseba_report.html"
    date_str  = datetime.now().strftime("%Y-%m-%d %H:%M")
    html = build_html(asm_df, ann_df, asm_size_b64, cumlen_b64,
                      gaqet_img_src, args.version, date_str)
    html_path.write_text(html)
    print(f"HTML report      -> {html_path}")


if __name__ == "__main__":
    main()
