#!/usr/bin/env python3
"""
compute_usage.py
Collect Snakemake benchmark data + disk usage and write computer_usage.log
with a carbon footprint estimate.

Carbon model:
  energy (kWh) = total_cpu_seconds / 3600 * tdp_per_core_W / 1000
  CO2 (g)      = energy * carbon_intensity (gCO2/kWh)

Default carbon intensity: 475 gCO2/kWh (IEA 2022 world average).
Default CPU TDP: 10 W/core (conservative server estimate).
Both values are overridable via config / CLI.
"""

import argparse
import os
import subprocess
from datetime import datetime
from pathlib import Path

import pandas as pd


def get_disk_usage(path):
    human = "unknown"
    nbytes = 0
    try:
        r = subprocess.run(["du", "-sh", path], capture_output=True, text=True)
        if r.returncode == 0:
            human = r.stdout.split()[0]
        r = subprocess.run(["du", "-sb", path], capture_output=True, text=True)
        if r.returncode == 0:
            nbytes = int(r.stdout.split()[0])
    except Exception:
        pass
    return human, nbytes


def read_benchmarks(files):
    total_cpu_s = total_wall_s = 0.0
    n_jobs = 0
    failed = []
    for f in files:
        p = Path(f)
        if not p.exists() or p.stat().st_size == 0:
            continue
        try:
            df = pd.read_csv(p, sep="\t")
            if "cpu_time" in df.columns:
                total_cpu_s += df["cpu_time"].fillna(0).sum()
            if "s" in df.columns:
                total_wall_s += df["s"].fillna(0).sum()
            n_jobs += len(df)
        except Exception as e:
            failed.append(f"{p.name}: {e}")
    return total_cpu_s, total_wall_s, n_jobs, failed


def format_duration(seconds):
    h = int(seconds // 3600)
    m = int((seconds % 3600) // 60)
    s = seconds % 60
    return f"{h}h {m}m {s:.0f}s"


def main():
    ap = argparse.ArgumentParser(description="Generate computer_usage.log")
    ap.add_argument("--outdir",           required=True)
    ap.add_argument("--benchmarks",       nargs="*", default=[])
    ap.add_argument("--cores",            type=int,   default=1)
    ap.add_argument("--carbon-intensity", type=float, default=475.0,
                    help="gCO2/kWh (default: IEA 2022 world avg = 475)")
    ap.add_argument("--tdp",              type=float, default=10.0,
                    help="CPU TDP per core in watts (default: 10)")
    ap.add_argument("--version",          default="unknown")
    args = ap.parse_args()

    report_dir = Path(args.outdir) / "report"
    report_dir.mkdir(parents=True, exist_ok=True)

    disk_human, disk_bytes = get_disk_usage(args.outdir)
    total_cpu_s, total_wall_s, n_jobs, bench_errors = read_benchmarks(args.benchmarks)

    cpu_hours      = total_cpu_s  / 3600
    energy_kwh     = cpu_hours * args.tdp / 1000
    co2_g          = energy_kwh * args.carbon_intensity

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_path = report_dir / "computer_usage.log"

    with open(log_path, "w") as fh:
        fh.write(f"annotseba v{args.version} — Resource & Carbon Usage\n")
        fh.write(f"Generated : {now}\n")
        fh.write("=" * 60 + "\n\n")

        fh.write("[ Disk Usage ]\n")
        fh.write(f"  Output directory  : {args.outdir}\n")
        fh.write(f"  Total space used  : {disk_human}  ({disk_bytes:,} bytes)\n\n")

        fh.write("[ Compute Time  (from Snakemake benchmarks) ]\n")
        fh.write(f"  Jobs benchmarked  : {n_jobs}\n")
        fh.write(f"  Total CPU time    : {format_duration(total_cpu_s)}  ({total_cpu_s:.1f} s)\n")
        fh.write(f"  Total wall time   : {format_duration(total_wall_s)}  ({total_wall_s:.1f} s)\n\n")

        fh.write("[ Carbon Footprint (estimate) ]\n")
        fh.write(f"  Model             : CPU_time(h) x TDP(W/core) / 1000 x carbon_intensity\n")
        fh.write(f"  CPU TDP assumed   : {args.tdp:.1f} W / core\n")
        fh.write(f"  Carbon intensity  : {args.carbon_intensity:.0f} gCO2/kWh\n")
        fh.write(f"  Energy consumed   : {energy_kwh:.4f} kWh\n")
        fh.write(f"  CO2 equivalent    : {co2_g:.2f} g  ({co2_g/1000:.4f} kg)\n\n")

        fh.write("[ Notes ]\n")
        fh.write("  This is an estimate. Actual emissions depend on hardware\n")
        fh.write("  efficiency and local electricity source.  For location-\n")
        fh.write("  specific carbon intensity see https://app.electricitymaps.com\n")
        fh.write("  Override defaults with carbon_intensity and cpu_tdp_per_core\n")
        fh.write("  in config.yaml.\n")
        if bench_errors:
            fh.write("\n[ Benchmark parse warnings ]\n")
            for e in bench_errors:
                fh.write(f"  {e}\n")

    print(f"Usage log -> {log_path}")


if __name__ == "__main__":
    main()
