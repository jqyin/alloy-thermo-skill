#!/usr/bin/env python3
"""Analyze alloy-mc output: thermodynamics + order-disorder transition.

Reads one or more `thermo_run*.csv` files produced by the engine, averages over
independent runs (with error bars), produces plots, and estimates the
order-disorder transition temperature Tc from the specific-heat peak (and,
independently, the susceptibility peak).

Usage:
    python analyze.py --dir run/ --out run/analysis
    python analyze.py run/thermo_run0.csv run/thermo_run1.csv

Outputs (in --out, default ./analysis):
    thermo.png        E, C, chi, Binder vs T
    order.png         short-range-order parameters vs T
    summary.json      averaged data + estimated Tc
    summary.md        human-readable summary
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")  # headless / HPC-friendly
import matplotlib.pyplot as plt  # noqa: E402


def load_runs(paths):
    """Load CSVs sharing a header; return (header, stack[n_runs, n_T, n_col])."""
    header = None
    runs = []
    for p in paths:
        with open(p) as f:
            cols = f.readline().strip().split(",")
        if header is None:
            header = cols
        elif cols != header:
            raise ValueError(f"column mismatch in {p}")
        data = np.loadtxt(p, delimiter=",", skiprows=1, ndmin=2)
        # Sort rows by temperature so all runs align.
        data = data[np.argsort(data[:, 0])]
        runs.append(data)
    shapes = {r.shape for r in runs}
    if len(shapes) != 1:
        raise ValueError(f"runs have differing shapes: {shapes}")
    return header, np.stack(runs, axis=0)


def peak_temperature(T, y):
    """Estimate the peak position of y(T) with a local quadratic refinement."""
    i = int(np.argmax(y))
    if 0 < i < len(T) - 1:
        # Parabolic interpolation through the three points around the max.
        x0, x1, x2 = T[i - 1], T[i], T[i + 1]
        y0, y1, y2 = y[i - 1], y[i], y[i + 1]
        denom = (x0 - x1) * (x0 - x2) * (x1 - x2)
        if abs(denom) > 0:
            A = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom
            B = (x2**2 * (y0 - y1) + x1**2 * (y2 - y0) + x0**2 * (y1 - y2)) / denom
            if A < 0:  # concave -> real maximum
                xv = -B / (2 * A)
                if x0 <= xv <= x2:
                    return float(xv), float(y[i]), i
    return float(T[i]), float(y[i]), i


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csvs", nargs="*", help="thermo_run*.csv files")
    p.add_argument("--dir", help="directory to glob thermo_run*.csv from")
    p.add_argument("--out", default="analysis", help="output directory")
    args = p.parse_args()

    paths = list(args.csvs)
    if args.dir:
        paths += sorted(glob.glob(os.path.join(args.dir, "thermo_run*.csv")))
    paths = sorted(set(paths))
    if not paths:
        sys.exit("no input CSVs found (pass files or --dir)")

    header, stack = load_runs(paths)
    n_runs = stack.shape[0]
    col = {name: k for k, name in enumerate(header)}

    mean = stack.mean(axis=0)
    # Standard error of the mean across independent runs.
    err = (stack.std(axis=0, ddof=1) / np.sqrt(n_runs)
           if n_runs > 1 else np.zeros_like(mean))

    T = mean[:, col["T"]]
    os.makedirs(args.out, exist_ok=True)

    def c(name):
        return mean[:, col[name]]

    def e(name):
        return err[:, col[name]]

    # --- Transition estimates -------------------------------------------
    Tc_C, Cmax, _ = peak_temperature(T, c("C_per_site"))
    Tc_chi, chimax, _ = peak_temperature(T, c("chi"))

    # --- Plot 1: thermodynamics -----------------------------------------
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    ax[0, 0].errorbar(T, c("E_per_site"), yerr=e("E_per_site"), fmt="o-", ms=4)
    ax[0, 0].set(xlabel="T (K)", ylabel="Energy per site", title="Energy")

    ax[0, 1].errorbar(T, c("C_per_site"), yerr=e("C_per_site"), fmt="o-", ms=4,
                      color="crimson")
    ax[0, 1].axvline(Tc_C, ls="--", color="gray")
    ax[0, 1].set(xlabel="T (K)", ylabel="C per site (units of k_B)",
                 title=f"Specific heat (peak ~ {Tc_C:.0f} K)")

    ax[1, 0].errorbar(T, c("chi"), yerr=e("chi"), fmt="o-", ms=4, color="teal")
    ax[1, 0].axvline(Tc_chi, ls="--", color="gray")
    ax[1, 0].set(xlabel="T (K)", ylabel="Susceptibility",
                 title=f"Susceptibility (peak ~ {Tc_chi:.0f} K)")

    ax[1, 1].errorbar(T, c("binder"), yerr=e("binder"), fmt="o-", ms=4,
                      color="purple")
    ax[1, 1].set(xlabel="T (K)", ylabel="Binder cumulant", title="Binder cumulant")
    for a in ax.flat:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "thermo.png"), dpi=130)
    plt.close(fig)

    # --- Plot 2: short-range order --------------------------------------
    alpha_cols = [h for h in header if h.startswith("alpha_") and h != "alpha_mean"]
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for name in alpha_cols:
        ax.plot(T, c(name), "o-", ms=3, label=name.replace("alpha_", ""))
    ax.plot(T, c("alpha_mean"), "k--", lw=2, label="mean |alpha|")
    ax.axhline(0, color="gray", lw=0.8)
    ax.set(xlabel="T (K)", ylabel="Warren-Cowley SRO parameter",
           title="Short-range order vs temperature")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "order.png"), dpi=130)
    plt.close(fig)

    # --- Summaries -------------------------------------------------------
    summary = {
        "n_runs": int(n_runs),
        "n_temperatures": int(len(T)),
        "T_range": [float(T.min()), float(T.max())],
        "Tc_from_specific_heat": Tc_C,
        "specific_heat_peak": Cmax,
        "Tc_from_susceptibility": Tc_chi,
        "susceptibility_peak": chimax,
        "columns": header,
        "T": T.tolist(),
        "mean": {h: mean[:, col[h]].tolist() for h in header},
        "sem": {h: err[:, col[h]].tolist() for h in header},
    }
    with open(os.path.join(args.out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    spread = abs(Tc_C - Tc_chi)
    consistent = spread < 0.15 * max(Tc_C, 1.0)
    md = [
        "# Alloy-MC analysis summary",
        "",
        f"- Independent runs averaged: **{n_runs}**",
        f"- Temperature replicas: **{len(T)}** spanning "
        f"{T.min():.0f}-{T.max():.0f} K",
        "",
        "## Order-disorder transition",
        f"- Specific-heat peak: **Tc ~ {Tc_C:.0f} K** (C_max = {Cmax:.4g})",
        f"- Susceptibility peak: **Tc ~ {Tc_chi:.0f} K** (chi_max = {chimax:.4g})",
        "",
    ]
    if consistent:
        md.append(f"The two estimators agree to within {spread:.0f} K, "
                  "supporting a transition near "
                  f"**{0.5 * (Tc_C + Tc_chi):.0f} K**.")
    else:
        md.append(f"The estimators differ by {spread:.0f} K. The peak may be "
                  "broad (finite-size rounding) or there may be no sharp "
                  "transition; consider larger N, more samples, or a denser "
                  "temperature ladder near the peak.")
    md += [
        "",
        "## Caveats",
        "- A single system size cannot pin Tc precisely; the specific-heat "
        "peak shifts and sharpens with N. For a rigorous Tc use Binder-cumulant "
        "crossings across several N.",
        "- Ensure the temperature ladder brackets the peak and that runs are "
        "equilibrated (check swap/move acceptance in the CSV).",
        "",
        "Figures: `thermo.png`, `order.png`. Data: `summary.json`.",
    ]
    with open(os.path.join(args.out, "summary.md"), "w") as f:
        f.write("\n".join(md) + "\n")

    print(f"Analysis written to {args.out}/")
    print(f"  Tc (specific heat)   ~ {Tc_C:.0f} K")
    print(f"  Tc (susceptibility)  ~ {Tc_chi:.0f} K")
    print("  Figures: thermo.png, order.png ; summary.json, summary.md")


if __name__ == "__main__":
    main()
