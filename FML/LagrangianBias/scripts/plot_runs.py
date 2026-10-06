#!/usr/bin/env python3
"""Plot the 15 spectra of any runs: one line per group of runs, averaged over the group with a band of +- 1 standard
deviation if it has several runs. One page per redshift, and per 5 groups if there are more.

Usage: python3 scripts/plot_runs.py OUTPUT GROUP [GROUP ...] [--ratio] [--z Z ...]
    OUTPUT    PDF to write in figures/, e.g. ensemble/all_spectra.pdf
    GROUP     a run folder, or a quoted glob pattern matching several, relative to the LagrangianBias folder
    --ratio   divide every group by the first (interpolated to the same k)
    --z       redshifts to plot (default 0)
Examples:
    plot_runs.py ensemble/all_spectra.pdf 'results/ensemble_100seeds/*'
    plot_runs.py pofk_nmesh/all_spectra.pdf results/convergence_runs_v1/f1024_n1024_s40_v4 output --ratio
    plot_runs.py resolution_scan/all_spectra.pdf results/resolution_scan/f512_b1024_n512 results/resolution_scan/f1024_b1024_n1024
"""
import argparse

import numpy as np

from spectra import load_runs, plot_spectra, save

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("output")
parser.add_argument("groups", nargs="+")
parser.add_argument("--ratio", action="store_true")
parser.add_argument("--z", type=float, nargs="+", default=[0.0])
args = parser.parse_args()

figures = []
for z in args.z:
    curves = []
    for group in args.groups:
        k, P = load_runs(group, z)
        if k is None:
            print(f"No finished runs in {group} at z = {z}")
            continue
        error = P.std(axis=0) if len(P) > 1 else None
        curves.append({"label": f"{group} ({len(P)} runs)", "k": k, "P": P.mean(axis=0), "error": error})
    if args.ratio:
        reference = curves.pop(0)
        for curve in curves:
            P_ref = np.array([np.interp(curve["k"], reference["k"], P) for P in reference["P"]])
            curve["P"] = curve["P"] / P_ref
            curve["error"] = None if curve["error"] is None else curve["error"] / np.abs(P_ref)
    title = f"z = {z}" + (f", relative to {reference['label']}" if args.ratio else "")
    for start in range(0, len(curves), 5):  # At most 5 lines per page
        figures.append(plot_spectra(curves[start:start + 5], title, ratio=args.ratio))
save(figures, args.output)
