#!/usr/bin/env python3
"""Convergence of the spectra with force mesh f, particles n and time steps s, from results/convergence_runs_v2_fixedseeds
(runs f<f>_n<n>_s<s>_v<version>, averaged over the 10 versions = seeds of each setup).

Writes figures/convergence/all_spectra.pdf, with pages where one of f, n and s varies:
    first the spectra (mean and +- 1 standard deviation over the seeds), with the 2048^3 reference run in black
    then the same divided by the best setup (the highest resolution)

Usage: python3 scripts/plot/plot_convergence.py
"""
import glob
import os

from spectra import ROOT, load_runs, pages_varying_one_parameter, plot_spectra, save

CAMPAIGN = "results/convergence_runs_v2_fixedseeds"

# Setups (f, n, s) from the run folder names, leaving out f or n = 256
setups = set()
for run_dir in glob.glob(os.path.join(ROOT, CAMPAIGN, "f*_n*_s*_v*")):
    f, n, s = (int(part[1:]) for part in os.path.basename(run_dir).split("_")[:3])
    if f != 256 and n != 256:
        setups.add((f, n, s))

spectra = {}
for f, n, s in setups:
    k, P = load_runs(f"{CAMPAIGN}/f{f}_n{n}_s{s}_v*")
    spectra[(f, n, s)] = {"k": k, "P": P.mean(axis=0), "error": P.std(axis=0), "info": f" ({len(P)} seeds)"}
best = max(spectra)  # Highest f, then n, then s
ratios = {setup: {**curve, "P": curve["P"] / spectra[best]["P"], "error": None} for setup, curve in spectra.items()}
k, P = load_runs("results/reference_2048")
reference = {"label": "f = n = 2048 (reference)", "k": k, "P": P[0], "color": "black"}

names = ["f", "n", "s"]
figures = [plot_spectra(curves + [reference], title) for title, curves in pages_varying_one_parameter(spectra, names)]
figures += [plot_spectra(curves, f"{title}, relative to f = {best[0]}, n = {best[1]}, s = {best[2]}", ratio=True)
            for title, curves in pages_varying_one_parameter(ratios, names)]
save(figures, "convergence/all_spectra.pdf")
