#!/usr/bin/env python3
"""Plot the f(R) boost P_ij^f(R) / P_ij^GR of all 15 bias spectra, saved as a PDF.

Usage: python3 scripts/plot/plot_fofr_boost.py [campaign folder] [redshift] [filters ...]
       filters: setup values (L=1024, f=1024, n=1024, s=30), f(R) models (F5) and seed=1001, e.g.
       plot_fofr_boost.py results/weekend 0 L=1024 s=30   or   plot_fofr_boost.py results/weekend 0 F5 seed=1001
       (default: results/gr_vs_fofr at z = 0, all setups). Unfinished runs are skipped. For each seed, GR and
       f(R) are averaged over the phases both have finished; with several seeds the band is the standard error.
"""
import glob
import os
import re
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
campaign = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results", "gr_vs_fofr")
z = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
filters = sys.argv[3:]
seed_filter = [f.split("=")[1] for f in filters if f.startswith("seed=")]
model_filter = [f for f in filters if f.startswith("F")]  # GR is always kept as the reference
setup_filter = [f for f in filters if "=" in f and not f.startswith("seed=")]
PAIRS = [(i, j) for i in range(5) for j in range(i, 5)]
KMIN, KMAX = 0.1, 1.0  # plotted k range in h/Mpc
NAMES = [r"1", r"\delta_L", r"\delta_L^2", r"s^2", r"\nabla^2\delta_L"]

# spectra[(model, setup)][(seed, reversed)] = {pair: [k, P(k)]} for every finished run
spectra = defaultdict(dict)
for param_file in glob.glob(os.path.join(campaign, "*", "parameterfile.lua")):
    snap = glob.glob(os.path.join(os.path.dirname(param_file), f"snapshot_*_z{z:.3f}"))
    if not snap or not os.path.exists(os.path.join(snap[0], "pofk_bias_info.txt")):
        continue
    # First assignment of each key in the Lua file, without trailing comments and quotes
    p = dict(re.findall(r'^\s*(\w+)\s*=\s*"?(.*?)"?\s*(?:--.*)?$', open(param_file).read(), re.M)[::-1])
    model = "GR" if p["gravity_model"] == "GR" else f"F{-np.log10(float(p['gravity_model_fofr_fofr0'])):.0f}"
    nsteps = sum(map(int, re.findall(r"\d+", p["timestep_nsteps"])))
    setup = f"L={p['simulation_boxsize']}, f={p['force_nmesh']}, n={p['particle_Npart_1D']}, s={nsteps}"
    if not all(f in setup.split(", ") for f in setup_filter) or (seed_filter and p["ic_random_seed"] not in seed_filter) \
            or (model_filter and model != "GR" and model not in model_filter):
        continue
    run = (p["ic_random_seed"], p["ic_reverse_phases"])
    spectra[(model, setup)][run] = {ij: np.loadtxt(f"{snap[0]}/pofk_{ij[0]}{ij[1]}.txt") for ij in PAIRS}

# 16:10 figure with fixed, equal spacing between all panels
fig, axes = plt.subplots(5, 5, figsize=(16, 10), sharex=True)
fig.subplots_adjust(left=0.05, right=0.99, bottom=0.07, top=0.97, wspace=0.3, hspace=0.35)
for ax in axes.flat:
    ax.axis("off")
for n, (model, setup) in enumerate(sorted(k for k in spectra if k[0] != "GR" and ("GR", k[1]) in spectra)):
    runs = sorted(set(spectra[("GR", setup)]) & set(spectra[(model, setup)]))  # same seeds and phases
    seeds = sorted({seed for seed, _ in runs})
    k = spectra[("GR", setup)][runs[0]][(0, 0)][:, 0]
    sel = (k >= KMIN) & (k <= KMAX)  # plot only this range, so the y axis scales to it
    for i, j in PAIRS:
        # Boost per seed (GR and f(R) averaged over that seed's phases), then mean and standard error over seeds
        boosts = []
        for seed in seeds:
            pg, pf = [np.mean([spectra[(m, setup)][r][(i, j)][:, 1] for r in runs if r[0] == seed], axis=0)
                      for m in ("GR", model)]
            boosts.append(np.where(np.sign(pf) == np.sign(pg), pf / pg, np.nan))  # meaningless where signs differ
        boosts = np.array(boosts)
        with np.errstate(invalid="ignore"), np.testing.suppress_warnings() as sup:
            sup.filter(RuntimeWarning)
            mean = np.nanmean(boosts, axis=0)
            err = np.nanstd(boosts, axis=0, ddof=1) / np.sqrt(len(seeds)) if len(seeds) > 1 else None
        ax = axes[i, j]
        ax.axis("on")
        ax.plot(k[sel], mean[sel], color=f"C{n}", label=f"{model}: {setup} ({len(seeds)} seeds)")
        if err is not None:
            ax.fill_between(k[sel], (mean - err)[sel], (mean + err)[sel], color=f"C{n}", alpha=0.25, lw=0)
        ax.axhline(1, color="k", lw=0.8)
        ax.set_xscale("log")
        ax.set_xlim(KMIN, KMAX)
        ax.xaxis.set_major_locator(FixedLocator([0.1, 0.2, 0.5, 1.0]))
        ax.xaxis.set_major_formatter(lambda x, _: f"{x:g}")
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_title(rf"${NAMES[i]} \times {NAMES[j]}$", fontsize=10, pad=3)
        ax.tick_params(labelbottom=(i == j), labelsize=8)  # x labels only on the bottom panel of each column

# y range of each panel from the 2-98th percentile of its curves, so single near-zero spikes don't flatten it
for i, j in PAIRS:
    y = np.concatenate([line.get_ydata() for line in axes[i, j].get_lines()[::2]])
    y = y[np.isfinite(y)]
    if y.size:
        lo, hi = np.percentile(y, [2, 98])
        lo, hi = min(lo, 1.0), max(hi, 1.0)
        axes[i, j].set_ylim(lo - 0.05 * (hi - lo), hi + 0.05 * (hi - lo))

# Legend in the empty lower-left triangle, outside the layout so it does not push the panels apart
fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="lower left", bbox_to_anchor=(0.05, 0.07), fontsize=10)
fig.supxlabel(r"$k\ [h/\mathrm{Mpc}]$", y=0.01)
fig.supylabel(r"$P_{ij}^{f(R)} / P_{ij}^{\rm GR}$", x=0.01)
tag = "".join("_" + f.replace("=", "") for f in filters)
out = os.path.join(ROOT, "figures", "fofr_boost", f"boost_{os.path.basename(os.path.normpath(campaign))}_z{z:.3f}{tag}.pdf")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out)
print(f"Saved {out}")
