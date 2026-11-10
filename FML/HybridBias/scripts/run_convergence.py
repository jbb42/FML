#!/usr/bin/env python3
"""Convergence in a 1024 Mpc/h box: force mesh f x particles n x time steps s, 10 seeds each, compared with the 2048^3
reference run. Writes figures/convergence/all_spectra.pdf and figures/convergence/P_ij.pdf: the spectra (mean and
standard deviation over the seeds), then divided by the best setup.

    python3 scripts/run_convergence.py [--dry-run] [--plot-only]
"""
from plot_page import save_pages
from plot_single import save_single
from runs import run, summary
from spectra import ROOT, load_run, mean_curve

CAMPAIGN = "convergence_runs_v2_fixedseeds"
SETUPS = [(f, n, s) for f in [512, 1024] for n in [512, 1024] for s in [10, 20, 40]] + [(1280, 1280, 40)]
VERSIONS = range(1, 11)  # Seeds 1231, ..., 1240

folders = {}  # (f, n, s) -> run folders
for f, n, s in SETUPS:
    folders[(f, n, s)] = [run(f"{CAMPAIGN}/f{f}_n{n}_s{s}_v{version}", ntasks=128 if n > 1024 else 64,
                              force_nmesh=f, particle_Npart_1D=n, timestep_nsteps=[s], ic_random_seed=1230 + version)
                          for version in VERSIONS]
summary()

# Also in the plots: a 2048^3 run with 40 steps and the 2048^3 reference, both run by hand
folders[(2048, 2048, 40)] = [f"{ROOT}/results/{CAMPAIGN}/f2048_n2048_s40_v1"]
spectra = {setup: curve for setup, runs in folders.items() if (curve := mean_curve(runs))}
best = max(spectra)  # Highest f, then n, then s
ratios = {setup: {**curve, "P": curve["P"] / spectra[best]["P"], "error": None} for setup, curve in spectra.items()}
k, P = load_run(f"{ROOT}/results/reference_2048")
reference = {"label": "f = n = 2048 (reference)", "k": k, "P": P, "color": "black", "ls": "-"}

names = ["f", "n", "s"]
figures = [{"title": "", "curves": spectra, "names": names, "extra": [reference]},
           {"title": f"relative to f = {best[0]}, n = {best[1]}, s = {best[2]}", "curves": ratios, "names": names,
            "ratio": True}]
save_pages(figures, "convergence/all_spectra.pdf")
save_single(figures, "convergence")
