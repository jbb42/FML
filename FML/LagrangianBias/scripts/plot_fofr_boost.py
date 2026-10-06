#!/usr/bin/env python3
"""The f(R) boost P_ij^f(R) / P_ij^GR of the 15 spectra for a campaign of GR and f(R) runs with the same seeds.

Writes figures/fofr_boost/all_boost_<campaign>.pdf. For each redshift there are pages where one parameter varies: the
f(R) model (F4 = fR0 1e-4, ...), box size L, force mesh f, particles n or time steps s. For each seed, GR and f(R) are
averaged over the phases (normal/reversed) both have finished; a line is the mean boost over the seeds and, with
several seeds, the band is its standard error. Unfinished runs are skipped, so this also works during a campaign.

Usage: python3 scripts/plot_fofr_boost.py [campaign]    (a folder in results/, default gr_vs_fofr)
"""
import os
import re
import sys
import warnings

import numpy as np

from spectra import ROOT, load_run, pages_varying_one_parameter, plot_spectra, save

campaign = sys.argv[1] if len(sys.argv) > 1 else "gr_vs_fofr"
campaign_dir = os.path.join(ROOT, "results", campaign)


def read_parameters(path):
    """{name: value} of the 'name = value' lines of a Lua file, values without quotes and comments (first one wins)."""
    parameters = {}
    for name, value in re.findall(r'^\s*(\w+)\s*=\s*"?(.*?)"?\s*(?:--.*)?$', open(path).read(), re.M):
        parameters.setdefault(name, value)
    return parameters


# Every run, with its model ("GR", "F4", ...), setup (L, f, n, s), seed, phase and output redshifts
runs = []
for run in sorted(os.listdir(campaign_dir)):
    if os.path.exists(os.path.join(campaign_dir, run, "parameterfile.lua")):
        p = read_parameters(os.path.join(campaign_dir, run, "parameterfile.lua"))
        model = "GR" if p["gravity_model"] == "GR" else f"F{-np.log10(float(p['gravity_model_fofr_fofr0'])):.0f}"
        nsteps = sum(int(steps) for steps in re.findall(r"\d+", p["timestep_nsteps"]))  # Summed over output intervals
        setup = (int(float(p["simulation_boxsize"])), int(p["force_nmesh"]), int(p["particle_Npart_1D"]), nsteps)
        runs.append({"dir": os.path.join(campaign_dir, run), "model": model, "setup": setup, "seed": p["ic_random_seed"],
                     "reversed": p["ic_reverse_phases"], "redshifts": [float(z) for z in re.findall(r"[\d.]+", p["output_redshifts"])]})


def finished_runs(model, setup, z):
    """{(seed, reversed): spectra} of the runs of this model and setup that have finished redshift z."""
    found = {}
    for r in runs:
        if r["model"] == model and r["setup"] == setup and load_run(r["dir"], z)[0] is not None:
            found[(r["seed"], r["reversed"])] = load_run(r["dir"], z)
    return found


def boost(model, setup, z):
    """Mean boost over the seeds that finished both GR and f(R), with its standard error, or None if there are none."""
    gr, fr = finished_runs("GR", setup, z), finished_runs(model, setup, z)
    both = set(gr) & set(fr)  # (seed, reversed) finished for both
    seeds = sorted({seed for seed, _ in both})
    if not seeds:
        return None
    boosts = []
    for seed in seeds:
        P_gr = np.mean([gr[run][1] for run in both if run[0] == seed], axis=0)  # Averaged over the phases
        P_fr = np.mean([fr[run][1] for run in both if run[0] == seed], axis=0)
        boosts.append(np.where(np.sign(P_fr) == np.sign(P_gr), P_fr / P_gr, np.nan))  # Meaningless where signs differ
    with warnings.catch_warnings():  # k where every seed is NaN give NaN, which is fine
        warnings.simplefilter("ignore", RuntimeWarning)
        mean = np.nanmean(boosts, axis=0)
        error = np.nanstd(boosts, axis=0, ddof=1) / np.sqrt(len(seeds)) if len(seeds) > 1 else None
    k = next(iter(gr.values()))[0]
    return {"k": k, "P": mean, "error": error, "info": f" ({len(seeds)} seeds)"}


figures = []
for z in sorted({z for r in runs for z in r["redshifts"]}):
    curves = {}  # (model, L, f, n, s) -> boost
    for model, setup in sorted({(r["model"], r["setup"]) for r in runs if r["model"] != "GR"}):
        curve = boost(model, setup, z)
        if curve is not None:
            curves[(model, *setup)] = curve
    for title, page in pages_varying_one_parameter(curves, ["model", "L", "f", "n", "s"]):
        figures.append(plot_spectra(page, f"f(R) / GR at z = {z}: {title}", ratio=True, kmin=0.1, kmax=1.0))
save(figures, f"fofr_boost/all_boost_{campaign}.pdf")
