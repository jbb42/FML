"""Scan parameters of parameterfile.lua, with both phases and the given seeds, then plot the spectra. A scan is a small
file calling scan(), see scan_cosmology.py, run as
    python3 scripts/scan_<name>.py [--dry-run] [--plot-only]
Runs go to results/scans/<name>/<values>[_<gravity>]_seed<seed>_<phase>/, where <values> are the parameters that
differ from the fiducial (or "fiducial"), and finished ones are skipped, so a scan can be restarted or extended with more
values, seeds or gravity models. figures/scans/<name>.pdf has, for each output redshift, pages varying one parameter
(mean over seeds and phases, with a 1 sigma band), each followed by the same divided by the fiducial run (or, if there
is none, by the first value), and the f(R)/GR boosts if gravity models are compared.
"""
import glob
import itertools
import os
import sys

import numpy as np

from input_power import write_input_power
from simulations import ROOT, base_parameter, run, summary
from spectra import load_run, pages_varying_one_parameter, plot_spectra, save


def fiducial(parameter):
    """The value of a parameter in parameterfile.lua: a number, or a string such as "GR"."""
    value = base_parameter(parameter)
    try:
        return float(value)
    except ValueError:
        return value


def around_fiducial(parameter, fractions):
    """The fiducial value changed by these fractions, e.g. fractions [-0.5, 0.1] give fiducial x 0.5 and x 1.1."""
    return [float(f"{fiducial(parameter) * (1 + fraction):.6g}") for fraction in fractions]


def short(parameter):
    """The parameter name as used in file names and labels: cosmology_OmegaCDM -> OmegaCDM."""
    return parameter.removeprefix("cosmology_")


def gravity_parameters(model):
    """The parameters of a gravity model: "GR", or "F5" for f(R) with fofr0 = 1e-5."""
    if model == "GR":
        return {"gravity_model": "GR"}
    return {"gravity_model": "f(R)", "gravity_model_fofr_fofr0": 10.0 ** -int(model[1:])}


def scan(name, vary=None, one_at_a_time=None, fixed={}, gravity=None, seeds=None, phases=(False, True), ntasks=64):
    """Give either
        vary: {parameter: [values]}, to run every combination of the values, or
        one_at_a_time: {parameter: [values]}, to change one parameter at a time while the others keep their fiducial
                       value (from fixed, or else parameterfile.lua), plus one fiducial run.
    fixed: {parameter: value} for all runs; seeds: the random seeds (default the one in parameterfile.lua); phases:
    False for normal, True for reversed phases. When cosmology_* parameters change, the input P(k) of each cosmology is
    computed with CLASS (input_power.py).
    gravity: f(R) models, e.g. ["F5"] or ["F4", "F5", "F6"] (fofr0 = 1e-4, ...): every combination is then also run in
    these, and in GR, all from the same initial conditions (ic_use_gravity_model_GR = true), and the figure gets the
    boosts f(R)/GR."""
    if vary:
        combinations = list(itertools.product(*vary.values()))
        reference = None
    else:
        vary = one_at_a_time
        reference = tuple(fixed.get(parameter, fiducial(parameter)) for parameter in vary)
        combinations = [reference] + [reference[:i] + (value,) + reference[i + 1:]
                                      for i, values in enumerate(vary.values()) for value in values if value != reference[i]]
    models = ["GR"] + gravity if gravity else [None]
    runs = {}  # (values of the varied parameters, gravity model) -> their run folders
    for values in combinations:
        changed = {p: v for p, v, ref in zip(vary, values, reference or [None] * len(vary)) if v != ref}
        label = "_".join(f"{short(p)}{v}" for p, v in changed.items()) or "fiducial"
        parameters = {**fixed, **dict(zip(vary, values))}
        if any(p.startswith("cosmology_") for p in parameters):  # A new input P(k), computed once per cosmology
            cosmology_label = "_".join(f"{short(p)}{v}" for p, v in changed.items() if p.startswith("cosmology_"))
            input_file = os.path.join(ROOT, "results", "scans", name, "input", f"pofk_{cosmology_label or 'fiducial'}.txt")
            # FML needs P(k) up to the corner of the IC grid (ic_nmesh = particle_Npart_1D), sqrt(3) k_Nyquist
            npart, box = (float(parameters.get(p, base_parameter(p))) for p in ("particle_Npart_1D", "simulation_boxsize"))
            write_input_power(parameters, input_file, kmax=3**0.5 * np.pi * npart / box)
            parameters.update(ic_type_of_input="powerspectrum", ic_input_filename=input_file)
        for model in models:
            model_parameters = {**gravity_parameters(model), "ic_use_gravity_model_GR": True} if model else {}
            runs[(values, model)] = []
            for seed in seeds or [int(base_parameter("ic_random_seed"))]:
                for reverse in phases:
                    run_name = f"scans/{name}/{label}{'_' + model if model else ''}_seed{seed}_{'reversed' if reverse else 'normal'}"
                    runs[(values, model)].append(os.path.join(ROOT, "results", run_name))
                    if "--plot-only" not in sys.argv:
                        run(run_name, ntasks, **parameters, **model_parameters, ic_random_seed=seed, ic_reverse_phases=reverse)
    summary()
    if "--dry-run" not in sys.argv:
        plot(name, [short(p) for p in vary], runs, reference, models)


def seed_boosts(run_dirs, gr, fr):
    """[(k, boost)] for each seed: the f(R) spectra divided by the GR ones, each averaged over the phases that finished
    in both. Averaging the phases first cancels the leading cosmic variance, also in spectra that are odd in the
    initial field (which a ratio per phase would not)."""
    by_seed = {}
    for run_dir, g, f in zip(run_dirs, gr, fr):
        if g[0] is not None and f[0] is not None:
            by_seed.setdefault(run_dir.split("_seed")[1].split("_")[0], []).append((g, f))
    return [(pairs[0][0][0], np.mean([f[1] for g, f in pairs], axis=0) / np.mean([g[1] for g, f in pairs], axis=0))
            for pairs in by_seed.values()]


def plot(name, parameters, runs, reference, models):
    """figures/scans/<name>.pdf from the runs that have finished: per redshift, the spectra (in GR if gravity models
    are compared) with ratios to the reference values, then the boost of each f(R) model."""
    snapshots = glob.glob(os.path.join(ROOT, "results", "scans", name, "*", "snapshot_*_z*"))
    figures = []
    for z in sorted({float(snapshot.rsplit("_z", 1)[1]) for snapshot in snapshots}, reverse=True):
        spectra = {model: {} for model in models}  # model -> {values: [(k, P) of each finished run]}
        for (values, model), run_dirs in runs.items():
            spectra[model][values] = [load_run(run_dir, z) for run_dir in run_dirs]

        curves = {}  # The spectra of the first model (GR, or the only one): mean over seeds and phases
        for values, finished in spectra[models[0]].items():
            P = np.array([P for k, P in finished if k is not None])
            if len(P):
                k = next(k for k, P in finished if k is not None)
                curves[values] = {"k": k, "P": P.mean(axis=0), "error": P.std(axis=0) if len(P) > 1 else None,
                                  "info": f" ({len(P)} runs)"}
        for title, page in pages_varying_one_parameter(curves, parameters):
            figures.append(plot_spectra(page, f"{name}, z = {z}: {title}"))
            first = next((curve for curve in page if curve["values"] == reference), page[0])
            ratios = [{**curve, "error": None, "P": curve["P"] / [np.interp(curve["k"], first["k"], P) for P in first["P"]]}
                      for curve in page if curve is not first]
            if ratios:
                figures.append(plot_spectra(ratios, f"{name}, z = {z}: {title}, relative to {first['label']}", ratio=True))

        for model in models[1:]:  # Boosts f(R)/GR
            boosts = {}
            for values in spectra[model]:
                boost = seed_boosts(runs[(values, "GR")], spectra["GR"][values], spectra[model][values])
                if boost:
                    k, boost = boost[0][0], np.array([b for k, b in boost])
                    boosts[values] = {"k": k, "P": boost.mean(axis=0), "info": f" ({len(boost)} seeds)",
                                      "error": boost.std(axis=0, ddof=1) / np.sqrt(len(boost)) if len(boost) > 1 else None}
            for title, page in pages_varying_one_parameter(boosts, parameters):
                figures.append(plot_spectra(page, f"{name}, z = {z}: {model} / GR, {title}", ratio=True))
    if figures:
        save(figures, f"scans/{name}.pdf")
