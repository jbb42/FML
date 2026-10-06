"""Scan parameters of parameterfile.lua, with both phases and the given seeds, then plot the spectra. A scan is a small
file calling scan(), see scan_cosmology.py, run as
    python3 scripts/scan_<name>.py [--dry-run] [--plot-only]
Runs go to results/scans/<name>/<values>_seed<seed>_<phase>/, and finished ones are skipped, so a scan can be restarted
or extended with more values or seeds. figures/scans/<name>.pdf has, for each output redshift, pages varying one
parameter (mean over seeds and phases, with a 1 sigma band), each followed by the same divided by the fiducial run
(or, if there is none, by the first value).
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
    """The value of a parameter in parameterfile.lua, as a number."""
    return float(base_parameter(parameter))


def around_fiducial(parameter, fractions):
    """The fiducial value changed by these fractions, e.g. fractions [-0.5, 0.1] give fiducial x 0.5 and x 1.1."""
    return [float(f"{fiducial(parameter) * (1 + fraction):.6g}") for fraction in fractions]


def scan(name, vary=None, one_at_a_time=None, fixed={}, seeds=None, phases=(False, True), ntasks=64):
    """Give either
        vary: {parameter: [values]}, to run every combination of the values, or
        one_at_a_time: {parameter: [values]}, to change one parameter at a time while the others keep their fiducial
                       value (from parameterfile.lua), plus one fiducial run.
    fixed: {parameter: value} for all runs; seeds: the random seeds (default the one in parameterfile.lua); phases:
    False for normal, True for reversed phases. When cosmology_* parameters change, the input P(k) of each cosmology is
    computed with CLASS (input_power.py)."""
    if vary:
        combinations = list(itertools.product(*vary.values()))
        reference = None
    else:
        vary = one_at_a_time
        reference = tuple(fiducial(parameter) for parameter in vary)
        combinations = [reference] + [reference[:i] + (value,) + reference[i + 1:]
                                      for i, values in enumerate(vary.values()) for value in values if value != reference[i]]
    runs = {}  # Values of the varied parameters -> their run folders
    for values in combinations:
        label = "_".join(f"{parameter}{value}" for parameter, value in zip(vary, values))
        parameters = {**fixed, **dict(zip(vary, values))}
        cosmology = {p: v for p, v in parameters.items() if p.startswith("cosmology_")}
        if cosmology:  # A new input P(k), computed once per cosmology
            cosmology_label = "_".join(f"{parameter}{value}" for parameter, value in sorted(cosmology.items()))
            input_file = os.path.join(ROOT, "results", "scans", name, "input", f"pofk_{cosmology_label}.txt")
            # FML needs P(k) up to the corner of the IC grid (ic_nmesh = particle_Npart_1D), sqrt(3) k_Nyquist
            npart, box = (float(parameters.get(p, base_parameter(p))) for p in ("particle_Npart_1D", "simulation_boxsize"))
            write_input_power(parameters, input_file, kmax=3**0.5 * np.pi * npart / box)
            parameters.update(ic_type_of_input="powerspectrum", ic_input_filename=input_file)
        runs[values] = []
        for seed in seeds or [int(base_parameter("ic_random_seed"))]:
            for reverse in phases:
                run_name = f"scans/{name}/{label}_seed{seed}_{'reversed' if reverse else 'normal'}"
                runs[values].append(os.path.join(ROOT, "results", run_name))
                if "--plot-only" not in sys.argv:
                    run(run_name, ntasks, **parameters, ic_random_seed=seed, ic_reverse_phases=reverse)
    summary()
    if "--dry-run" not in sys.argv:
        plot(name, list(vary), runs, reference)


def plot(name, parameters, runs, reference=None):
    """figures/scans/<name>.pdf from the runs that have finished, with ratios to the runs with the reference values."""
    snapshots = glob.glob(os.path.join(ROOT, "results", "scans", name, "*", "snapshot_*_z*"))
    figures = []
    for z in sorted({float(snapshot.rsplit("_z", 1)[1]) for snapshot in snapshots}, reverse=True):
        curves = {}
        for values, run_dirs in runs.items():
            finished = [load_run(run_dir, z) for run_dir in run_dirs if load_run(run_dir, z)[0] is not None]
            if finished:
                P = np.array([P for k, P in finished])
                curves[values] = {"k": finished[0][0], "P": P.mean(axis=0), "info": f" ({len(P)} runs)",
                                  "error": P.std(axis=0) if len(P) > 1 else None}
        for title, page in pages_varying_one_parameter(curves, parameters):
            figures.append(plot_spectra(page, f"{name}, z = {z}: {title}"))
            first = next((curve for curve in page if curve["values"] == reference), page[0])
            ratios = [{**curve, "error": None, "P": curve["P"] / [np.interp(curve["k"], first["k"], P) for P in first["P"]]}
                      for curve in page if curve is not first]
            if ratios:
                figures.append(plot_spectra(ratios, f"{name}, z = {z}: {title}, relative to {first['label']}", ratio=True))
    if figures:
        save(figures, f"scans/{name}.pdf")
