#!/usr/bin/env python3
"""Scan one or more parameters: run every combination (with both phases and any number of seeds), then plot.

Usage: python3 scripts/scan.py NAME PARAMETER=VALUES... [options]
    NAME                 the scan; runs go to results/scans/NAME/, the figure to figures/scans/NAME.pdf
    PARAMETER=VALUES     a parameter of parameterfile.lua and comma-separated values, e.g. cosmology_OmegaCDM=0.2,0.3,0.4
    --set PARAMETER=VALUE ...  parameters to change for every run, e.g. --set force_nmesh=1024 timestep_nsteps={30}
    --seeds SEED ...     random seeds (default: the one in parameterfile.lua)
    --normal-only        only normal phases (default: normal and reversed, which cancels most cosmic variance)
    --ntasks N           MPI tasks per run (default 64)
    --plot-only          only (re)make the figure from the runs that have finished
    --dry-run            only write the parameter files and print the commands
Examples:
    python3 scripts/scan.py omega_m cosmology_OmegaCDM=0.2,0.3,0.4 --seeds 1001 1002
    python3 scripts/scan.py mesh force_nmesh=512,1024,1280 particle_Npart_1D=512,1024

Runs are skipped once finished, so a scan can be restarted or extended with more values or seeds. When cosmology_*
parameters change, the input P(k) is computed for each cosmology with CAMB or CLASS, see input_power.py. The figure has,
for each output redshift, pages varying one parameter (mean over seeds and phases, with a +- 1 sigma band if there
are several), and the same divided by the first value of that parameter.
"""
import argparse
import glob
import itertools
import os
import sys

import numpy as np

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [SCRIPTS, os.path.join(SCRIPTS, "run"), os.path.join(SCRIPTS, "plot")]
from simulations import ROOT, base_parameter, run, summary  # noqa: E402
from spectra import load_run, pages_varying_one_parameter, plot_spectra, save  # noqa: E402
from input_power import write_input_power  # noqa: E402


def parse_value(text):
    """Command line text as a Python value: "1024" -> 1024, "0.3" -> 0.3, "true" -> True, "{30}" -> [30], else text."""
    if text.startswith("{") and text.endswith("}"):
        return [parse_value(t.strip()) for t in text[1:-1].split(",")]
    if text in ("true", "false"):
        return text == "true"
    for convert in (int, float):
        try:
            return convert(text)
        except ValueError:
            pass
    return text


def short_name(parameter):
    """The parameter without its group prefix, for run names and labels: cosmology_OmegaCDM -> OmegaCDM."""
    return parameter.split("_", 1)[1] if parameter.startswith(("cosmology_", "simulation_", "particle_")) else parameter


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("name")
parser.add_argument("scan", nargs="+", metavar="PARAMETER=VALUES")
parser.add_argument("--set", nargs="+", default=[], metavar="PARAMETER=VALUE")
parser.add_argument("--seeds", nargs="+", type=int, default=[int(base_parameter("ic_random_seed"))])
parser.add_argument("--normal-only", action="store_true")
parser.add_argument("--ntasks", type=int, default=64)
parser.add_argument("--plot-only", action="store_true")
parser.add_argument("--dry-run", action="store_true")
args = parser.parse_args()

# The cosmology_* parameters that set the input P(k); only LCDM is set up
COSMOLOGY = ["cosmology_h", "cosmology_Omegab", "cosmology_OmegaCDM", "cosmology_OmegaMNu", "cosmology_OmegaK",
             "cosmology_Neffective", "cosmology_TCMB_kelvin", "cosmology_As", "cosmology_ns", "cosmology_kpivot_mpc"]
scanned = {p.split("=")[0]: [parse_value(v) for v in p.split("=", 1)[1].split(",")] for p in args.scan}
fixed = {p.split("=")[0]: parse_value(p.split("=", 1)[1]) for p in args.set}
phases = [False] if args.normal_only else [False, True]

# Run every combination of the scanned values, seeds and phases in results/scans/NAME/<values>_seed<seed>_<phase>
runs = {}  # Scanned values -> the run folders with those values
for values in itertools.product(*scanned.values()):
    parameters = {**fixed, **dict(zip(scanned, values))}
    label = "_".join(f"{short_name(p)}{v}" for p, v in zip(scanned, values))
    if any(p.startswith("cosmology_") for p in parameters):
        input_file = os.path.join(ROOT, "results", "scans", args.name, "input", f"pofk_{label}.txt")
        if not args.plot_only and not os.path.exists(input_file):  # Computed once per cosmology
            if parameters.get("cosmology_model", base_parameter("cosmology_model")) != "LCDM":
                sys.exit("The input P(k) is only set up for cosmology_model = LCDM")
            cosmology = {name[len("cosmology_"):]: float(parameters.get(name, base_parameter(name))) for name in COSMOLOGY}
            os.makedirs(os.path.dirname(input_file), exist_ok=True)
            write_input_power(cosmology, float(parameters.get("ic_input_redshift", base_parameter("ic_input_redshift"))), input_file)
        parameters.update(ic_type_of_input="powerspectrum", ic_input_filename=input_file)
    runs[values] = []
    for seed in args.seeds:
        for reverse in phases:
            name = f"scans/{args.name}/{label}_seed{seed}_{'reversed' if reverse else 'normal'}"
            runs[values].append(os.path.join(ROOT, "results", name))
            if not args.plot_only:
                run(name, args.ntasks, **parameters, ic_random_seed=seed, ic_reverse_phases=reverse)
if not args.plot_only:
    summary()
if args.dry_run:
    sys.exit()

# For each redshift, pages varying one scanned parameter (mean over seeds and phases), then relative to its first value
names = [short_name(p) for p in scanned]
snapshots = [s for run_dirs in runs.values() for run_dir in run_dirs for s in glob.glob(f"{run_dir}/snapshot_*_z*")]
redshifts = sorted({float(s.rsplit("_z", 1)[1]) for s in snapshots}, reverse=True)
figures = []
for z in redshifts:
    curves = {}
    for values, run_dirs in runs.items():
        finished = [load_run(run_dir, z) for run_dir in run_dirs if load_run(run_dir, z)[0] is not None]
        if finished:
            P = np.array([P for k, P in finished])
            curves[values] = {"k": finished[0][0], "P": P.mean(axis=0), "error": P.std(axis=0) if len(P) > 1 else None,
                              "info": f" ({len(P)} runs)"}
    for title, page in pages_varying_one_parameter(curves, names):
        figures.append(plot_spectra(page, f"{args.name}, z = {z}: {title}"))
        reference = page[0]
        ratios = [{**c, "P": c["P"] / np.array([np.interp(c["k"], reference["k"], P) for P in reference["P"]]),
                   "error": None} for c in page[1:]]
        if ratios:
            figures.append(plot_spectra(ratios, f"{args.name}, z = {z}: {title}, relative to {reference['label']}", ratio=True))
if figures:
    save(figures, f"scans/{args.name}.pdf")
else:
    print("No finished runs to plot yet")
