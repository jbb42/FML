"""Load the 15 spectra of runs and combine them into figures for plot_page.py and plot_single.py.

A run folder has snapshot_*_z<redshift>/pofk_ij.txt files. Its 15 spectra P_ij(k), i <= j, are of the fields
0 = matter, 1 = delta_L, 2 = delta_L^2, 3 = s^2 and 4 = nabla^2 delta_L.

A curve is a dict with "k", "P" (array [15, len(k)]), optionally "error" (same shape, drawn as a band), "label",
"info" (added to the label) and matplotlib options such as "color".
A figure is a dict with
    "title"
    "curves": {values: curve}, with the values of the parameters "names", or a list of labelled curves
    "names", "fiducial" (the fiducial values, so labels only show what differs from them)
    "extra": curves drawn on every page, e.g. a reference
    "kind": "spectra", "relative" or "boost" for the figures of scan_figures()
    "ratio", "kmin", "kmax", "line", "ylabel", "ylim", "band": options of plot_page() and plot_single()
"""
import functools
import glob
import os
import warnings
from collections import defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # The LagrangianBias folder
PAIRS = [(i, j) for i in range(5) for j in range(i, 5)]  # The 15 spectra, in the order of the arrays
NAMES = [r"\delta_m", r"\delta_L", r"\delta_L^2", r"s^2", r"\nabla^2\delta_L"]


@functools.lru_cache(maxsize=None)  # Each run is read only once
def load_run(folder, z=0.0):
    """k and the 15 spectra of a run as an array [15, len(k)], or (None, None) if it has no output at redshift z."""
    snapshots = glob.glob(os.path.join(folder, f"snapshot_*_z{z:.3f}"))
    if not snapshots or not os.path.exists(os.path.join(snapshots[0], "pofk_44.txt")):
        return None, None
    spectra = [np.loadtxt(os.path.join(snapshots[0], f"pofk_{i}{j}.txt")) for i, j in PAIRS]
    return spectra[0][:, 0], np.array([s[:, 1] for s in spectra])


def redshifts(folders):
    """The output redshifts of these runs, highest first."""
    snapshots = [s for folder in folders for s in glob.glob(os.path.join(folder, "snapshot_*_z*"))]
    return sorted({float(s.rsplit("_z", 1)[1]) for s in snapshots}, reverse=True)


def mean_curve(folders, z=0.0):
    """The mean of the runs in these folders that have finished redshift z, with their standard deviation as the
    error if there are several, or None if none have finished."""
    finished = [(k, P) for k, P in (load_run(folder, z) for folder in folders) if k is not None]
    if not finished:
        return None
    P = np.array([P for k, P in finished])
    return {"k": finished[0][0], "P": P.mean(axis=0), "error": P.std(axis=0) if len(P) > 1 else None,
            "info": f" ({len(P)} runs)"}


def paired_ratios(pairs):
    """{values: curve} from pairs {values: {seed: [((k, P_denominator), (k, P_numerator)), ...]}}, one pair per phase.
    The ratio of the means over all phases and seeds: averaging the phases cancels the leading cosmic variance, also
    of spectra that are odd in the initial field, and averaging the seeds before dividing keeps the ratio stable where
    a single seed is close to zero. The error is a jackknife over the seeds (the ratio of the means without one seed
    at a time), so it includes that numerator and denominator share their initial modes. Where the two have opposite
    signs the ratio is meaningless and left out."""
    def ratio(numerator, denominator):
        return np.where(np.sign(numerator) == np.sign(denominator), numerator / denominator, np.nan)

    curves = {}
    for values, seeds in pairs.items():
        k = next(iter(seeds.values()))[0][1][0]
        # Per seed, the phase means [seed, spectrum, k]
        den = np.array([np.mean([[np.interp(k, d[0], P) for P in d[1]] for d, n in runs], axis=0)
                        for runs in seeds.values()])
        num = np.array([np.mean([n[1] for d, n in runs], axis=0) for runs in seeds.values()])
        n = len(num)
        with warnings.catch_warnings():  # Ratios of zero
            warnings.simplefilter("ignore", RuntimeWarning)
            error = None
            if n > 1:
                leave_one_out = np.array([ratio(np.delete(num, i, 0).mean(0), np.delete(den, i, 0).mean(0))
                                          for i in range(n)])
                error = np.sqrt((n - 1) / n * ((leave_one_out - leave_one_out.mean(0)) ** 2).sum(0))
            curves[values] = {"k": k, "P": ratio(num.mean(0), den.mean(0)), "error": error,
                              "info": f" ({n} seed{'s' if n > 1 else ''})"}
    return curves


def boosts(runs, model, z=0.0):
    """{values: f(R)/GR curve} for the runs (runs.Run) in this model, each divided by the GR run with the same values,
    seed and phase (see paired_ratios)."""
    gr = {(r.values, r.seed, r.reverse): r.folder for r in runs if r.model == "GR"}
    pairs = defaultdict(lambda: defaultdict(list))  # values -> seed -> [(GR, f(R)) spectra]
    for r in runs:
        if r.model == model:
            P_gr, P_fr = load_run(gr[(r.values, r.seed, r.reverse)], z), load_run(r.folder, z)
            if P_gr[0] is not None and P_fr[0] is not None:
                pairs[r.values][r.seed].append((P_gr, P_fr))
    return paired_ratios(pairs)


def relative(runs, model, reference, z=0.0):
    """{values: curve} for the runs (runs.Run) in this model, each divided by the run with the reference values and
    the same seed and phase (see paired_ratios). Runs with the same seed share their initial modes, so this cancels
    cosmic variance and the error shows how much the ratio changes from seed to seed."""
    ref = {(r.seed, r.reverse): r.folder for r in runs if r.values == reference and r.model == model}
    pairs = defaultdict(lambda: defaultdict(list))  # values -> seed -> [(reference, run) spectra]
    for r in runs:
        if r.model == model and (r.seed, r.reverse) in ref:
            P_ref, P = load_run(ref[(r.seed, r.reverse)], z), load_run(r.folder, z)
            if P_ref[0] is not None and P[0] is not None:
                pairs[r.values][r.seed].append((P_ref, P))
    return paired_ratios(pairs)


def scan_figures(scan):
    """The figures of a scan (runs.scan), for each redshift: the spectra (mean over seeds and phases, in GR if gravity
    models are compared), the same divided by the fiducial run (or else the first) seed by seed, and the f(R)/GR
    boosts."""
    figures = []
    for z in redshifts([r.folder for r in scan.runs]):
        title = f"{scan.name}, z = {z}"
        names = {"names": scan.parameters, "fiducial": scan.reference}
        spectra = {}
        for values in dict.fromkeys(r.values for r in scan.runs):
            curve = mean_curve([r.folder for r in scan.runs if r.values == values and r.model == scan.models[0]], z)
            if curve:
                spectra[values] = curve
        figures.append({"title": title, "kind": "spectra", "curves": spectra, **names})
        if len(spectra) > 1:
            reference = scan.reference if scan.reference in spectra else min(spectra)
            relative_to = "the fiducial" if reference == scan.reference and len(scan.parameters) > 2 else \
                ", ".join(f"{n} = {v}" for n, v in zip(scan.parameters, reference))
            figures.append({"title": f"{title}, relative to {relative_to}", "kind": "relative", "ratio": True,
                            "curves": relative(scan.runs, scan.models[0], reference, z), **names})
        for model in scan.models[1:]:
            figures.append({"title": f"{title}: {model} / GR", "kind": "boost", "ratio": True,
                            "curves": boosts(scan.runs, model, z), **names})
    return figures

def pages_varying_one_parameter(curves, names, max_lines=6):
    """Split {values: curve} into pages where only one parameter varies, with at most max_lines lines each, so they
    stay readable. Returns a list of (title, curves), each curve labelled by the value that varies on its page.
    Curves that differ from all others in more than one parameter are put on extra pages at the end."""
    pages, shown = [], set()
    for i, name in enumerate(names):
        others = names[:i] + names[i + 1:]
        groups = defaultdict(list)  # Values of the other parameters -> the curves with those
        for values in sorted(curves):
            groups[values[:i] + values[i + 1:]].append(values)
        for fixed, members in groups.items():
            if len(members) < 2:
                continue
            for start in range(0, len(members), max_lines):
                chunk = members[start:start + max_lines]
                title = ", ".join([f"Varying {name}"] + [f"{n} = {v}" for n, v in zip(others, fixed)])
                pages.append((title, [{**curves[v], "label": f"{name} = {v[i]}" + curves[v].get("info", "")}
                                      for v in chunk]))
                shown.update(chunk)
    rest = [v for v in sorted(curves) if v not in shown]
    for start in range(0, len(rest), max_lines):
        pages.append(("Other", [{**curves[v], "label": label(v, names) + curves[v].get("info", "")}
                                for v in rest[start:start + max_lines]]))
    return pages


def label(values, names, fiducial=None):
    """"f = 512, n = 1024", or with fiducial values only the ones that differ ("fiducial" if none)."""
    fiducial = fiducial or [None] * len(values)
    return ", ".join(f"{n} = {v}" for n, v, f in zip(names, values, fiducial) if v != f) or "fiducial"


def group_figures(groups, zs=(0.0,), ratio=False):
    """Figures for the command lines of plot_page.py and plot_single.py: a curve per group of run folders (a glob
    pattern relative to the LagrangianBias folder), the mean of its runs, for each redshift; with ratio, each divided
    by the first group."""
    figures = []
    for z in zs:
        curves = []
        for group in groups:
            curve = mean_curve(sorted(glob.glob(os.path.join(ROOT, group))), z)
            if curve is None:
                print(f"No finished runs in {group} at z = {z}")
            else:
                curves.append({**curve, "label": group + curve.pop("info")})
        if ratio and curves:
            reference = curves.pop(0)
            for curve in curves:
                P_ref = np.array([np.interp(curve["k"], reference["k"], P) for P in reference["P"]])
                curve["P"] = curve["P"] / P_ref
                curve["error"] = None if curve["error"] is None else curve["error"] / np.abs(P_ref)
        title = f"z = {z}" + (f", relative to {reference['label']}" if ratio and curves else "")
        figures.append({"title": title, "curves": curves, "ratio": ratio})
    return figures
