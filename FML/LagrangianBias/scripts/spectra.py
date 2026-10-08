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
    "ratio", "kmin", "kmax", "line", "ylabel": options of plot_page() and plot_single()
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


def boosts(runs, model, z=0.0):
    """{values: f(R)/GR curve} for the runs (runs.Run) in this model, each divided by the GR run with the same values,
    seed and phase. Per seed, both are averaged over the phases that finished in both, which cancels the leading
    cosmic variance, also of spectra that are odd in the initial field. The curve is the mean over the seeds, with
    their standard error; where f(R) and GR have opposite signs the ratio is meaningless and left out."""
    gr = {(r.values, r.seed, r.reverse): r.folder for r in runs if r.model == "GR"}
    pairs = defaultdict(lambda: defaultdict(list))  # values -> seed -> [(GR, f(R)) spectra]
    for r in runs:
        if r.model == model:
            P_gr, P_fr = load_run(gr[(r.values, r.seed, r.reverse)], z), load_run(r.folder, z)
            if P_gr[0] is not None and P_fr[0] is not None:
                pairs[r.values][r.seed].append((P_gr, P_fr))
    curves = {}
    for values, seeds in pairs.items():
        per_seed = []
        for runs_of_seed in seeds.values():
            P_gr = np.mean([g[1] for g, f in runs_of_seed], axis=0)
            P_fr = np.mean([f[1] for g, f in runs_of_seed], axis=0)
            per_seed.append(np.where(np.sign(P_fr) == np.sign(P_gr), P_fr / P_gr, np.nan))
        with warnings.catch_warnings():  # k where every seed is NaN give NaN
            warnings.simplefilter("ignore", RuntimeWarning)
            n = len(per_seed)
            curves[values] = {"k": runs_of_seed[0][0][0], "P": np.nanmean(per_seed, axis=0),
                              "error": np.nanstd(per_seed, axis=0, ddof=1) / np.sqrt(n) if n > 1 else None,
                              "info": f" ({n} seed{'s' if n > 1 else ''})"}
    return curves


def scan_figures(scan):
    """The figures of a scan (runs.scan), for each redshift: the spectra (mean over seeds and phases, in GR if gravity
    models are compared), the same divided by the fiducial run (or else the first), and the f(R)/GR boosts."""
    figures = []
    for z in redshifts([r.folder for r in scan.runs]):
        title = f"{scan.name}, z = {z}"
        names = {"names": scan.parameters, "fiducial": scan.reference}
        spectra = {}
        for values in dict.fromkeys(r.values for r in scan.runs):
            curve = mean_curve([r.folder for r in scan.runs if r.values == values and r.model == scan.models[0]], z)
            if curve:
                spectra[values] = curve
        figures.append({"title": title, "curves": spectra, **names})
        if len(spectra) > 1:
            reference = scan.reference if scan.reference in spectra else min(spectra)
            ref = spectra[reference]
            ratios = {values: {**c, "error": None, "P": c["P"] / [np.interp(c["k"], ref["k"], P) for P in ref["P"]]}
                      for values, c in spectra.items()}
            relative_to = "the fiducial" if reference == scan.reference else \
                ", ".join(f"{n} = {v}" for n, v in zip(scan.parameters, reference))
            figures.append({"title": f"{title}, relative to {relative_to}", "curves": ratios, "ratio": True, **names})
        for model in scan.models[1:]:
            figures.append({"title": f"{title}: {model} / GR", "curves": boosts(scan.runs, model, z), "ratio": True,
                            **names})
    return figures


def pages_varying_one_parameter(curves, names, max_lines=5):
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
