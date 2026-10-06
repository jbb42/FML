"""Shared by the plot scripts: load the 15 spectra of runs, and draw them as one figure with a panel per spectrum.

A run is a folder with snapshot_*_z<redshift>/pofk_ij.txt files, as written by the run scripts. Its 15 spectra P_ij(k),
i <= j, are of the fields 0 = matter, 1 = delta_L, 2 = delta_L^2, 3 = s^2 and 4 = nabla^2 delta_L.
"""
import functools
import glob
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

plt.rcParams["figure.max_open_warning"] = 0  # The pages of a PDF are all kept open until it is saved
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # The LagrangianBias folder
PAIRS = [(i, j) for i in range(5) for j in range(i, 5)]  # The 15 spectra, in the order of the arrays below
NAMES = [r"\delta_m", r"\delta_L", r"\delta_L^2", r"s^2", r"\nabla^2\delta_L"]


@functools.lru_cache(maxsize=None)  # Each run is read only once
def load_run(run_dir, z=0.0):
    """k and the 15 spectra of a run as an array [15, len(k)], or (None, None) if it has no output at redshift z."""
    snapshots = glob.glob(os.path.join(run_dir, f"snapshot_*_z{z:.3f}"))
    if not snapshots or not os.path.exists(os.path.join(snapshots[0], "pofk_44.txt")):
        return None, None
    spectra = [np.loadtxt(os.path.join(snapshots[0], f"pofk_{i}{j}.txt")) for i, j in PAIRS]
    return spectra[0][:, 0], np.array([s[:, 1] for s in spectra])


def load_runs(pattern, z=0.0):
    """k and the spectra of every run whose folder matches the glob pattern (relative to the LagrangianBias folder),
    as an array [run, 15, len(k)], or (None, None) if there are none."""
    runs = [load_run(run_dir, z) for run_dir in sorted(glob.glob(os.path.join(ROOT, pattern)))]
    runs = [(k, P) for k, P in runs if k is not None]
    if not runs:
        return None, None
    return runs[0][0], np.array([P for k, P in runs])


def plot_spectra(curves, title="", ratio=False, kmin=0.0, kmax=np.inf):
    """A figure with one panel per spectrum, the upper triangle of a 5 x 5 grid, and one line per curve.

    curves: list of dicts with "label", "k", "P" (array [15, len(k)]), optionally "error" (same shape, drawn as a band
    around P) and any other matplotlib line options, such as "color" or "ls". The spectra are drawn as |P| on log
    axes, or with ratio=True on a linear axis with a line at 1, which is the right choice for ratios of spectra.
    """
    fig, axes = plt.subplots(5, 5, figsize=(16, 10), sharex=True)
    fig.subplots_adjust(left=0.05, right=0.99, bottom=0.07, top=0.93, wspace=0.3, hspace=0.35)
    for ax in axes.flat:
        ax.axis("off")  # Panels below the diagonal stay empty

    for n, curve in enumerate(curves):
        options = {key: value for key, value in curve.items() if key not in ("k", "P", "error", "info")}
        options.setdefault("color", f"C{n % 10}")
        k = curve["k"]
        shown = (k >= kmin) & (k <= kmax)
        for p, (i, j) in enumerate(PAIRS):
            P = curve["P"][p] if ratio else np.abs(curve["P"][p])
            axes[i, j].plot(k[shown], P[shown], lw=1, **options)
            if curve.get("error") is not None:
                low, high = P - curve["error"][p], P + curve["error"][p]
                axes[i, j].fill_between(k[shown], low[shown], high[shown], color=options["color"], alpha=0.2, lw=0)

    for i, j in PAIRS:
        ax = axes[i, j]
        ax.axis("on")
        ax.set_xscale("log")
        ax.set_title(rf"${NAMES[i]} \times {NAMES[j]}$", fontsize=10, pad=3)
        ax.tick_params(labelbottom=(i == j), labelsize=8)  # x labels only on the bottom panel of each column
        if ratio:
            ax.axhline(1, color="k", lw=0.8)
            set_ratio_range(ax)
        else:
            ax.set_yscale("log")

    # Legend in the empty lower-left triangle
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="lower left", bbox_to_anchor=(0.05, 0.07), fontsize=10)
    fig.suptitle(title, y=0.99)
    fig.supxlabel(r"$k\ [h/\mathrm{Mpc}]$", y=0.01)
    fig.supylabel("ratio" if ratio else r"$|P_{ij}(k)|\ [(\mathrm{Mpc}/h)^3]$", x=0.01)
    return fig


def pages_varying_one_parameter(curves, names, max_lines=5):
    """Split curves for different parameter values into pages where only one parameter varies, with at most max_lines
    lines each, so they stay readable.

    curves: {parameter values (tuple): curve}, with the parameters named by names. Returns a list of (title, curves),
    where each curve is labelled by the value of the parameter that varies on its page (plus its "info", if any).
    Curves that differ from all others in more than one parameter are put on extra pages at the end.
    """
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
                pages.append((title, [{**curves[v], "label": f"{name} = {v[i]}" + curves[v].get("info", "")} for v in chunk]))
                shown.update(chunk)
    rest = [v for v in sorted(curves) if v not in shown]
    for start in range(0, len(rest), max_lines):
        chunk = rest[start:start + max_lines]
        pages.append(("Other", [{**curves[v], "label": ", ".join(f"{n} = {x}" for n, x in zip(names, v))
                                 + curves[v].get("info", "")} for v in chunk]))
    return pages


def set_ratio_range(ax):
    """y range from the 2-98th percentile of the lines (always including 1), so a single spike, e.g. where a
    spectrum crosses zero, does not flatten the panel."""
    y = np.concatenate([line.get_ydata() for line in ax.get_lines()[:-1]])  # The last line is the one at 1
    y = y[np.isfinite(y)]
    if y.size:
        low, high = np.percentile(y, [2, 98])
        low, high = min(low, 1.0), max(high, 1.0)
        ax.set_ylim(low - 0.05 * (high - low), high + 0.05 * (high - low))


def save(figures, filename):
    """Save the figures as the pages of one PDF, figures/<filename>."""
    path = os.path.join(ROOT, "figures", filename)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with PdfPages(path) as pdf:
        for fig in figures:
            pdf.savefig(fig)
            plt.close(fig)
    print(f"Saved {path} ({len(figures)} pages)")
