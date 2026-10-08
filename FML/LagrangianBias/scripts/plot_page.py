"""All 15 spectra on one page, a panel per spectrum, with at most 5 lines per page: figures with more curves are split
into pages where one parameter varies.

From a run script:
    save_pages(figures, "scans/planck.pdf")    # figures from spectra.py, written to figures/scans/planck.pdf

From the command line, a line per group of runs (a run folder, or a quoted glob pattern matching several, averaged):
    python3 scripts/plot_page.py OUTPUT GROUP [GROUP ...] [--ratio] [--z Z ...]
    e.g. plot_page.py ensemble/all_spectra.pdf 'results/ensemble_100seeds/*'
         plot_page.py scans/mesh.pdf results/scans/planck/fiducial_GR_seed1001_normal results/scans/planck_mesh/* --ratio
    --ratio divides every group by the first, --z gives the redshifts (default 0).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from spectra import NAMES, PAIRS, ROOT, group_figures, pages_varying_one_parameter

plt.rcParams["figure.max_open_warning"] = 0  # The pages of a PDF are all kept open until it is saved
OPTIONS = ("ratio", "kmin", "kmax", "line", "ylabel")  # Of a figure, passed on to the plot functions


def plot_page(curves, title="", ratio=False, kmin=0.0, kmax=np.inf, line=1.0, ylabel="ratio"):
    """A page with one panel per spectrum, the upper triangle of a 5 x 5 grid, and one line per curve. The spectra
    are drawn as |P| on log axes, or with ratio=True on a linear axis with a horizontal line at line (1 for ratios,
    0 for differences)."""
    fig, axes = plt.subplots(5, 5, figsize=(16, 10), sharex=True)
    fig.subplots_adjust(left=0.05, right=0.99, bottom=0.07, top=0.93, wspace=0.3, hspace=0.35)
    for ax in axes.flat:
        ax.axis("off")  # Panels below the diagonal stay empty

    for n, curve in enumerate(curves):
        style = {"color": f"C{n % 10}", **{key: curve[key] for key in ("color", "ls") if key in curve}}
        k = curve["k"]
        shown = (k >= kmin) & (k <= kmax)
        for p, (i, j) in enumerate(PAIRS):
            P = curve["P"][p] if ratio else np.abs(curve["P"][p])
            axes[i, j].plot(k[shown], P[shown], lw=1, label=curve["label"], **style)
            if curve.get("error") is not None:
                low, high = P - curve["error"][p], P + curve["error"][p]
                axes[i, j].fill_between(k[shown], low[shown], high[shown], color=style["color"], alpha=0.2, lw=0)

    for i, j in PAIRS:
        ax = axes[i, j]
        ax.axis("on")
        ax.set_xscale("log")
        ax.set_title(rf"${NAMES[i]} \times {NAMES[j]}$", fontsize=10, pad=3)
        ax.tick_params(labelbottom=(i == j), labelsize=8)  # x labels only on the bottom panel of each column
        if ratio:
            ax.axhline(line, color="k", lw=0.8)
            set_ratio_range(ax, line)
        else:
            ax.set_yscale("log")

    # Legend in the empty lower-left triangle
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="lower left", bbox_to_anchor=(0.05, 0.07), fontsize=10)
    fig.suptitle(title, y=0.99)
    fig.supxlabel(r"$k\ [h/\mathrm{Mpc}]$", y=0.01)
    fig.supylabel(ylabel if ratio else r"$|P_{ij}(k)|\ [(\mathrm{Mpc}/h)^3]$", x=0.01)
    return fig


def set_ratio_range(ax, line=1.0):
    """y range from the 2-98th percentile of the lines (always including the horizontal line), so a single spike,
    e.g. where a spectrum crosses zero, does not flatten the panel."""
    y = np.concatenate([l.get_ydata() for l in ax.get_lines()[:-1]])  # The last line is the horizontal one
    y = y[np.isfinite(y)]
    if y.size:
        low, high = np.percentile(y, [2, 98])
        low, high = min(low, line), max(high, line)
        ax.set_ylim(low - 0.05 * (high - low), high + 0.05 * (high - low))


def save(figures, filename, verbose=True):
    """Save matplotlib figures as the pages of one PDF, figures/<filename> (nothing if there are none)."""
    if not figures:
        return
    path = os.path.join(ROOT, "figures", filename)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with PdfPages(path) as pdf:
        for fig in figures:
            pdf.savefig(fig)
            plt.close(fig)
    if verbose:
        print(f"Saved {path} ({len(figures)} pages)")


def save_pages(figures, filename):
    """figures/<filename>: the figures (from spectra.py), each split into pages of at most 5 lines, where only one
    parameter varies."""
    pages = []
    for figure in figures:
        curves = figure["curves"]
        if isinstance(curves, dict):
            split = pages_varying_one_parameter(curves, figure["names"])
        else:
            split = [("", curves[start:start + 5]) for start in range(0, len(curves), 5)]
        for subtitle, page in split:
            title = ": ".join(text for text in (figure["title"], subtitle) if text)
            pages.append(plot_page(page + figure.get("extra", []), title,
                                   **{key: figure[key] for key in OPTIONS if key in figure}))
    save(pages, filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output")
    parser.add_argument("groups", nargs="+")
    parser.add_argument("--ratio", action="store_true")
    parser.add_argument("--z", type=float, nargs="+", default=[0.0])
    args = parser.parse_args()
    save_pages(group_figures(args.groups, args.z, args.ratio), args.output)
