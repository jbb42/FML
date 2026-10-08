"""One figure per spectrum, with every curve of a figure in it (more than fit on the 15-panel pages), as
figures/<folder>/P_ij.pdf, with one page per figure.

From a run script:
    save_single(figures, "scans/planck")    # figures from spectra.py, written to figures/scans/planck/P_00.pdf, ...

From the command line, as plot_page.py, but with a folder instead of a PDF:
    python3 scripts/plot_single.py OUTPUT_FOLDER GROUP [GROUP ...] [--ratio] [--z Z ...]
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np

from plot_page import OPTIONS, save, set_ratio_range
from spectra import NAMES, PAIRS, group_figures, label

STYLE = {"font.family": "serif", "mathtext.fontset": "cm", "font.size": 12, "axes.labelsize": 14}
LINESTYLES = ["-", "--", ":", "-."]  # After every 10 colours


def plot_single(curves, p, title="", ratio=False, kmin=0.0, kmax=np.inf, line=1.0, ylabel="ratio"):
    """Spectrum number p (of spectra.PAIRS) of every curve, as |P| on log axes, or with ratio=True on a linear axis
    with a horizontal line at line."""
    i, j = PAIRS[p]
    fig, ax = plt.subplots(figsize=(10, 7))
    for n, curve in enumerate(curves):
        style = {"color": f"C{n % 10}", "ls": LINESTYLES[n // 10 % 4],
                 **{key: curve[key] for key in ("color", "ls") if key in curve}}
        k = curve["k"]
        shown = (k >= kmin) & (k <= kmax)
        P = curve["P"][p] if ratio else np.abs(curve["P"][p])
        ax.plot(k[shown], P[shown], lw=1, label=curve["label"], **style)
        if curve.get("error") is not None:
            error = curve["error"][p]
            ax.fill_between(k[shown], (P - error)[shown], (P + error)[shown], color=style["color"], alpha=0.05, lw=0)
    ax.set_xscale("log")
    if ratio:
        ax.axhline(line, color="k", lw=0.8)
        set_ratio_range(ax, line)
    else:
        ax.set_yscale("log")
    ax.set_title(rf"$P_{{{i}{j}}}$: ${NAMES[i]} \times {NAMES[j]}$" + (f"\n{title}" if title else ""))
    ax.set_xlabel(r"Wavenumber $k\ [h/\mathrm{Mpc}]$")
    ax.set_ylabel(ylabel if ratio else rf"$|P_{{{i}{j}}}(k)|\ [(\mathrm{{Mpc}}/h)^3]$")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    ax.grid(True, which="both", ls=":", alpha=0.4)
    fig.tight_layout()
    return fig


def save_single(figures, folder):
    """figures/<folder>/P_ij.pdf for each spectrum, with a page per figure (from spectra.py) holding all its curves."""
    pages = []
    for figure in figures:
        curves = figure["curves"]
        if isinstance(curves, dict):  # Label by the values that differ from the fiducial
            curves = [{**curve, "label": label(values, figure["names"], figure.get("fiducial")) + curve.get("info", "")}
                      for values, curve in curves.items()]
        if curves:
            pages.append((figure["title"], curves + figure.get("extra", []),
                          {key: figure[key] for key in OPTIONS if key in figure}))
    with plt.rc_context(STYLE):
        for p, (i, j) in enumerate(PAIRS):
            save([plot_single(curves, p, title, **options) for title, curves, options in pages],
                 f"{folder}/P_{i}{j}.pdf", verbose=False)
    if pages:
        print(f"Saved figures/{folder}/P_ij.pdf (15 files, {len(pages)} pages each)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output")
    parser.add_argument("groups", nargs="+")
    parser.add_argument("--ratio", action="store_true")
    parser.add_argument("--z", type=float, nargs="+", default=[0.0])
    args = parser.parse_args()
    save_single(group_figures(args.groups, args.z, args.ratio), args.output)
