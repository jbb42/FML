#!/usr/bin/env python3
"""The fixes for the noisy low-k and nabla^2 spectra that need no new simulations, applied to the fiducial cosmology of
a scan (default planck) and compared with the spectra as measured, all averaged over the two phases:
1. GR spectra: the nabla^2 delta_L spectra as measured vs <X, nabla^2 delta_L> = -k^2 <X, 1> (and k^4 <1, 1> for the
   auto spectrum), as in Kokron et al. 2021 (arXiv:2101.11014), since the measured field is dominated by grid scales.
2. f(R)/GR boosts as measured vs with the -k^2 rule and smoothed with a Savitzky-Golay filter (order 3, 11 points),
   as Kokron et al. 2021 do for their N-body/LPT ratios.
3. (P_f(R) - P_GR) / sqrt(P_ii P_jj) (GR), measured and smoothed: unlike the ratio, this stays finite where a cross
   spectrum crosses zero, and shows how much each change matters relative to the auto spectra.
Writes figures/scans/<scan>_fixes.pdf and figures/scans/<scan>_fixes/P_ij.pdf.

    python3 scripts/plot_fixes.py [scan]
"""
import os
import sys

import numpy as np

from plot_page import save_pages
from plot_single import save_single
from spectra import PAIRS, ROOT, mean_curve

SCAN = sys.argv[1] if len(sys.argv) > 1 else "planck"


def phase_mean(model):
    """k and the 15 spectra of the fiducial run in this gravity model, averaged over the phases."""
    curve = mean_curve([os.path.join(ROOT, "results", "scans", SCAN, f"fiducial_{model}_seed1001_{phase}")
                        for phase in ("normal", "reversed")])
    return curve["k"], curve["P"]


def with_laplacian_rule(k, P):
    """The spectra with <X, nabla^2 delta_L> replaced by -k^2 <X, 1>, and <nabla^2, nabla^2> by k^4 <1, 1>."""
    P = P.copy()
    for p, (i, j) in enumerate(PAIRS):
        if j == 4:
            P[p] = k**4 * P[0] if i == 4 else -k**2 * P[PAIRS.index((0, i))]
    return P


def savgol(y, window=11, order=3):
    """Savitzky-Golay filter: each point replaced by a least-squares polynomial of this order through the window of
    points around it (at the ends, the first or last window)."""
    smooth = np.empty_like(y)
    for i in range(len(y)):
        start = min(max(i - window // 2, 0), len(y) - window)
        x = np.arange(start, start + window) - i
        smooth[i] = np.polyfit(x, y[start:start + window], order)[-1]  # The polynomial at x = 0
    return smooth


def normalised_difference(fr, gr):
    """(P_f(R) - P_GR) / sqrt(P_ii P_jj), with the GR auto spectra."""
    auto = [gr[PAIRS.index((i, i))] for i in range(5)]
    return np.array([(fr[p] - gr[p]) / np.sqrt(auto[i] * auto[j]) for p, (i, j) in enumerate(PAIRS)])


k, gr = phase_mean("GR")
_, fr = phase_mean("F5")
gr_fixed, fr_fixed = with_laplacian_rule(k, gr), with_laplacian_rule(k, fr)
smoothed = r"$-k^2$ rule + Savitzky-Golay"


def measured_and_fixed(measured, fixed, fixed_label=smoothed):
    return [{"label": "measured", "color": "C0", "k": k, "P": measured},
            {"label": fixed_label, "color": "C3", "k": k, "P": fixed}]


figures = [
    {"title": f"{SCAN} fiducial, GR, z = 0: measured vs " + r"$\langle X, \nabla^2\delta_L\rangle = -k^2\langle X, 1\rangle$",
     "curves": measured_and_fixed(gr, gr_fixed, r"$-k^2$ rule for $\nabla^2\delta_L$")},
    {"title": f"{SCAN} fiducial, z = 0: F5 / GR", "ratio": True,
     "curves": measured_and_fixed(fr / gr, np.array([savgol(b) for b in fr_fixed / gr_fixed]))},
    {"title": f"{SCAN} fiducial, z = 0: " + r"$(P^{F5}_{ij} - P^{GR}_{ij}) / \sqrt{P^{GR}_{ii} P^{GR}_{jj}}$",
     "ratio": True, "line": 0.0, "ylabel": "normalised difference",
     "curves": measured_and_fixed(normalised_difference(fr, gr),
                                  np.array([savgol(d) for d in normalised_difference(fr_fixed, gr_fixed)]))},
]
save_pages(figures, f"scans/{SCAN}_fixes.pdf")
save_single(figures, f"scans/{SCAN}_fixes")
