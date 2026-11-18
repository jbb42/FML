#!/usr/bin/env python3
"""How much does each cosmological parameter matter for the mu - GR change of the spectra, for each gravity model?
From the full-scale cosmology test (results/mu_cosmology, run_all.py and run_extra.py): every parameter at -5 and +5
sigma (wa: +3), joint extremes, Omega_b, and three seeds of the fiducial cosmology for the noise.

D_ij = (P_mu - P_GR) / sqrt(P_ii^GR P_jj^GR) (GR and mu from the same seed in each cosmology, no control variate).
The effect of a parameter: |D(parameter at +-5 sigma) - D(fiducial)|, the larger of -5 and +5, averaged over a k range,
in % of max_k |D(fiducial)|, the worst of the main spectra (11, dd, d2d2, s2s2). The noise: the same for one seed
against another, divided by sqrt(2) (the noise of one run).

Writes figures/mu_boost/cosmology.pdf:
  page 1-2  overview (z = 0, z = 1): parameters x gravity models, colour = effect / noise (relevance), number = effect
            in %; large scales (k < 0.3) and small scales (0.6 < k < 1)
  page 3-4  k dependence (z = 0, z = 1): rows = parameters, columns = main spectra; signed effect
            (D_hi - D_lo) / 2 in % of max|D_fid| per model (colour: n, solid: A > 0, dashed: A < 0), grey: +-noise
    python3 scripts/plot_mu_cosmology.py
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

from runs import ROOT  # noqa: E402

sys.path.insert(0, os.path.join(ROOT, "results", "mu_cosmology"))
import analyze as A  # noqa: E402

MAIN = [(0, 0), (1, 1), (2, 2), (3, 3)]
NAMES = {(0, 0): r"$\delta_m\delta_m$", (1, 1): r"$\delta_L\delta_L$", (2, 2): r"$\delta_L^2\delta_L^2$",
         (3, 3): r"$s^2 s^2$"}
# One at a time (lo, hi), then joint extremes (single runs)
PARAMS = [(r"$\Omega_b$", "Ob_lo", "Ob_hi"), (r"$n_s$", "ns_lo", "ns_hi"), (r"$A_s$", "As_lo", "As_hi"),
          (r"$h$", "h_lo", "h_hi"), (r"$\Omega_{cdm}$", "OmegaCDM_lo", "OmegaCDM_hi"),
          (r"$w_a$", "wa_lo", "wa_hi"), (r"$w_0$", "w0_lo", "w0_hi")]
JOINT = [(r"$w_0$-, $w_a$-", "w0lo_walo"), (r"$w_0$+, $w_a$-", "w0hi_walo"), (r"$w_0$-, $w_a$+", "w0lo_wahi"),
         ("all +", "all_hi"), ("all -", "all_lo"), ("mixed a", "mix_a"), ("mixed b", "mix_b")]
RANGES = [("large scales, 0.01 < k < 0.3", 0.01, 0.3), ("small scales, 0.6 < k < 1", 0.6, 1.0)]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e0"


def models():
    """The mu models ordered by n, then A, with their parameters."""
    out = []
    for m in A.run_all.MODELS:
        if m == "GR":
            continue
        a, n = (float(x[1:]) for x in m.split("_"))
        out.append((n, a, m))
    return [(m, a, n) for n, a, m in sorted(out)]


def effect(z, label, model, lo, hi, which="max"):
    """Mean |D_label - D_fid| over lo < k <= hi in % of max|D_fid|, the worst main spectrum."""
    A.Z = z
    k, Df = A.D("fiducial", model)
    _, Dc = A.D(label, model)
    m = (k > lo) & (k <= hi)
    vals = []
    for p in MAIN:
        i = A.PAIRS.index(p)
        scale = np.max(np.abs(Df[i][(k > 0.01) & (k < 1.0)]))
        vals.append(100 * np.mean(np.abs(Dc[i][m] - Df[i][m])) / scale)
    return max(vals)


def noise(z, model, lo, hi):
    return np.mean([effect(z, s, model, lo, hi) for s in ("fiducial_s1002", "fiducial_s1003")]) / np.sqrt(2)


def model_label(m, a, n):
    return f"A = {a:g}\nn = {n:g}"


def overview(pdf, z):
    ms = models()
    rows = [(name, [lo_l, hi_l]) for name, lo_l, hi_l in PARAMS] + [(name, [l]) for name, l in JOINT]
    fig, axes = plt.subplots(1, 2, figsize=(15, 8.2), sharey=True)
    fig.subplots_adjust(left=0.1, right=0.9, top=0.84, bottom=0.12, wspace=0.06)
    vmax = 10.0
    for ax, (title, lo, hi) in zip(axes, RANGES):
        N = np.array([noise(z, m, lo, hi) for m, _, _ in ms])
        E = np.array([[max(effect(z, l, m, lo, hi) for l in labels) for m, _, _ in ms] for _, labels in rows])
        R = E / N
        im = ax.imshow(np.clip(R, 0.1, vmax), cmap="Blues", norm=matplotlib.colors.LogNorm(0.3, vmax), aspect="auto")
        for r in range(len(rows)):
            for c in range(len(ms)):
                ax.text(c, r, f"{E[r, c]:.1f}", ha="center", va="center", fontsize=8.5,
                        color="white" if R[r, c] > 3 else INK)
        ax.axhline(len(PARAMS) - 0.5, color="white", lw=4)
        ax.set_xticks(range(len(ms)))
        ax.set_xticklabels([model_label(*x) for x in ms], fontsize=8.5, color=INK)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([r[0] for r in rows], fontsize=10, color=INK)
        ax.set_title(f"{title}\nnoise (one seed): " + " ".join(f"{v:.1f}" for v in N), fontsize=10, color=INK, loc="left")
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        # Separate the models by n
        ns = [n for _, _, n in ms]
        for c in range(1, len(ms)):
            if ns[c] != ns[c - 1]:
                ax.axvline(c - 0.5, color="white", lw=3)
    cax = fig.add_axes([0.915, 0.12, 0.012, 0.72])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("effect / noise of one run (relevance)", color=MUTED)
    cb.outline.set_visible(False)
    fig.suptitle(f"How much does each parameter change the mu - GR change D of the spectra?  z = {z:g}", x=0.1,
                 ha="left", fontsize=14, color=INK)
    fig.text(0.1, 0.885, "Number: |D(cosmology) - D(fiducial)| in % of max|D(fiducial)|, worst of the main spectra "
             "(11, dd, d2d2, s2s2) and of -5 / +5 sigma (wa: +3). Upper block: one parameter at a time; lower block: "
             "joint extremes.\nColour: the number divided by the realisation noise of one run (from seeds 1001-1003); "
             "pale (< 1): not distinguishable from noise. Models ordered by n.", fontsize=9, color=MUTED)
    pdf.savefig(fig)
    plt.close(fig)


def kdependence(pdf, z):
    ms = models()
    nvals = sorted({n for _, _, n in ms})
    blues = ["#86b6ef", "#4a90e2", "#2a78d6", "#0d3e80"]  # Sequential by n (light: small n)
    colour = {n: blues[int(round(i * (len(blues) - 1) / max(1, len(nvals) - 1)))] for i, n in enumerate(nvals)}
    fig, axes = plt.subplots(len(PARAMS), len(MAIN), figsize=(15, 17), sharex=True)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.94, bottom=0.06, hspace=0.25, wspace=0.22)
    A.Z = z
    for r, (name, lo_l, hi_l) in enumerate(PARAMS):
        for c, p in enumerate(MAIN):
            ax = axes[r, c]
            i = A.PAIRS.index(p)
            noise_curves = []
            for m, a, n in ms:
                k, Df = A.D("fiducial", m)
                scale = np.max(np.abs(Df[i][(k > 0.01) & (k < 1.0)]))
                _, Dlo = A.D(lo_l, m)
                _, Dhi = A.D(hi_l, m)
                ax.plot(k, 100 * (Dhi[i] - Dlo[i]) / 2 / scale, color=colour[n], lw=1.3, ls="-" if a > 0 else "--",
                        label=f"A = {a:g}, n = {n:g}")
                seeds = np.array([A.D(s, m)[1][i] for s in ("fiducial", "fiducial_s1002", "fiducial_s1003")])
                noise_curves.append(100 * np.std(seeds, axis=0, ddof=1) / scale)
            nz = np.max(noise_curves, axis=0)
            nz = np.convolve(nz, np.ones(5) / 5, mode="same")  # Smooth over a few bins
            ax.fill_between(k, -nz, nz, color=GRID, lw=0, zorder=0)
            ax.axhline(0, color=MUTED, lw=0.6)
            ax.set_xscale("log")
            ax.set_xlim(0.01, 1.0)
            ax.grid(True, color=GRID, lw=0.5)
            ax.tick_params(labelsize=8, colors=MUTED)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            if r == 0:
                ax.set_title(NAMES[p], fontsize=12, color=INK)
            if c == 0:
                ax.set_ylabel(name + "\n(% of max|D|)", fontsize=10, color=INK)
            if r == len(PARAMS) - 1:
                ax.set_xlabel(r"$k$ [$h$/Mpc]", color=INK)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=7, fontsize=9, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle(f"k dependence, z = {z:g}: (D(+5 sigma) - D(-5 sigma)) / 2 in % of max|D(fiducial)|; "
                 "grey: +-noise of one run (largest over models)", x=0.07, ha="left", fontsize=12, color=INK, y=0.975)
    pdf.savefig(fig)
    plt.close(fig)


if __name__ == "__main__":
    plt.rcParams.update({"figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "font.size": 10})
    path = os.path.join(ROOT, "figures", "mu_boost", "cosmology.pdf")
    with PdfPages(path) as pdf:
        for z in (0.0, 1.0):
            overview(pdf, z)
        for z in (0.0, 1.0):
            kdependence(pdf, z)
    print("Saved", path)
