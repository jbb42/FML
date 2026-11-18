#!/usr/bin/env python3
"""How much does each cosmological parameter matter for the mu - GR change of the spectra, for each gravity model?
From the full-scale cosmology test (results/mu_cosmology, run_all.py and run_extra.py): every parameter at -5 and +5
sigma (wa: +3), joint extremes, Omega_b, and three seeds of the fiducial cosmology for the noise.

D_ij = (P_mu - P_GR) / sqrt(P_ii^GR P_jj^GR) (GR and mu from the same seed in each cosmology, no control variate).
The effect of a parameter: |D(parameter at +-5 sigma) - D(fiducial)|, the larger of -5 and +5, averaged over a k range,
in % of max_k |D(fiducial)|, the worst of the main spectra (11, dd, d2d2, s2s2). The noise: the same for one seed
against another, divided by sqrt(2) (the noise of one run).

Writes figures/mu_boost/cosmology.pdf:
  page 1    verdict (z = 0): which parameters can be fixed (1-sigma shift against the noise); the pages below
            are then numbered 2-7
  page 1-2  overview (z = 0, z = 1): per parameter one dot per gravity model (colour: n), the noise of one run as a
            grey band; large (k < 0.3), intermediate (0.3 < k < 0.6) and small scales (0.6 < k < 1)
  page 5-6  all 15 spectra (z = 0, z = 1) in the 15-spectra layout: per parameter the change of D at large,
            intermediate and small scales (worst gravity model) and the noise of one run
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
RANGES = [("large scales, 0.01 < k < 0.3", 0.01, 0.3), ("intermediate scales, 0.3 < k < 0.6", 0.3, 0.6),
          ("small scales, 0.6 < k < 1", 0.6, 1.0)]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
# One colour rule for every heatmap: the effect in %, in four steps of one hue
CLASS_EDGES = [0.0, 1.0, 3.0, 10.0, 1e9]
CLASS_COLOURS = ["#eaf2fc", "#a9cdf3", "#4a90e2", "#0d3e80"]
CLASS_LABELS = ["< 1 %", "1 - 3 %", "3 - 10 %", "> 10 %"]
CMAP = matplotlib.colors.ListedColormap(CLASS_COLOURS)
NORM = matplotlib.colors.BoundaryNorm(CLASS_EDGES, CMAP.N)


def draw_cells(ax, E, R, fontsize):
    """Colour and number: the effect E in %; (number) in brackets: not larger than the noise of one run (R = E / noise < 1)."""
    ax.imshow(E, cmap=CMAP, norm=NORM, aspect="auto")
    for r in range(E.shape[0]):
        for c in range(E.shape[1]):
            txt = f"{E[r, c]:.1f}" if E[r, c] < 99.5 else f"{E[r, c]:.0f}"
            ax.text(c, r, f"({txt})" if R[r, c] < 1 else txt, ha="center", va="center", fontsize=fontsize,
                    color="white" if E[r, c] >= 3 else INK)


def class_legend(fig, y):
    handles = [matplotlib.patches.Patch(color=c, label=l) for c, l in zip(CLASS_COLOURS, CLASS_LABELS)]
    handles.append(matplotlib.patches.Patch(color="none", label="(x): within the noise of one run"))
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.99, y), ncol=5, frameon=False, fontsize=10,
               title="change of D in % of its size", title_fontsize=10)


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
    fig, axes = plt.subplots(1, len(RANGES), figsize=(21, 8.2), sharey=True)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.8, bottom=0.12, wspace=0.06)
    for ax, (title, lo, hi) in zip(axes, RANGES):
        N = np.array([noise(z, m, lo, hi) for m, _, _ in ms])
        E = np.array([[max(effect(z, l, m, lo, hi) for l in labels) for m, _, _ in ms] for _, labels in rows])
        R = E / N
        draw_cells(ax, E, R, 8.5)
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
    class_legend(fig, 0.995)
    fig.suptitle(f"How much does each parameter change the mu - GR change D of the spectra?  z = {z:g}", x=0.1,
                 ha="left", fontsize=14, color=INK)
    fig.text(0.1, 0.885, "Colour and number: |D(cosmology) - D(fiducial)| in % of max|D(fiducial)|, worst of the main "
             "spectra (11, dd, d2d2, s2s2) and of -5 / +5 sigma (wa: +3). Upper block: one parameter at a time; lower "
             "block: joint extremes.\n(Number) in brackets: not larger than the realisation noise of one run (seeds "
             "1001-1003). Models ordered by n.", fontsize=9, color=MUTED)
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


_DCACHE = {}


def D_cached(z, label, model):
    key = (z, label, model)
    if key not in _DCACHE:
        A.Z = z
        _DCACHE[key] = A.D(label, model)
    return _DCACHE[key]


OPNAMES = [r"1", r"\delta", r"\delta^2", r"s^2", r"\nabla^2\delta"]


def all_spectra(pdf, z):
    """All 15 spectra: rows = parameters and joint extremes, columns = spectra; per spectrum, the worst gravity model.
    Colour: effect / noise (both for that spectrum); number: effect in % of max|D_fid| of that spectrum."""
    ms = [m for m, _, _ in models()]
    rows = [(name, [lo_l, hi_l]) for name, lo_l, hi_l in PARAMS] + [(name, [l]) for name, l in JOINT]
    seeds = ("fiducial_s1002", "fiducial_s1003")
    fig, axes = plt.subplots(len(RANGES), 1, figsize=(16, 22))
    fig.subplots_adjust(left=0.13, right=0.98, top=0.93, bottom=0.05, hspace=0.42)
    for ax, (title, lo, hi) in zip(axes, RANGES):
        E = np.zeros((len(rows), 15))
        R = np.zeros((len(rows), 15))
        for p in range(15):
            for mi, m in enumerate(ms):
                k, Df = D_cached(z, "fiducial", m)
                sel = (k > lo) & (k <= hi)
                scale = np.max(np.abs(Df[p][(k > 0.01) & (k < 1.0)]))
                ch = lambda lab: 100 * np.mean(np.abs(D_cached(z, lab, m)[1][p][sel] - Df[p][sel])) / scale
                nz = np.mean([ch(s_) for s_ in seeds]) / np.sqrt(2)
                for r, (_, labels) in enumerate(rows):
                    e = max(ch(l) for l in labels)
                    E[r, p] = max(E[r, p], e)
                    R[r, p] = max(R[r, p], e / nz)
        draw_cells(ax, E, R, 7.5)
        ax.axhline(len(PARAMS) - 0.5, color="white", lw=4)
        for c in range(1, 15):
            if A.PAIRS[c][0] != A.PAIRS[c - 1][0]:
                ax.axvline(c - 0.5, color="white", lw=3)
        ax.set_xticks(range(15))
        ax.set_xticklabels([f"${OPNAMES[i]} \\times {OPNAMES[j]}$" for i, j in A.PAIRS], fontsize=8.5, color=INK,
                           rotation=35, ha="right")
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([r[0] for r in rows], fontsize=9, color=INK)
        ax.set_title(title, fontsize=11, color=INK, loc="left")
        ax.tick_params(length=0)
        for s_ in ax.spines.values():
            s_.set_visible(False)
    class_legend(fig, 0.999)
    fig.suptitle(f"All 15 spectra: how much each parameter changes the mu - GR change D, z = {z:g}", x=0.02, ha="left",
                 fontsize=14, color=INK, y=0.985)
    fig.text(0.02, 0.955, "Colour and number: |D(cosmology) - D(fiducial)| in % of max|D(fiducial)| of that spectrum, worst "
             "of -5 / +5 sigma and of the 7 gravity models.\n(Number) in brackets: the effect is not larger than the "
             "realisation noise of one run (seeds 1001-1003) for that spectrum, i.e. not distinguishable from noise.",
             fontsize=9, color=MUTED)
    pdf.savefig(fig)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------- dot plots
N_COLOURS = {1.0: "#86b6ef", 2.0: "#3987e5", 3.0: "#1c5cab", 4.0: "#0d366b"}  # Sequential blue by n (light: early MG)
SCALE_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a"]  # Large, intermediate, small scales (categorical slots 1-3)
NOISE_COLOUR = "#d9d8d2"


def _style(ax, xlim=(0.03, 100)):
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.grid(True, axis="x", color=GRID, lw=0.6)
    ax.tick_params(colors=MUTED, labelsize=8, length=0)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))


def overview_dots(pdf, z):
    """Main spectra: per parameter, one dot per gravity model (colour: n), the noise of one run as a grey band."""
    ms = models()
    rows = [(name, [lo_l, hi_l]) for name, lo_l, hi_l in PARAMS] + [(name, [l]) for name, l in JOINT]
    ypos = list(range(len(PARAMS))) + [len(PARAMS) + 0.8 + i for i in range(len(JOINT))]
    fig, axes = plt.subplots(1, len(RANGES), figsize=(18, 8.5), sharey=True)
    fig.subplots_adjust(left=0.13, right=0.99, top=0.83, bottom=0.09, wspace=0.08)
    for ax, (title, lo, hi) in zip(axes, RANGES):
        N = np.array([noise(z, m, lo, hi) for m, _, _ in ms])
        ax.axvspan(0.05, N.max(), color="#ecebe6", lw=0, zorder=0)
        ax.axvspan(N.min(), N.max(), color=NOISE_COLOUR, lw=0, zorder=0)
        ax.text(np.sqrt(0.05 * N.max()), -0.9, "drowned in noise", ha="center", va="center", fontsize=8.5,
                color=MUTED, style="italic")
        for (name, labels), y in zip(rows, ypos):
            E = np.array([max(effect(z, l, m, lo, hi) for l in labels) for m, _, _ in ms])
            ax.plot([E.min(), E.max()], [y, y], color=MUTED, lw=1.2, zorder=1)
            for (m, a, n), e in zip(ms, E):
                ax.plot(e, y, "o" if a > 0 else "s", ms=8, color=N_COLOURS[n], mec="white", mew=1.2, zorder=2)
            ax.text(min(E.max() * 1.3, 40), y, f"{E.max():.1f}%", va="center", fontsize=8.5, color=INK)
        ax.axhline(len(PARAMS) - 0.1, color=GRID, lw=1)
        _style(ax, (0.05, 60))
        ax.set_ylim(len(PARAMS) + 0.8 + len(JOINT) - 0.4, -1.4)
        ax.set_title(title, fontsize=11.5, color=INK, loc="left")
        ax.set_xlabel("change of D, % of its size (log scale)", color=MUTED, fontsize=9)
    axes[0].set_yticks(ypos)
    axes[0].set_yticklabels([r[0] for r in rows], fontsize=10, color=INK)
    handles = [matplotlib.lines.Line2D([], [], ls="", marker="o", ms=7, color=N_COLOURS[n], label=f"n = {n:g}")
               for n in sorted(N_COLOURS)]
    handles += [matplotlib.lines.Line2D([], [], ls="", marker="o", ms=7, color=MUTED, label="A > 0"),
                matplotlib.lines.Line2D([], [], ls="", marker="s", ms=7, color=MUTED, label="A < 0"),
                matplotlib.patches.Patch(color=NOISE_COLOUR, label="noise of one run (range over models)")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.13, 0.905), ncol=7, frameon=False, fontsize=9.5)
    fig.suptitle(f"Which parameters change the mu - GR change D, for which gravity models?  z = {z:g}", x=0.01,
                 ha="left", fontsize=14, color=INK, y=0.985)
    fig.text(0.01, 0.935, "One dot per gravity model (worst of the main spectra 11, dd, d2d2, s2s2 and of -5 / +5 sigma; "
             "wa: +3); the number is the largest. Inside or left of the grey band: drowned in the noise of one run. "
             "Upper block: one parameter at a time; lower block: joint extremes.", fontsize=9, color=MUTED)
    pdf.savefig(fig)
    plt.close(fig)


def all_spectra_dots(pdf, z):
    """All 15 spectra in the 15-spectra layout: per parameter the change of D (worst gravity model) at large,
    intermediate and small scales, and the noise of one run (grey tick, small scales)."""
    ms = [m for m, _, _ in models()]
    names = [r"1", r"\delta", r"\delta^2", r"s^2", r"\nabla^2\delta"]
    seeds = ("fiducial_s1002", "fiducial_s1003")
    fig, axes = plt.subplots(5, 5, figsize=(19, 13))
    fig.subplots_adjust(left=0.06, right=0.99, top=0.9, bottom=0.05, wspace=0.12, hspace=0.45)
    for ax in axes.flat:
        ax.axis("off")
    for p, (i, j) in enumerate(A.PAIRS):
        ax = axes[i, j]
        ax.axis("on")
        for r, (name, lo_l, hi_l) in enumerate(PARAMS):
            vals, nzs = [], []
            for _, lo, hi in RANGES:
                E, Nz = 0.0, 0.0
                for m in ms:
                    k, Df = D_cached(z, "fiducial", m)
                    sel = (k > lo) & (k <= hi)
                    scale = np.max(np.abs(Df[p][(k > 0.01) & (k < 1.0)]))
                    ch = lambda lab: 100 * np.mean(np.abs(D_cached(z, lab, m)[1][p][sel] - Df[p][sel])) / scale
                    E = max(E, ch(lo_l), ch(hi_l))
                    Nz = max(Nz, np.mean([ch(s_) for s_ in seeds]) / np.sqrt(2))
                vals.append(E)
                nzs.append(Nz)
            ax.plot([min(vals), max(vals)], [r, r], color=GRID, lw=1.2, zorder=1)
            ax.fill_betweenx([r - 0.5, r + 0.5], 0.03, nzs[2], color="#e3e2dc", lw=0, zorder=0)
            for c, v in zip(SCALE_COLOURS, vals):
                ax.plot(v, r, "o", ms=5, color=c, mec="white", mew=0.8, zorder=2)
        _style(ax, (0.03, 60))
        ax.set_ylim(len(PARAMS) - 0.5, -0.5)
        ax.set_yticks(range(len(PARAMS)))
        ax.set_yticklabels([n for n, _, _ in PARAMS] if j == i else [], fontsize=8, color=INK)
        ax.set_title(rf"${names[i]} \times {names[j]}$", fontsize=10, color=INK, pad=3)
    handles = [matplotlib.lines.Line2D([], [], ls="", marker="o", ms=7, color=c, label=t.split(",")[0] + ":" +
                                       t.split(",")[1]) for c, (t, _, _) in zip(SCALE_COLOURS, RANGES)]
    handles.append(matplotlib.patches.Patch(color="#e3e2dc", label="drowned in the noise of one run (small scales)"))
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.06, 0.08), ncol=1, frameon=False, fontsize=10)
    fig.suptitle(f"All 15 spectra: change of the mu - GR change D per parameter (+-5 sigma), z = {z:g}", x=0.01,
                 ha="left", fontsize=14, color=INK, y=0.985)
    fig.text(0.01, 0.94, "x: change of D in % of its size (log scale), worst of -5 / +5 sigma and of the 7 gravity models. "
             "Dots in the grey area are drowned in the noise of one run. Joint extremes: pages 1-2.",
             fontsize=9.5, color=MUTED)
    pdf.savefig(fig)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------- verdict
SIGMAS = {"Ob": (5, 5), "ns": (5, 5), "As": (5, 5), "h": (5, 5), "OmegaCDM": (5, 5), "wa": (5, 3), "w0": (5, 5)}
FIX_COLOUR, VARY_COLOUR = "#2a78d6", "#eb6834"


def effect_1sigma(z, lo_l, hi_l, model, lo, hi, pairs=None):
    """The change of D for a 1 sigma shift (linear in the +-5 sigma runs; wa high: +3 sigma), worst side and worst
    of the given spectra (default: the main ones), in % of max|D_fid| of each spectrum."""
    key = lo_l.rsplit("_", 1)[0]
    s_lo, s_hi = SIGMAS[key]
    k, Df = D_cached(z, "fiducial", model)
    sel = (k > lo) & (k <= hi)
    out = 0.0
    for p in pairs if pairs is not None else [A.PAIRS.index(q) for q in MAIN]:
        scale = np.max(np.abs(Df[p][(k > 0.01) & (k < 1.0)]))
        for lab, sg in ((lo_l, s_lo), (hi_l, s_hi)):
            out = max(out, 100 * np.mean(np.abs(D_cached(z, lab, model)[1][p][sel] - Df[p][sel])) / scale / sg)
    return out


def verdict_page(pdf, z):
    ms = models()
    ranges = [(t.split(",")[0], lo, hi) for t, lo, hi in RANGES]
    noise_max = [max(noise(z, m, lo, hi) for m, _, _ in ms) for _, lo, hi in ranges]
    rows = []
    for name, lo_l, hi_l in PARAMS:
        main = [max(effect_1sigma(z, lo_l, hi_l, m, lo, hi) for m, _, _ in ms) for _, lo, hi in ranges]
        fix = all(e <= nz for e, nz in zip(main, noise_max))
        # The worst of all 15 spectra (small scales), skipping spectra where the effect is within their own noise
        worst, worst_p = 0.0, 0
        for p in range(15):
            e = max(effect_1sigma(z, lo_l, hi_l, m, 0.6, 1.0, [p]) for m, _, _ in ms)
            if e > worst:
                worst, worst_p = e, p
        rows.append((name, main, fix, worst, worst_p))
    names = [r"1", r"\delta", r"\delta^2", r"s^2", r"\nabla^2\delta"]

    fig = plt.figure(figsize=(17, 8.5))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.9], left=0.08, right=0.99, top=0.8, bottom=0.1, wspace=0.08)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    for ax, (title, lo, hi), nz in zip(axes, ranges, noise_max):
        ax.axvspan(0.01, nz, color="#e3e2dc", lw=0, zorder=0)
        ax.text(np.sqrt(0.01 * nz), -0.85, "within the noise\nof one run", ha="center", va="center", fontsize=8.5,
                color=MUTED, style="italic")
        for y, (name, main, fix, _, _) in enumerate(rows):
            e = main[ranges.index((title, lo, hi))]
            col = FIX_COLOUR if fix else VARY_COLOUR
            ax.plot([0.01, e], [y, y], color=col, lw=2.5, alpha=0.35, solid_capstyle="butt", zorder=1)
            ax.plot(e, y, "o", ms=10, color=col, mec="white", mew=1.5, zorder=2)
            ax.text(e * 1.35, y, f"{e:.2f}%" if e < 1 else f"{e:.1f}%", va="center", fontsize=9, color=INK)
        _style(ax, (0.01, 30))
        ax.set_ylim(len(rows) - 0.5, -1.5)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([r[0] for r in rows] if ax is axes[0] else [], fontsize=12, color=INK)
        ax.set_title(f"{title}, {lo:g} < k < {hi:g}", fontsize=11.5, color=INK, loc="left")
        ax.set_xlabel("change of D for a 1-sigma shift, % of its size", color=MUTED, fontsize=9)
    vax = fig.add_subplot(gs[0, 3])
    vax.axis("off")
    vax.set_ylim(len(rows) - 0.5, -1.5)
    vax.set_xlim(0, 1)
    vax.text(0.02, -0.85, "verdict", fontsize=11.5, color=INK, va="center")
    for y, (name, main, fix, worst, wp) in enumerate(rows):
        i, j = A.PAIRS[wp]
        vax.text(0.02, y - 0.12, "FIX" if fix else "VARY", fontsize=13, fontweight="bold", va="center",
                 color=FIX_COLOUR if fix else VARY_COLOUR)
        vax.text(0.22, y - 0.12, "within the noise at all scales" if fix else "above the noise", fontsize=9.5,
                 va="center", color=INK)
        vax.text(0.22, y + 0.22, rf"worst of all 15: ${names[i]}\times{names[j]}$ {worst:.1f}% (k > 0.6)", fontsize=8.5,
                 va="center", color=MUTED)
    fig.suptitle(f"Which cosmological parameters can be fixed in the mu(a) emulator?  z = {z:g}", x=0.01, ha="left",
                 fontsize=15, color=INK, y=0.975)
    fig.text(0.01, 0.885, "Dot: how much the mu - GR change D moves if the parameter is off by 1 sigma (Planck 2018), "
             "worst of the 7 gravity models and of the main spectra (11, dd, d2d2, s2s2); linear in the +-5 sigma runs "
             "(wa high: +3 sigma).\nGrey: the noise of one run (largest over the models). FIX: within the noise at "
             "every scale, so fixing it costs less than the noise every production run has anyway. VARY: otherwise.",
             fontsize=9.5, color=MUTED)
    pdf.savefig(fig)
    plt.close(fig)


if __name__ == "__main__":
    plt.rcParams.update({"figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "font.size": 10})
    path = os.path.join(ROOT, "figures", "mu_boost", "cosmology.pdf")
    with PdfPages(path) as pdf:
        verdict_page(pdf, 0.0)
        for z in (0.0, 1.0):
            overview_dots(pdf, z)
        for z in (0.0, 1.0):
            kdependence(pdf, z)
        for z in (0.0, 1.0):
            all_spectra_dots(pdf, z)
    print("Saved", path)
