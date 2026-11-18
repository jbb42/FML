#!/usr/bin/env python3
"""How much does each cosmological parameter matter for the mu - GR change of the spectra, for each gravity model?
From the full-scale cosmology test (results/mu_cosmology, run_all.py and run_extra.py): every parameter at -5 and +5
sigma (wa: +3), joint extremes, Omega_b, and three seeds of the fiducial cosmology for the noise.

D_ij = (P_mu - P_GR) / sqrt(P_ii^GR P_jj^GR) (GR and mu from the same seed in each cosmology, no control variate).
The effect of a parameter: |D(parameter at +-5 sigma) - D(fiducial)|, the larger of -5 and +5, averaged over a k range,
in % of max_k |D(fiducial)|, the worst of the main spectra (11, dd, d2d2, s2s2). The noise: the same for one seed
against another, divided by sqrt(2) (the noise of one run).

Writes figures/mu_boost/cosmology.pdf (everything on one scale: the change of D for a realistic 1-sigma shift):
  page 1    summary (z = 0): which parameters can be fixed (worst model and main spectrum, against the noise)
  page 2-4  z = 0: per gravity model (and the joint extremes as run); all 15 spectra in the 15-spectra layout;
            k dependence of the main spectra
  page 5-7  the same for z = 1
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


_DCACHE = {}


def D_cached(z, label, model):
    key = (z, label, model)
    if key not in _DCACHE:
        A.Z = z
        _DCACHE[key] = A.D(label, model)
    return _DCACHE[key]


# ---------------------------------------------------------------------------------------------- dot plots
N_COLOURS = {1.0: "#86b6ef", 2.0: "#3987e5", 3.0: "#1c5cab", 4.0: "#0d366b"}  # Sequential blue by n (light: early MG)
SCALE_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a"]  # Large, intermediate, small scales (categorical slots 1-3)


def _style(ax, xlim=(0.03, 100)):
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.grid(True, axis="x", color=GRID, lw=0.6)
    ax.tick_params(colors=MUTED, labelsize=8, length=0)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))


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


# ---------------------------------------------------------------------------------------------- detail pages
OPS = [r"1", r"\delta", r"\delta^2", r"s^2", r"\nabla^2\delta"]
NOISE_FILL = "#e3e2dc"


def noise_area(ax, nz, y_text, xmin=0.01):
    ax.axvspan(xmin, nz, color=NOISE_FILL, lw=0, zorder=0)
    ax.text(np.sqrt(xmin * nz), y_text, "within the noise\nof one run", ha="center", va="center", fontsize=8,
            color=MUTED, style="italic")


def pct(v):
    return f"{v:.2f}%" if v < 1 else f"{v:.1f}%"


def model_page(pdf, z):
    """Per gravity model: the change of D for a 1-sigma shift (main spectra), and the joint extremes as run."""
    ms = models()
    ypos = list(range(len(PARAMS))) + [len(PARAMS) + 1.0 + i for i in range(len(JOINT))]
    fig, axes = plt.subplots(1, len(RANGES), figsize=(18, 9.5), sharey=True)
    fig.subplots_adjust(left=0.12, right=0.99, top=0.82, bottom=0.08, wspace=0.08)
    for ax, (title, lo, hi) in zip(axes, RANGES):
        noise_area(ax, max(noise(z, m, lo, hi) for m, _, _ in ms), -1.0)
        for y, (name, lo_l, hi_l) in zip(ypos, PARAMS):
            E = np.array([effect_1sigma(z, lo_l, hi_l, m, lo, hi) for m, _, _ in ms])
            ax.plot([E.min(), E.max()], [y, y], color=MUTED, lw=1.2, zorder=1)
            for (m, a, n), e in zip(ms, E):
                ax.plot(e, y, "o" if a > 0 else "s", ms=9, color=N_COLOURS[n], mec="white", mew=1.3, zorder=2)
            ax.text(E.max() * 1.35, y, pct(E.max()), va="center", fontsize=9, color=INK)
        for y, (name, lab) in zip(ypos[len(PARAMS):], JOINT):
            E = np.array([effect(z, lab, m, lo, hi) for m, _, _ in ms])
            ax.plot([E.min(), E.max()], [y, y], color=MUTED, lw=1.2, zorder=1)
            for (m, a, n), e in zip(ms, E):
                ax.plot(e, y, "o" if a > 0 else "s", ms=7, color=N_COLOURS[n], mec="white", mew=1.0, alpha=0.8,
                        zorder=2)
            ax.text(min(E.max() * 1.35, 25), y, pct(E.max()), va="center", fontsize=8.5, color=MUTED)
        ax.axhline(len(PARAMS) + 0.15, color=GRID, lw=1)
        ax.text(0.012, len(PARAMS) + 0.45, "joint extremes, as run (several parameters at +-5 sigma; not rescaled)",
                fontsize=8.5, color=MUTED, va="center")
        _style(ax, (0.01, 40))
        ax.set_ylim(ypos[-1] + 0.6, -1.7)
        ax.set_title(f"{title.split(',')[0]}, {lo:g} < k < {hi:g}", fontsize=11.5, color=INK, loc="left")
        ax.set_xlabel("change of D, % of its size", color=MUTED, fontsize=9)
    axes[0].set_yticks(ypos)
    axes[0].set_yticklabels([n for n, _, _ in PARAMS] + [n for n, _ in JOINT], fontsize=10.5, color=INK)
    handles = [matplotlib.lines.Line2D([], [], ls="", marker="o", ms=8, color=N_COLOURS[n], label=f"n = {n:g}")
               for n in sorted(N_COLOURS)]
    handles += [matplotlib.lines.Line2D([], [], ls="", marker="o", ms=8, color=MUTED, label="A > 0"),
                matplotlib.lines.Line2D([], [], ls="", marker="s", ms=8, color=MUTED, label="A < 0")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.12, 0.895), ncol=6, frameon=False, fontsize=10)
    fig.suptitle(f"Per gravity model: how much a 1-sigma shift changes the mu - GR change D  (z = {z:g})", x=0.01,
                 ha="left", fontsize=15, color=INK, y=0.975)
    fig.text(0.01, 0.925, "One dot per gravity model (colour: n, light = early-acting; shape: sign of A), worst of the "
             "main spectra (11, dd, d2d2, s2s2) and of the two signs; the number is the largest. Same scale as the "
             "summary.", fontsize=9.5, color=MUTED)
    pdf.savefig(fig)
    plt.close(fig)


def spectra_page(pdf, z):
    """All 15 spectra in the 15-spectra layout: change of D for a 1-sigma shift (worst gravity model) at the three
    scales; grey: the noise of one run for that spectrum (small scales)."""
    ms = [m for m, _, _ in models()]
    seeds = ("fiducial_s1002", "fiducial_s1003")
    fig, axes = plt.subplots(5, 5, figsize=(19, 14))
    fig.subplots_adjust(left=0.06, right=0.99, top=0.89, bottom=0.04, wspace=0.12, hspace=0.42)
    for ax in axes.flat:
        ax.axis("off")
    for p, (i, j) in enumerate(A.PAIRS):
        ax = axes[i, j]
        ax.axis("on")
        nz = 0.0
        for m in ms:
            k, Df = D_cached(z, "fiducial", m)
            sel = (k > 0.6) & (k <= 1.0)
            scale = np.max(np.abs(Df[p][(k > 0.01) & (k < 1.0)]))
            nz = max(nz, np.mean([100 * np.mean(np.abs(D_cached(z, s_, m)[1][p][sel] - Df[p][sel])) / scale
                                  for s_ in seeds]) / np.sqrt(2))
        ax.axvspan(0.01, nz, color=NOISE_FILL, lw=0, zorder=0)
        for r, (name, lo_l, hi_l) in enumerate(PARAMS):
            vals = [max(effect_1sigma(z, lo_l, hi_l, m, lo, hi, [p]) for m in ms) for _, lo, hi in RANGES]
            ax.plot([0.01, max(vals)], [r, r], color=GRID, lw=1, zorder=1)
            for c, v in zip(SCALE_COLOURS, vals):
                ax.plot(v, r, "o", ms=6, color=c, mec="white", mew=0.8, zorder=2)
        _style(ax, (0.01, 40))
        ax.set_ylim(len(PARAMS) - 0.5, -0.5)
        ax.set_yticks(range(len(PARAMS)))
        ax.set_yticklabels([n for n, _, _ in PARAMS] if j == i else [], fontsize=9, color=INK)
        ax.set_title(rf"${OPS[i]} \times {OPS[j]}$", fontsize=11, color=INK, pad=3)
    handles = [matplotlib.lines.Line2D([], [], ls="", marker="o", ms=8, color=c,
                                       label=f"{t.split(',')[0]}, {lo:g} < k < {hi:g}")
               for c, (t, lo, hi) in zip(SCALE_COLOURS, RANGES)]
    handles.append(matplotlib.patches.Patch(color=NOISE_FILL, label="within the noise of one run (small scales)"))
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.06, 0.06), frameon=False, fontsize=11)
    fig.suptitle(f"All 15 spectra: how much a 1-sigma shift changes the mu - GR change D  (z = {z:g})", x=0.01,
                 ha="left", fontsize=15, color=INK, y=0.975)
    fig.text(0.01, 0.935, "x: change of D in % of its size (log scale), worst gravity model and sign. Same scale as the "
             "summary.", fontsize=10, color=MUTED)
    pdf.savefig(fig)
    plt.close(fig)


def k_page(pdf, z):
    """k dependence: signed change of D per 1-sigma shift, (D_hi - D_lo) / (sigma_lo + sigma_hi), per gravity model;
    grey: +- the noise of one run; rows share the y axis so the spectra compare directly."""
    ms = models()
    fig, axes = plt.subplots(len(PARAMS), len(MAIN), figsize=(16, 18), sharex=True, sharey="row")
    fig.subplots_adjust(left=0.08, right=0.99, top=0.93, bottom=0.07, hspace=0.3, wspace=0.08)
    for r, (name, lo_l, hi_l) in enumerate(PARAMS):
        s_lo, s_hi = SIGMAS[lo_l.rsplit("_", 1)[0]]
        for c, q in enumerate(MAIN):
            ax = axes[r, c]
            i = A.PAIRS.index(q)
            nzs = []
            for m, a, n in ms:
                k, Df = D_cached(z, "fiducial", m)
                scale = np.max(np.abs(Df[i][(k > 0.01) & (k < 1.0)]))
                y = 100 * (D_cached(z, hi_l, m)[1][i] - D_cached(z, lo_l, m)[1][i]) / (s_lo + s_hi) / scale
                ax.plot(k, y, color=N_COLOURS[n], lw=1.4, ls="-" if a > 0 else "--")
                seeds = np.array([D_cached(z, s_, m)[1][i] for s_ in ("fiducial", "fiducial_s1002", "fiducial_s1003")])
                nzs.append(100 * np.std(seeds, axis=0, ddof=1) / scale)
            nz = np.convolve(np.max(nzs, axis=0), np.ones(5) / 5, mode="same")
            ax.fill_between(k, -nz, nz, color=NOISE_FILL, lw=0, zorder=0)
            ax.axhline(0, color=MUTED, lw=0.6)
            ax.set_xscale("log")
            ax.set_xlim(0.01, 1.0)
            ax.grid(True, color=GRID, lw=0.5)
            ax.tick_params(labelsize=8.5, colors=MUTED, length=0)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            if r == 0:
                ax.set_title(NAMES[q], fontsize=13, color=INK)
            if c == 0:
                ax.set_ylabel(name, fontsize=13, color=INK, rotation=0, ha="right", va="center")
            if r == len(PARAMS) - 1:
                ax.set_xlabel(r"$k$ [$h$/Mpc]", color=INK)
    handles = [matplotlib.lines.Line2D([], [], color=N_COLOURS[n], lw=2, label=f"n = {n:g}") for n in sorted(N_COLOURS)]
    handles += [matplotlib.lines.Line2D([], [], color=MUTED, lw=1.5, ls="-", label="A > 0"),
                matplotlib.lines.Line2D([], [], color=MUTED, lw=1.5, ls="--", label="A < 0"),
                matplotlib.patches.Patch(color=NOISE_FILL, label="+- noise of one run")]
    fig.legend(handles=handles, loc="lower center", ncol=7, fontsize=10.5, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle(f"Where in k: signed change of D for a +1-sigma shift, in % of its size  (z = {z:g})", x=0.01,
                 ha="left", fontsize=15, color=INK, y=0.985)
    fig.text(0.01, 0.955, "One line per gravity model; each row shares its y axis, so the spectra compare directly. "
             "Same scale as the summary (there: the absolute value averaged over each k range, worst model).",
             fontsize=10, color=MUTED)
    pdf.savefig(fig)
    plt.close(fig)


if __name__ == "__main__":
    plt.rcParams.update({"figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "font.size": 10})
    path = os.path.join(ROOT, "figures", "mu_boost", "cosmology.pdf")
    with PdfPages(path) as pdf:
        verdict_page(pdf, 0.0)
        for z in (0.0, 1.0):
            model_page(pdf, z)
            spectra_page(pdf, z)
            k_page(pdf, z)
    print("Saved", path)
