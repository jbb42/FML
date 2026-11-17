#!/usr/bin/env python3
"""The Zel'dovich and one-loop (CLEFT, EdS kernels) predictions against the simulations of run_mu_convergence.py,
for GR and mu(a) = 1 + A a^n. Simulation: the fiducial normal-phase run with the control variate. Zel'dovich: from the
same run. One-loop: results/mu_convergence/theory/<model>_EdS (results/mu_convergence/theory/run_theory.py).

Writes (z = 0)
    figures/mu_boost/theory_vs_sim.pdf          all 15 spectra, per model: spectra, (simulation - theory), and (mu
                                                models) the change mu - GR from the simulation and from both theories
    figures/mu_boost/theory_vs_sim_summary.pdf  one page: up to which k each theory reproduces the simulation (five key
                                                spectra within 2%, the change mu - GR within 0.005), all models
Differences are normalised, D_ij = Delta P_ij / sqrt(P_ii P_jj) with the diagonal GR simulation spectra (as in
plot_fixes.py), so spectra that cross zero (e.g. delta s^2) give no spikes, unlike ratios.

    python3 scripts/plot_mu_theory.py
"""
import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

from plot_page import save_pages  # noqa: E402
from runs import ROOT  # noqa: E402

MODELS = ["GR", "A0.1_n3", "A0.3_n3", "A1.0_n3", "A0.1_n1", "A0.3_n1", "A1.0_n1"]
PAIRS = [(i, j) for i in range(5) for j in range(i, 5)]
NAMES = [r"\delta_m", r"\delta_L", r"\delta_L^2", r"s^2", r"\nabla^2\delta_L"]
KEY = [(0, 0), (1, 1), (2, 2), (3, 3), (0, 4)]  # The five spectra of the model pages
KMIN, KMAX = 0.01, 1.0
TOL_SPECTRUM, TOL_BOOST = 0.02, 0.005  # In units of sqrt(P_ii P_jj)
# Reference palette (dataviz skill), fixed order: slot 1 blue, slot 2 orange; ink for the simulation
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
STYLE = {"sim": dict(color=INK, lw=1.6, label="simulation (normal + CV)"),
         "oneloop": dict(color="#2a78d6", lw=1.8, label="one-loop (CLEFT, EdS kernels)"),
         "zel": dict(color="#eb6834", lw=1.8, ls="--", label="Zel'dovich")}
FIGURES = os.path.join(ROOT, "figures", "mu_boost")


def title(model):
    return "GR" if model == "GR" else rf"$\mu = 1 + {model.split('_')[0][1:]}\,a^{{{model.split('_')[1][1:]}}}$"


def pair_name(p):
    i, j = PAIRS[p] if isinstance(p, int) else p
    return rf"$P_{{{NAMES[i]}\,{NAMES[j]}}}$"


def load(model, source, z):
    """k and [15, nk]: source 'sim' (with the control variate), 'zel' (Zel'dovich, from the run), 'oneloop'."""
    prefix = "gr" if model == "GR" else f"mu_{model}"
    run = os.path.join(ROOT, "results", "mu_convergence", f"{prefix}_fiducial_L1024_f2048_n1024_s30_normal")
    path = {"sim": os.path.join(run, f"snapshot_*_z{z:.3f}", f"hybridbias_cv_z{z:.3f}.txt"),
            "zel": os.path.join(run, f"snapshot_*_z{z:.3f}", f"hybridbias_zeldovich_z{z:.3f}.txt"),
            "oneloop": os.path.join(ROOT, "results", "mu_convergence", "theory", f"{model}_EdS",
                                    f"snapshot_*_z{z:.3f}", f"hybridbias_oneloop_z{z:.3f}.txt")}[source]
    files = glob.glob(path)
    if not files:
        return None, None
    data = np.loadtxt(files[0])
    return data[:, 0], data[:, 2:17].T


def norm(P):
    """[15, nk]: sqrt(|P_ii P_jj|) for each pair, from the diagonal spectra of P."""
    d = [np.abs(P[PAIRS.index((i, i))]) for i in range(5)]
    return np.array([np.sqrt(d[i] * d[j]) for i, j in PAIRS])


def agreement_k(k, x, tol):
    """The largest k such that |x| <= tol in every bin from KMIN up to it (NaN: not even at KMIN)."""
    sel = k >= KMIN
    ok = np.abs(np.nan_to_num(x[sel], nan=np.inf)) <= tol
    if not ok[0]:
        return np.nan
    bad = np.where(~ok)[0]
    return k[sel][bad[0] - 1] if bad.size else k[sel][-1]


def style_axis(ax):
    ax.set_xscale("log")
    ax.set_xlim(KMIN, KMAX)
    ax.grid(True, color=GRID, lw=0.6)
    ax.tick_params(colors=MUTED, labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)


def summary_page(pdf, z):
    """Two tables, models x key spectra: the k up to which each theory reproduces the spectrum and the boost."""
    rows = {"spectrum": [], "boost": []}
    for model in MODELS:
        k, sim = load(model, "sim", z)
        _, sim_gr = load("GR", "sim", z)
        cells = {"spectrum": [], "boost": []}
        for p in KEY:
            i = PAIRS.index(p)
            for kind in ("spectrum", "boost"):
                entry = []
                for source in ("zel", "oneloop"):
                    _, th = load(model, source, z)
                    _, th_gr = load("GR", source, z)
                    if th is None or (kind == "boost" and (model == "GR" or th_gr is None)):
                        entry.append(None)
                    elif kind == "spectrum":
                        entry.append(agreement_k(k, (th[i] - sim[i]) / norm(sim_gr)[i], TOL_SPECTRUM))
                    else:
                        entry.append(agreement_k(k, ((th[i] - th_gr[i]) - (sim[i] - sim_gr[i])) / norm(sim_gr)[i],
                                                 TOL_BOOST))
                cells[kind].append(entry)
        for kind in rows:
            rows[kind].append(cells[kind])

    fig, axes = plt.subplots(2, 1, figsize=(11, 8.5))
    fig.suptitle(f"Up to which k (h/Mpc) does the theory reproduce the simulation? z = {z}", fontsize=13, x=0.06,
                 ha="left", color=INK)
    titles = {"spectrum": f"Spectrum: |P_theory - P_sim| / sqrt(P_ii P_jj) <= {TOL_SPECTRUM} for all k up to ...",
              "boost": f"Change mu - GR: |D_theory - D_sim| <= {TOL_BOOST}, D = (P_mu - P_GR) / sqrt(P_ii P_jj), "
                       "for all k up to ..."}
    for ax, kind in zip(axes, ("spectrum", "boost")):
        ax.axis("off")
        ax.set_title(titles[kind], loc="left", fontsize=10.5, color=INK)
        cell_text, cell_colour = [], []
        models = [m for m in MODELS if not (kind == "boost" and m == "GR")]
        for model, entries in zip(MODELS, rows[kind]):
            if model not in models:
                continue
            text, colour = [], []
            for zel, one in entries:
                fmt = lambda v: "-" if v is None else ("< 0.01" if np.isnan(v) else f"{v:.2f}")
                text.append(f"{fmt(zel)}  |  {fmt(one)}")
                best = max([v for v in (zel, one) if v is not None and not np.isnan(v)], default=KMIN)
                shade = np.clip(np.log10(best / KMIN) / np.log10(KMAX / KMIN), 0, 1)  # Sequential blue by reach
                colour.append(matplotlib.colors.to_hex((1 - 0.35 * shade, 1 - 0.2 * shade, 1)))
            cell_text.append(text)
            cell_colour.append(colour)
        labels = [m if m == "GR" else f"A = {m.split('_')[0][1:]}, n = {m.split('_')[1][1:]}" for m in models]
        table = ax.table(cellText=[[lab] + t for lab, t in zip(labels, cell_text)],
                         cellColours=[["white"] + c for c in cell_colour], colLabels=["model"] + [pair_name(p) for p in KEY],
                         colWidths=[0.16] + [0.168] * len(KEY), loc="upper left", cellLoc="center",
                         bbox=[0.0, 0.0, 1.0, 0.88])
        table.auto_set_font_size(False)
        table.set_fontsize(9.5)
        for cell in table.get_celld().values():
            cell.set_edgecolor("white")
    fig.text(0.06, 0.02, "Each cell: Zel'dovich  |  one-loop. '< 0.01': off already at k = 0.01; '-': not available. "
             "Darker: the better of the two reaches further.", fontsize=9, color=MUTED)
    pdf.savefig(fig)
    plt.close(fig)


def all15_pages(z):
    """The 15-spectra layout of plot_page.py: per model a page with the spectra, one with simulation / theory and (mu
    models) one with the boosts."""
    figures = []
    _, sim_gr = load("GR", "sim", z)
    th_gr = {s: load("GR", s, z)[1] for s in ("oneloop", "zel")}
    for model in MODELS:
        k, sim = load(model, "sim", z)
        th = {s: load(model, s, z)[1] for s in ("oneloop", "zel")}
        name = "GR" if model == "GR" else "mu = 1 + " + model.replace("_n", " a^").replace("A", "")
        curves = [{"k": k, "P": sim, "label": STYLE["sim"]["label"], "color": INK}]
        nrm = norm(sim_gr)
        ratios, boosts = [], [{"k": k, "P": (sim - sim_gr) / nrm, "label": "simulation", "color": INK}]
        for s in ("oneloop", "zel"):
            if th[s] is None:
                continue
            look = {"color": STYLE[s]["color"], "ls": STYLE[s].get("ls", "-")}
            curves.append({"k": k, "P": th[s], "label": STYLE[s]["label"], **look})
            ratios.append({"k": k, "P": (sim - th[s]) / nrm, "label": "simulation - " + STYLE[s]["label"], **look})
            if th_gr[s] is not None:
                boosts.append({"k": k, "P": (th[s] - th_gr[s]) / nrm, "label": STYLE[s]["label"], **look})
        figures.append({"title": f"{name}, z = {z}: simulation (normal phase + control variate) and theory",
                        "kmin": KMIN, "kmax": KMAX, "curves": curves})
        figures.append({"title": f"{name}, z = {z}: (simulation - theory) / sqrt(P_ii P_jj), GR simulation diagonals "
                                 "(grey band: +-0.02; y range per panel)",
                        "ratio": True, "line": 0.0, "band": 0.02,
                        "ylabel": r"$(P_{\rm sim} - P_{\rm theory}) / \sqrt{P_{ii} P_{jj}}$", "kmin": KMIN,
                        "kmax": KMAX, "curves": ratios})
        if model != "GR":
            figures.append({"title": f"{name}, z = {z}: change mu - GR, D = (P_mu - P_GR) / sqrt(P_ii P_jj), "
                                     "simulation and theory (y range per panel)",
                            "ratio": True, "line": 0.0, "ylabel": r"$(P_\mu - P_{\rm GR}) / \sqrt{P_{ii} P_{jj}}$",
                            "kmin": KMIN, "kmax": KMAX, "curves": boosts})
    save_pages(figures, "mu_boost/theory_vs_sim.pdf")


if __name__ == "__main__":
    plt.rcParams.update({"font.size": 10, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})
    os.makedirs(FIGURES, exist_ok=True)
    with PdfPages(os.path.join(FIGURES, "theory_vs_sim_summary.pdf")) as pdf:
        summary_page(pdf, 0.0)
    print("Saved", os.path.join(FIGURES, "theory_vs_sim_summary.pdf"))
    all15_pages(0.0)
