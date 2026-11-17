#!/usr/bin/env python3
"""Stitch theory (large scales) and simulation (small scales) hybrid bias spectra with one tanh weight in ln k:

    P_ab(k) = [1 - w(k)] P_ab^theory(k) + w(k) P_ab^sim(k),   w(k) = 1/2 [1 + tanh(ln(k / kt) / width)]

The same w(k) for all 15 spectra keeps the 5x5 matrix P_ab(k) positive semi-definite if both inputs are (a convex
combination), so every tracer spectrum sum_ab b_a b_b P_ab stays >= 0. The blend is additive, so spectra that cross
zero (e.g. delta s^2) are no problem. The same for the change mu - GR, D = (P_mu - P_GR) / sqrt(P_ii P_jj).

Tests (z = 0, the fiducial runs of run_mu_convergence.py, simulation = normal phase + control variate):
  1. positive semi-definiteness of the 5x5 matrix: simulation, theory and stitched
  2. smoothness: the largest jump in d ln|P| / d ln k across the transition, stitched vs the simulation alone
  3. sensitivity: the stitched spectra and D for kt = 0.012, 0.03, 0.05 and Zel'dovich vs one-loop theory
Writes figures/mu_boost/stitch.pdf (15-spectra layout; per mu model the change D and the error of the stitched D)
and prints the tables. With K10 K90, the transition goes from
w = 0.1 at k = K10 to w = 0.9 at k = K90 (kt = sqrt(K10 K90), width = ln(K90 / K10) / 2.2), and the PDF is
figures/mu_boost/stitch_K10_K90.pdf:

    python3 scripts/stitch.py [K10 K90]
"""
import os

import numpy as np

from plot_mu_theory import MODELS, PAIRS, load, norm

KT, WIDTH = 0.012, 0.3  # Just above the fundamental mode 0.006: the CV spectra are already Zel'dovich + (N-body - Zel'dovich field), nearly free of cosmic variance, so theory only extends them
NABLA = [p for p, (i, j) in enumerate(PAIRS) if 4 in (i, j)]


def weight(k, kt=None, width=None):
    """0 (theory) at k << kt, 1 (simulation) at k >> kt (defaults: KT, WIDTH)."""
    kt, width = kt or KT, width or WIDTH
    return 0.5 * (1.0 + np.tanh(np.log(k / kt) / width))


def stitch(k, P_sim, P_theory, kt=None, width=None):
    """P_sim, P_theory: [..., nk] on the same k bins. Returns the stitched spectra (same shape)."""
    w = weight(np.asarray(k), kt, width)
    return (1.0 - w) * P_theory + w * P_sim


def matrix(P):
    """[15, nk] -> [nk, 5, 5]"""
    M = np.zeros((P.shape[1], 5, 5))
    for p, (i, j) in enumerate(PAIRS):
        M[:, i, j] = M[:, j, i] = P[p]
    return M


def min_eigenvalue(P):
    """The smallest eigenvalue of the 5x5 matrix, relative to the largest diagonal element, per k."""
    M = matrix(P)
    lam = np.linalg.eigvalsh(M)[:, 0]
    return lam / np.max(np.abs(np.diagonal(M, axis1=1, axis2=2)), axis=1)


def slope_jump(k, P, sel):
    """The largest |change of d ln|P| / d ln k| between neighbouring bins in sel (a measure of kinks and noise)."""
    s = np.diff(np.log(np.abs(P[:, sel]) + 1e-300), axis=1) / np.diff(np.log(k[sel]))
    return np.max(np.abs(np.diff(s, axis=1)), axis=1)


def report():
    k, _ = load("GR", "sim", 0.0)
    window = (k > 0.015) & (k < 0.15)  # The transition region
    low = (k > 0.006) & (k < 0.03)
    print(f"Stitch: kt = {KT}, width = {WIDTH} in ln k (10-90% between k = "
          f"{KT * np.exp(-1.1 * WIDTH):.3f} and {KT * np.exp(1.1 * WIDTH):.3f} h/Mpc)\n")

    print("1. Smallest eigenvalue of the 5x5 matrix / largest diagonal, minimum over k in (0.006, 0.3) "
          "(< 0: some biases give a negative tracer spectrum)")
    sel = (k > 0.006) & (k < 0.3)
    print(f"  {'model':10s} {'sim':>10s} {'zel':>10s} {'oneloop':>10s} {'stitch zel':>11s} {'stitch 1L':>10s}")
    for m in MODELS:
        _, s = load(m, "sim", 0.0)
        out = []
        for P in [s, load(m, "zel", 0.0)[1], load(m, "oneloop", 0.0)[1]]:
            out.append(np.nan if P is None else np.min(min_eigenvalue(P)[sel]))
        for src in ("zel", "oneloop"):
            th = load(m, src, 0.0)[1]
            out.append(np.nan if th is None else np.min(min_eigenvalue(stitch(k, s, th))[sel]))
        print(f"  {m:10s}" + "".join(f"{v:+11.1e}" for v in out))

    print("\n2. Smoothness in the transition region (0.015 < k < 0.15): largest jump in d ln|P| / d ln k between "
          "neighbouring bins, GR\n   (sim alone | stitched with Zel'dovich | theory alone); the stitch must not add jumps")
    _, s = load("GR", "sim", 0.0)
    _, th = load("GR", "zel", 0.0)
    js, jt, jth = slope_jump(k, s, window), slope_jump(k, stitch(k, s, th), window), slope_jump(k, th, window)
    for p, (i, j) in enumerate(PAIRS):
        print(f"   P_{i}{j}: {js[p]:6.2f} | {jt[p]:6.2f} | {jth[p]:6.2f}")

    nrm = norm(s)
    print("\n3a. (stitched - simulation) / sqrt(P_ii P_jj) x 1e3, mean over 0.006 < k < 0.03 (what the theory "
          "replaces), GR")
    for src in ("zel", "oneloop"):
        _, th = load("GR", src, 0.0)
        r = 1e3 * np.mean(((stitch(k, s, th) - s) / nrm)[:, low], axis=1)
        print(f"   {src:8s}" + " ".join(f"{v:+6.2f}" for v in r))

    print("\n3b. Sensitivity of the stitched change mu - GR, D = (P_mu - P_GR) / sqrt(P_ii P_jj), to kt and the theory:"
          "\n    max |D - D_ref| x 1e3 over 0.006 < k < 0.3, ref = Zel'dovich at kt = %g; nabla^2 spectra left out" % KT)
    keep = [p for p in range(15) if p not in NABLA]
    for m in MODELS[1:]:
        ref = D_stitched(m, "zel", KT)
        row = []
        for src, kt in [("zel", 0.03), ("zel", 0.05), ("oneloop", KT), ("oneloop", 0.05)]:
            d = D_stitched(m, src, kt)
            row.append(np.nan if d is None else 1e3 * np.max(np.abs(d - ref)[keep][:, sel]))
        print(f"   {m:8s} zel kt0.03 {row[0]:6.2f} | zel kt0.05 {row[1]:6.2f} | 1L kt{KT:g} {row[2]:6.2f} | "
              f"1L kt0.05 {row[3]:6.2f}")


def D_stitched(m, src, kt=None):
    """The stitched change D = (P_mu - P_GR) / sqrt(P_ii P_jj) (GR simulation diagonals), or None without theory."""
    k, sg = load("GR", "sim", 0.0)
    _, s = load(m, "sim", 0.0)
    th, thg = load(m, src, 0.0)[1], load("GR", src, 0.0)[1]
    if th is None or thg is None:
        return None
    return stitch(k, (s - sg) / norm(sg), (th - thg) / norm(sg), kt=kt or KT)


def plots():
    """Per mu model: the change D = (P_mu - P_GR) / sqrt(P_ii P_jj) (simulation and stitched) and the error of the
    stitched D, D_stitched - D_sim."""
    from plot_page import save_pages

    figures = []
    k, sg = load("GR", "sim", 0.0)
    span = f"10-90% from k = {K10:.3f} to {K90:.3f}"
    for m in MODELS[1:]:
        _, s = load(m, "sim", 0.0)
        name = "mu = 1 + " + m.replace("_n", " a^").replace("A", "")
        D_sim = (s - sg) / norm(sg)
        bc, err = [{"k": k, "P": D_sim, "label": "simulation (normal + CV)", "color": "#0b0b0b"}], []
        for src, colour in (("zel", "#eb6834"), ("oneloop", "#2a78d6")):
            d = D_stitched(m, src)
            if d is not None:
                bc.append({"k": k, "P": d, "label": f"stitched with {src}, {span}", "color": colour})
                err.append({"k": k, "P": d - D_sim, "label": f"stitched with {src} - simulation", "color": colour})
        figures.append({"title": f"{name}, z = 0: change mu - GR, D = (P_mu - P_GR) / sqrt(P_ii P_jj) "
                                 "(y range per panel)",
                        "ratio": True, "line": 0.0, "ylabel": r"$(P_\mu - P_{\rm GR}) / \sqrt{P_{ii} P_{jj}}$",
                        "kmin": 0.005, "kmax": 1.0, "curves": bc})
        figures.append({"title": f"{name}, z = 0: error of the stitched D, D_stitched - D_sim "
                                 "(grey band +-0.001; y range per panel)",
                        "ratio": True, "line": 0.0, "band": 0.001, "ylabel": r"$D_{\rm stitched} - D_{\rm sim}$",
                        "kmin": 0.005, "kmax": 0.2, "curves": err})
    save_pages(figures, f"mu_boost/{OUT}")


if __name__ == "__main__":
    import sys

    OUT = "stitch.pdf"
    if len(sys.argv) == 3:
        k10, k90 = (float(a) for a in sys.argv[1:])
        KT, WIDTH, OUT = np.sqrt(k10 * k90), np.log(k90 / k10) / 2.2, f"stitch_{sys.argv[1]}_{sys.argv[2]}.pdf"
    K10, K90 = KT * np.exp(-1.1 * WIDTH), KT * np.exp(1.1 * WIDTH)
    report()
    plots()
