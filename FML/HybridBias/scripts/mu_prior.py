#!/usr/bin/env python3
"""Realistic ranges of A and n in mu(a) = 1 + A a^n, from linear growth (Planck 2018 LCDM background, GR initial
conditions as in the simulations).

1. Equivalent mu0: the DESI 2024 analysis (arXiv:2411.12026) uses mu(a) = 1 + mu0 Omega_DE(a) / Omega_DE,0. For each
   (A, n) we find the mu0 that gives the same growth over the range the data measure: we match the mean of
   f sigma_8(z) / f sigma_8^GR(z) over z = 0.3-1.5 (the RSD redshifts), and also report the z = 0 growth D_mu / D_GR.
   Bounds on mu0 (68%): 0.05 +- 0.22 (DESI FS+BAO + CMB + DES Y3 + DESY5 SN, LCDM background),
   -0.24 +0.32 -0.28 (same, w0waCDM background).
2. Direct check: chi^2 of f sigma_8(z) against the eBOSS DR16 consensus RSD compilation (6dFGS, MGS, BOSS DR12, eBOSS
   LRG, ELG, QSO), with sigma_8 fixed by Planck (sigma_8^GR = 0.811 at z = 0; MG only changes the late growth).

Writes figures/mu_boost/mu_prior.pdf and prints the allowed region.
    python3 scripts/mu_prior.py
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

OM, SIGMA8 = 0.3153, 0.811  # Planck 2018 (TT,TE,EE+lowE+lensing)
OL = 1 - OM
# eBOSS DR16 consensus (Alam et al. 2021, Table 3) f sigma_8: z, value, error (diagonal; correlations neglected)
FS8 = np.array([[0.067, 0.423, 0.055], [0.15, 0.53, 0.16], [0.38, 0.497, 0.045], [0.51, 0.459, 0.038],
                [0.70, 0.473, 0.041], [0.85, 0.315, 0.095], [1.48, 0.462, 0.045]])
ZFIT = np.linspace(0.3, 1.5, 13)


def E2(a):
    return OM / a ** 3 + OL


def growth(mu, zs):
    """D(z) / D_GR(z = 0) and f(z) for mu(a), starting in matter domination at a = 1e-3 (D = a, GR IC)."""
    def rhs(lna, y):
        a = np.exp(lna)
        D, dD = y  # dD = dD / dln a
        dlnE = -1.5 * OM / a ** 3 / E2(a)  # dln E / dln a
        return [dD, -(2 + dlnE) * dD + 1.5 * OM / a ** 3 / E2(a) * mu(a) * D]

    a0 = 1e-3
    sol = solve_ivp(rhs, [np.log(a0), 0.0], [a0, a0], dense_output=True, rtol=1e-8, atol=1e-12)
    D, dD = sol.sol(np.log(1 / (1 + np.asarray(zs, float))))
    return D, dD / D, sol.sol(0.0)[0]


ZALL = np.concatenate([ZFIT, FS8[:, 0], [0.0]])
D_GR, f_GR, D0_GR = growth(lambda a: 1.0, ZALL)


def fs8(mu):
    D, f, _ = growth(mu, ZALL)
    return f * SIGMA8 * D / D0_GR  # sigma_8 normalised at early times (as the simulations: GR initial conditions)


FS8_GR = fs8(lambda a: 1.0)
nfit = len(ZFIT)


def summary(mu):
    """(mean fs8 / fs8_GR - 1 over the RSD range, D_mu / D_GR - 1 at z = 0, chi^2 of the RSD data)"""
    v = fs8(mu)
    _, _, D0 = growth(mu, [0.0])
    data = v[nfit:nfit + len(FS8)]
    return np.mean(v[:nfit] / FS8_GR[:nfit]) - 1, D0 / D0_GR - 1, np.sum(((data - FS8[:, 1]) / FS8[:, 2]) ** 2)


def mu_de(mu0):
    return lambda a: 1 + mu0 * (OL / (OM / a ** 3 + OL)) / OL


def powerlaw(A, n):
    return lambda a: 1 + A * a ** n


if __name__ == "__main__":
    # mu0 -> growth change, to invert
    mu0s = np.linspace(-1.5, 2.0, 71)
    g_mu0 = np.array([summary(mu_de(m))[0] for m in mu0s])
    to_mu0 = lambda g: np.interp(g, g_mu0, mu0s)  # Monotonic in mu0
    chi2_gr = summary(lambda a: 1.0)[2]

    As = np.linspace(-1.0, 2.0, 61)
    ns = np.linspace(0.25, 6.0, 47)
    G = np.zeros((len(ns), len(As)))
    D0 = np.zeros_like(G)
    C = np.zeros_like(G)
    for i, n in enumerate(ns):
        for j, A in enumerate(As):
            G[i, j], D0[i, j], C[i, j] = summary(powerlaw(A, n))
    M0 = to_mu0(G)
    dchi = C - C.min()

    print(f"chi^2 of the RSD data (7 points): GR {chi2_gr:.2f}, best (A, n) {C.min():.2f}")
    bounds = {"LCDM 2 sigma": (0.05 - 0.44, 0.05 + 0.44), "w0waCDM 2 sigma": (-0.24 - 0.56, -0.24 + 0.64),
              "union 2 sigma": (-0.80, 0.49), "union 3 sigma": (-1.08, 0.71)}
    for label, (lo, hi) in bounds.items():
        print(f"\nAllowed by mu0_equiv in [{lo:+.2f}, {hi:+.2f}] ({label}): A range for each n")
        for n in [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0]:
            i = np.argmin(np.abs(ns - n))
            ok = (M0[i] >= lo) & (M0[i] <= hi)
            print(f"   n = {ns[i]:4.2f}: A in [{As[ok].min():+.2f}, {As[ok].max():+.2f}]" if ok.any() else
                  f"   n = {ns[i]:4.2f}: none")
    print("\nCheck: mu0_equiv of the simulated models, with the z = 0 growth and the RSD Delta chi^2 (vs best)")
    for A, n in [(0.1, 3), (0.3, 3), (1.0, 3), (0.1, 1), (0.3, 1), (1.0, 1), (-0.5, 0.5), (-0.5, 4), (1, 0.5), (1, 4)]:
        g, d, c = summary(powerlaw(A, n))
        print(f"   A = {A:+.1f}, n = {n:3.1f}: mu0_equiv {to_mu0(g):+.2f}, D/D_GR - 1 = {100 * d:+5.1f}%, "
              f"fs8 change {100 * g:+5.1f}%, Delta chi^2 {c - C.min():5.1f}")

    # One panel: shaded = excluded (equivalent mu0 outside 3 sigma); lines = growth change today; box = prior
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    excluded = (M0 < -1.08) | (M0 > 0.71)
    ax.contourf(As, ns, excluded.astype(float), levels=[0.5, 1.5], colors=["#e6e5e0"])
    ax.contour(As, ns, M0, levels=[-1.08, 0.71], colors="#0b0b0b", linewidths=1.6)
    ax.contour(As, ns, M0, levels=[-0.39, 0.49], colors="#0b0b0b", linewidths=1.0, linestyles="--")
    cs = ax.contour(As, ns, 100 * D0, levels=[-10, -5, -2, -1, 1, 2, 5, 10], colors="#2a78d6", linewidths=0.8, negative_linestyles="solid")
    ax.clabel(cs, fmt=lambda v: f"{v:+g}%", fontsize=8)
    ax.add_patch(matplotlib.patches.Rectangle((-0.5, 1), 1.0, 3, fill=False, ec="#eb6834", lw=2))
    ax.text(0.52, 4.05, "proposed prior", color="#eb6834", fontsize=10, va="bottom")
    sims = [(0.3, 3), (-0.35, 1), (0.2, 1), (0.6, 2), (-0.5, 2), (-0.5, 4), (1, 4),  # cosmology test
            (0.1, 3), (1.0, 3), (0.1, 1), (0.3, 1), (1.0, 1), (-0.5, 0.5), (1, 0.5)]   # earlier runs
    for A, n in sims:
        out = excluded[np.argmin(np.abs(ns - n)), np.argmin(np.abs(As - A))]
        ax.plot(A, n, "x" if out else "*", color="#0b0b0b", ms=9 if not out else 7)
    ax.set_xlim(-1, 2)
    ax.set_xlabel("A  (strength today: mu(z = 0) = 1 + A)")
    ax.set_ylabel("n  (how late it switches on; large n = only at late times)")
    ax.set_title("Which mu = 1 + A a^n are allowed by the growth of structure (DESI 2024 + CMB + SN)?\n"
                 "grey: excluded (> 3 sigma); solid / dashed lines: 3 / 2 sigma edges; blue: growth today vs GR; "
                 "stars: allowed runs, x: excluded runs", fontsize=9, loc="left")
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures", "mu_boost",
                        "mu_prior.pdf")
    fig.savefig(path, bbox_inches="tight")
    print("Saved", path)
