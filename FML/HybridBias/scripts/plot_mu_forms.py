#!/usr/bin/env python3
"""mu(a) = 1 + A a^n (this project) against mu(a) = 1 + mu0 Omega_DE(a) / Omega_DE,0 (DES, Planck, DESI).
Left: the time dependence of mu - 1, normalised to 1 today. Right: its local power-law index n_eff = dln(mu - 1)/dln a
(constant for a^n). The Omega_DE form for LCDM and for the DESI w0wa best fit (w0 = -0.75, wa = -0.86, approximately).
Also prints the n for which A a^n with A = mu0 has the same growth as the Omega_DE form (scripts/mu_prior.py).
Writes figures/mu_boost/mu_forms.pdf.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

import mu_prior  # noqa: E402
from runs import ROOT  # noqa: E402

OM = mu_prior.OM


def omega_de_ratio(a, w0=-1.0, wa=0.0):
    """Omega_DE(a) / Omega_DE(a = 1) for a w0wa background (flat)."""
    rho_de = (1 - OM) * a ** (-3 * (1 + w0 + wa)) * np.exp(-3 * wa * (1 - a))
    return rho_de / (OM * a ** -3 + rho_de) / (1 - OM)


if __name__ == "__main__":
    a = np.logspace(-1.5, 0, 400)
    z = 1 / a - 1
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    blues = ["#86b6ef", "#4a90e2", "#2a78d6", "#0d3e80"]
    for n, c in zip([1, 2, 3, 4], blues):
        ax[0].plot(z, a ** n, color=c, lw=1.5, label=f"$a^{n}$ (ours, n = {n})")
        ax[1].plot(z, np.full_like(z, n), color=c, lw=1.5)
    for (w0, wa), ls, lab in [((-1.0, 0.0), "-", r"$\Omega_{DE}(a)/\Omega_{DE,0}$, $\Lambda$CDM"),
                              ((-0.75, -0.86), "--", r"$\Omega_{DE}(a)/\Omega_{DE,0}$, DESI $w_0w_a$")]:
        f = omega_de_ratio(a, w0, wa)
        ax[0].plot(z, f, color="#eb6834", lw=2, ls=ls, label=lab)
        ax[1].plot(z, np.gradient(np.log(f), np.log(a)), color="#eb6834", lw=2, ls=ls)
    for x in ax:
        x.set_xscale("log")
        x.set_xlim(0.01, 30)
        x.set_xlabel("z")
        x.grid(True, color="#e6e5e0", lw=0.6)
    ax[0].set_yscale("log")
    ax[0].set_ylim(1e-4, 1.5)
    ax[0].set_ylabel(r"$(\mu - 1) / (\mu - 1)_{z=0}$")
    ax[0].legend(fontsize=8.5, frameon=False)
    ax[0].set_title("time dependence of the modification", fontsize=10, loc="left")
    ax[1].set_ylim(0, 4.5)
    ax[1].set_ylabel(r"local index $n_{\rm eff} = d\ln(\mu-1)/d\ln a$")
    ax[1].set_title("how fast it switches on (constant for $A a^n$)", fontsize=10, loc="left")
    fig.tight_layout()
    path = os.path.join(ROOT, "figures", "mu_boost", "mu_forms.pdf")
    fig.savefig(path)
    print("Saved", path)

    # Which n gives A a^n (A = mu0) the same growth (fs8 over z = 0.3-1.5) as the Omega_DE form (LCDM background)?
    for mu0 in [-0.5, -0.2, 0.2, 0.5]:
        target = mu_prior.summary(mu_prior.mu_de(mu0))[0]
        n_eq = brentq(lambda n: mu_prior.summary(mu_prior.powerlaw(mu0, n))[0] - target, 0.3, 6)
        print(f"mu0 = {mu0:+.1f}: A a^n with A = mu0 has the same growth for n = {n_eq:.2f}")
