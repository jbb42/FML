#!/usr/bin/env python3
"""The Planck 2018 cosmology, and Omega_m, H_0, A_s, w0 and wa each changed on its own by -5, -1, +1 and +5 sigma while
the others keep their Planck values (wa: +3 instead of +5 sigma, see below). Every cosmology is run in GR and f(R)
(F5), from the same initial conditions, with one seed and both phases, output at z = 0:
21 cosmologies x 2 gravity models x 2 phases = 84 runs, about 18 hours (GR runs take about 12 min, f(R) about 14).

Planck 2018 VI (arXiv:1807.06209), Table 2, TT,TE,EE+lowE+lensing: Omega_m = 0.3153 +- 0.0073, H_0 = 67.36 +- 0.54,
ln(10^10 A_s) = 3.044 +- 0.014, Omega_b h^2 = 0.02237, n_s = 0.9649. w0 and wa are centred on LCDM (-1, 0), with the
Planck+BAO+SNe uncertainties 0.080 and 0.29 (Sec. 7.4) as steps; w0waCDM at (-1, 0) is identical to LCDM.
- Neutrinos are massless as in parameterfile.lua, so Omega_m = Omega_cdm + Omega_b; Omega_m is changed through Omega_cdm
  at fixed Omega_b, and H_0 at fixed Omega_m and Omega_b.
- wa = +5 sigma would give w0 + wa > 0, where the dark energy dominates at early times (CLASS cannot compute it), so the
  largest positive step is +3 sigma (w0 + wa = -0.13).
- f(R) on a w0wa background (the w0, wa runs) is COLASolver's f(R) force law on that expansion history.
- More f(R) models: gravity=["F4", "F5", "F6"] (about 10 more hours each); more seeds: seeds=[1001, 1002] (doubles
  the time). More redshifts: add "output_redshifts": [2.0, 1.0, 0.5, 0.0] and "timestep_nsteps":
  steps_per_interval(30, [2.0, 1.0, 0.5, 0.0]) (from runs) to fixed. Finished runs are skipped, so seeds and
  models can be added later.

Writes figures/scans/planck.pdf (15 spectra per page) and figures/scans/planck/P_ij.pdf (one spectrum per file).

    python3 scripts/scan_planck.py [--dry-run] [--plot-only]
"""
import numpy as np

from plot_page import save_pages
from plot_single import save_single
from runs import scan
from spectra import scan_figures

SIGMAS = [-5, -1, 1, 5]
h, Omega_m, Omega_b = 0.6736, 0.3153, round(0.02237 / 0.6736**2, 6)
ln_As = 3.044  # ln(10^10 A_s)

PLANCK = {"cosmology_model": "w0waCDM", "cosmology_w0": -1.0, "cosmology_wa": 0.0, "cosmology_h": h,
          "cosmology_Omegab": Omega_b, "cosmology_OmegaCDM": round(Omega_m - Omega_b, 6),
          "cosmology_As": float(f"{np.exp(ln_As) * 1e-10:.6g}"), "cosmology_ns": 0.9649}

if __name__ == "__main__":  # So that other scans can import PLANCK
    planck = scan("planck", fixed=PLANCK, gravity=["F5"], seeds=[1001], one_at_a_time={
        "cosmology_OmegaCDM": [round(Omega_m + n * 0.0073 - Omega_b, 6) for n in SIGMAS],
        "cosmology_h": [round(h + n * 0.0054, 6) for n in SIGMAS],
        "cosmology_As": [float(f"{np.exp(ln_As + n * 0.014) * 1e-10:.6g}") for n in SIGMAS],
        "cosmology_w0": [round(-1.0 + n * 0.080, 6) for n in SIGMAS],
        "cosmology_wa": [round(n * 0.29, 6) for n in [-5, -1, 1, 3]],
    })
    figures = scan_figures(planck)
    save_pages(figures, "scans/planck.pdf")
    save_single(figures, "scans/planck")
