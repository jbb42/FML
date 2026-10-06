#!/usr/bin/env python3
"""The Planck 2018 cosmology, and Omega_m, H_0, A_s, w0 and wa each changed on its own by -5, -1, +1 and +5 sigma while
the others keep their Planck values; one seed, both phases: 1 + 4 x 4 + 3 = 20 cosmologies, 40 runs.

Planck 2018 VI (arXiv:1807.06209), Table 2, TT,TE,EE+lowE+lensing: Omega_m = 0.3153 +- 0.0073, H_0 = 67.36 +- 0.54,
ln(10^10 A_s) = 3.044 +- 0.014, Omega_b h^2 = 0.02237, n_s = 0.9649. w0 and wa are centred on LCDM (-1, 0), with the
Planck+BAO+SNe uncertainties 0.080 and 0.29 (Sec. 7.4) as steps; w0waCDM at (-1, 0) is identical to LCDM.
- Neutrinos are massless as in parameterfile.lua, so Omega_m = Omega_cdm + Omega_b; Omega_m is changed through Omega_cdm
  at fixed Omega_b, and H_0 at fixed Omega_m and Omega_b.
- wa = +5 sigma is left out: with w0 + wa > 0 the dark energy dominates at early times, and CLASS cannot compute it.

Run as python3 scripts/scan_planck.py [--dry-run] [--plot-only]; writes figures/scans/planck.pdf.
"""
import numpy as np

from scan import scan

SIGMAS = [-5, -1, 1, 5]
h, Omega_m, Omega_b = 0.6736, 0.3153, round(0.02237 / 0.6736**2, 6)
ln_As = 3.044  # ln(10^10 A_s)

PLANCK = {"cosmology_model": "w0waCDM", "cosmology_w0": -1.0, "cosmology_wa": 0.0, "cosmology_h": h,
          "cosmology_Omegab": Omega_b, "cosmology_OmegaCDM": round(Omega_m - Omega_b, 6),
          "cosmology_As": float(f"{np.exp(ln_As) * 1e-10:.6g}"), "cosmology_ns": 0.9649}

scan("planck", one_at_a_time={
    "cosmology_OmegaCDM": [round(Omega_m + n * 0.0073 - Omega_b, 6) for n in SIGMAS],
    "cosmology_h": [round(h + n * 0.0054, 6) for n in SIGMAS],
    "cosmology_As": [float(f"{np.exp(ln_As + n * 0.014) * 1e-10:.6g}") for n in SIGMAS],
    "cosmology_w0": [round(-1.0 + n * 0.080, 6) for n in SIGMAS],
    "cosmology_wa": [round(n * 0.29, 6) for n in SIGMAS if n != 5],
}, fixed=PLANCK, seeds=[1001])
