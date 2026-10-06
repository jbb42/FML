#!/usr/bin/env python3
"""The fiducial cosmology of parameterfile.lua, and each parameter below changed on its own by -50%, -10%, +10% and +50%
while the others stay fiducial, with two seeds and both phases: 1 + 4 + 4 = 9 cosmologies, 36 runs.
Omega_m is changed through Omega_cdm at fixed Omega_b, so that Omega_m = Omega_cdm + Omega_b changes by these fractions.
Add more parameters the same way, e.g. "cosmology_ns": around_fiducial("cosmology_ns", FRACTIONS).

Run as python3 scripts/scan_cosmology.py [--dry-run] [--plot-only]; writes figures/scans/cosmology.pdf.
"""
from scan import around_fiducial, fiducial, scan

FRACTIONS = [-0.5, -0.1, 0.1, 0.5]
Omega_b = fiducial("cosmology_Omegab")
Omega_m = fiducial("cosmology_OmegaCDM") + Omega_b

scan("cosmology", one_at_a_time={
    "cosmology_OmegaCDM": [round(Omega_m * (1 + fraction) - Omega_b, 6) for fraction in FRACTIONS],
    "cosmology_As": around_fiducial("cosmology_As", FRACTIONS),
}, seeds=[1001, 1002])
