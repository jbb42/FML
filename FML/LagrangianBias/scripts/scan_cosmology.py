#!/usr/bin/env python3
"""The fiducial cosmology of parameterfile.lua, and Omega_m and A_s each changed on its own by -50%, -10%, +10% and
+50%, with two seeds and both phases: 1 + 4 + 4 = 9 cosmologies, 36 runs. Omega_m is changed through Omega_cdm at
fixed Omega_b. Add parameters the same way, e.g. "cosmology_ns": around_fiducial("cosmology_ns", FRACTIONS).
Writes figures/scans/cosmology.pdf and figures/scans/cosmology/P_ij.pdf.

    python3 scripts/scan_cosmology.py [--dry-run] [--plot-only]
"""
from plot_page import save_pages
from plot_single import save_single
from runs import around_fiducial, fiducial, scan
from spectra import scan_figures

FRACTIONS = [-0.5, -0.1, 0.1, 0.5]
Omega_b = fiducial("cosmology_Omegab")
Omega_m = fiducial("cosmology_OmegaCDM") + Omega_b

cosmology = scan("cosmology", seeds=[1001, 1002], one_at_a_time={
    "cosmology_OmegaCDM": [round(Omega_m * (1 + fraction) - Omega_b, 6) for fraction in FRACTIONS],
    "cosmology_As": around_fiducial("cosmology_As", FRACTIONS),
})
figures = scan_figures(cosmology)
save_pages(figures, "scans/cosmology.pdf")
save_single(figures, "scans/cosmology")
