#!/usr/bin/env python3
"""Smoothing of the linear field the bias weights are built from (pofk_bias_smoothing_scale_hmpc, Gaussian): the Planck
fiducial without smoothing and with R = 0.5, 1, 2, 3 and 4 Mpc/h, in GR and f(R) (F5), one seed and both phases:
6 x 2 x 2 = 24 runs, about 5 hours. The smoothing does not change the simulation, so the matter spectrum P_00 is the
same in all of them; the ratio pages are relative to R = 0.
Writes figures/scans/smoothing.pdf and figures/scans/smoothing/P_ij.pdf.

Convergence at one smoothing scale (commented out below): 1024^3 and 2048^3 particles x force mesh 1024 and 2048 at
R = R_CONVERGENCE, GR only. Without smoothing the quadratic fields depend on the particle grid; with it they should not.

    python3 scripts/scan_smoothing.py [--dry-run] [--plot-only]
"""
from plot_page import save_pages
from plot_single import save_single
from runs import scan
from scan_planck import PLANCK
from spectra import scan_figures

R = "pofk_bias_smoothing_scale_hmpc"
KRANGE = {"kmin": 0.01, "kmax": 1.0}

smoothing = scan("smoothing", fixed=PLANCK, gravity=["F5"], seeds=[1001], one_at_a_time={R: [0.5, 0.75, 1.0, 2.0, 3.0]})
figures = [{**f, **KRANGE} for f in scan_figures(smoothing)]
save_pages(figures, "scans/smoothing.pdf")

# Convergence at one smoothing scale: 4 runs, the 2048^3 ones take several hours each and use 128 tasks
# R_CONVERGENCE = 1.0
# convergence = scan("smoothing_convergence", seeds=[1001], phases=["normal"], ntasks=128,
#                    fixed={**PLANCK, R: R_CONVERGENCE, "gravity_model": "GR", "ic_use_gravity_model_GR": True},
#                    vary={"particle_Npart_1D": [1024, 2048], "force_nmesh": [1024, 2048]})
# figures = [{**f, **KRANGE} for f in scan_figures(convergence)]
# save_pages(figures, "scans/smoothing_convergence.pdf")
# save_single(figures, "scans/smoothing_convergence")
