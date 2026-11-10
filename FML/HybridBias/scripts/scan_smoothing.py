#!/usr/bin/env python3
"""Smoothing of the linear field the bias weights are built from (pofk_bias_smoothing_scale_hmpc, Gaussian): the Planck
fiducial without smoothing and with R = 0.5, 0.75, 1, 2 and 3 Mpc/h, in GR and f(R) (F5), one seed and both phases:
6 x 2 x 2 = 24 runs, about 3.5 hours (done; commented out below). The smoothing does not change the simulation, so the matter spectrum P_00 is the
same in all of them; the ratio pages are relative to R = 0.
Writes figures/scans/smoothing.pdf and figures/scans/smoothing/P_ij.pdf.

Convergence at one smoothing scale: 512^3, 1024^3 and 1280^3 particles x force mesh 512, 1024 and 1280 at
R = R_CONVERGENCE, in GR and F5 (for the boosts), one seed and the normal phase: 9 x 2 = 18 runs. The ratio pages are
relative to 1280/1280. Without smoothing the quadratic fields depend on the particle grid; with it they should not (at
512^3 the grid, 2 Mpc/h, still cuts off before R = 1 Mpc/h does).
Writes figures/scans/smoothing_convergence.pdf and figures/scans/smoothing_convergence/P_ij.pdf.

    python3 scripts/scan_smoothing.py [--dry-run] [--plot-only]
"""
from plot_page import save_pages
from plot_single import save_single
from runs import scan
from scan_planck import PLANCK
from spectra import scan_figures

R = "pofk_bias_smoothing_scale_hmpc"
KRANGE = {"kmin": 0.01, "kmax": 1.0}


def convergence_options(figure):
    """k range for every page, and on the pages relative to the best resolution the same +-5% y range in every panel,
    with a grey +-1% band, so the size of a deviation can be read off directly."""
    same_range = {"ylim": (0.95, 1.05), "band": 0.01} if figure["kind"] == "relative" else {}
    return {**figure, **KRANGE, **same_range}


#smoothing = scan("smoothing", fixed=PLANCK, gravity=["F5"], seeds=[1001], one_at_a_time={R: [0.5, 0.75, 1.0, 2.0, 3.0]})
#figures = [{**f, **KRANGE} for f in scan_figures(smoothing)]
#save_pages(figures, "scans/smoothing.pdf")

# Convergence at one smoothing scale: 9 setups in GR and F5, 18 runs. If a 1280^3 run runs out of memory, use more tasks
R_CONVERGENCE = 1.0
convergence = scan("smoothing_convergence", fixed={**PLANCK, R: R_CONVERGENCE}, gravity=["F5"], seeds=[1001], ntasks=64,
                   vary={"particle_Npart_1D": [1280, 1024, 512], "force_nmesh": [1280, 1024, 512]})  # Largest first
convergence.reference = (1280, 1280)  # Compare with the highest resolution (vary has no fiducial)
figures = [convergence_options(f) for f in scan_figures(convergence)]
save_pages(figures, "scans/smoothing_convergence.pdf")
save_single(figures, "scans/smoothing_convergence")
