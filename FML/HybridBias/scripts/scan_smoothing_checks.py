#!/usr/bin/env python3
"""Follow-up of the R = 1 Mpc/h convergence scan in scan_smoothing.py, about 10 hours, run in the order A, C, B.

A. Particle count without smoothing and with more: R = 0 and R = 2 Mpc/h at 1280^3, 1024^3 and 512^3 particles (force
   mesh 1280), GR,
   one seed and both phases: 2 x 3 x 2 = 12 runs, about 2 hours. With R = 1 (scan_smoothing.py) this shows how the
   smoothing removes the dependence on the particle grid, and whether the large-scale nabla^2 delta_L spectra
   (P_04, P_14, P_44 at k < 0.03), which differ between 1024^3 and 1280^3 at R = 1, converge with more smoothing.
B. Seeds 1002-1004 at R = 1 for 1280^3 and 1024^3 particles with force mesh 1280, GR and F5, both phases:
   3 x 2 x 2 x 2 = 24 runs, about 6 hours. Seed 1001 is reused from scan_smoothing.py. Shows whether the nabla^2 delta_L
   difference between 1024^3 and 1280^3 is larger than the scatter between seeds (different particle counts are
   different realizations), and gives the boosts error bars.
C. The same seeds for the production setup, 1024^3 particles and force mesh 1024: 3 x 2 x 2 = 12 runs, about 1.6 hours.
   Gives error bars on how well 1024/1024 agrees with 1280/1280, for the spectra and the boosts.

Writes figures/scans/smoothing_R0_convergence.pdf, smoothing_R2_convergence.pdf (A) and smoothing_seeds.pdf (B and C,
mean over the 4 seeds with their scatter, boosts with their standard error), all relative to 1280^3 / 1280.

    python3 scripts/scan_smoothing_checks.py [--dry-run] [--plot-only]
"""
from plot_page import save_pages
from runs import Scan, scan
from scan_planck import PLANCK
from spectra import scan_figures

R = "pofk_bias_smoothing_scale_hmpc"
KRANGE = {"kmin": 0.01, "kmax": 1.0}


def convergence_options(figure):
    """k range for every page, and on the pages relative to the best resolution the same +-5% y range in every panel,
    with a grey +-1% band, so the size of a deviation can be read off directly."""
    same_range = {"ylim": (0.95, 1.05), "band": 0.01} if figure["kind"] == "relative" else {}
    return {**figure, **KRANGE, **same_range}


GR ={"gravity_model": "GR", "ic_use_gravity_model_GR": True}
PARTICLES = {"particle_Npart_1D": [1280, 1024, 512], "force_nmesh": [1280]}  # Largest first (memory)

# A. R = 0 and R = 2 (R = 1 is in results/scans/smoothing_convergence). R is fixed in each scan, so it is not in the
# folder names, and each R needs its own scan name
for r in [0.0, 2.0]:
    particles = scan(f"smoothing_R{r:g}_convergence", fixed={**PLANCK, **GR, R: r}, seeds=[1001], vary=PARTICLES)
    particles.reference = (1280, 1280)
    save_pages([convergence_options(f) for f in scan_figures(particles)], f"scans/smoothing_R{r:g}_convergence.pdf")

# C and B. More seeds at R = 1, in the folders of the convergence scan, so seed 1001 is reused. The two are one scan
# in the figures, as vary can only give all combinations
SEEDS = {"fixed": {**PLANCK, R: 1.0}, "gravity": ["F5"], "seeds": [1001, 1002, 1003, 1004]}
production = scan("smoothing_convergence", **SEEDS, vary={"particle_Npart_1D": [1024], "force_nmesh": [1024]})
seeds = scan("smoothing_convergence", **SEEDS, vary={"particle_Npart_1D": [1280, 1024], "force_nmesh": [1280]})
combined = Scan("smoothing_convergence, 4 seeds", seeds.parameters, (1280, 1280), seeds.models,
                seeds.runs + production.runs)
save_pages([convergence_options(f) for f in scan_figures(combined)], "scans/smoothing_seeds.pdf")
