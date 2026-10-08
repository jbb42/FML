#!/usr/bin/env python3
"""Test of pofk_nmesh = 1024 (instead of 512): the Planck fiducial GR run with the normal phase, otherwise identical to
results/scans/planck/fiducial_GR_seed1001_normal, which is the 512 reference. Compare the two with
    python3 scripts/plot_page.py scans/planck_mesh_vs_512.pdf results/scans/planck/fiducial_GR_seed1001_normal \
        results/scans/planck_mesh/pofk_nmesh1024_seed1001_normal --ratio

    python3 scripts/scan_planck_mesh.py [--dry-run] [--plot-only]
"""
from plot_page import save_pages
from runs import scan
from scan_planck import PLANCK
from spectra import scan_figures

mesh = scan("planck_mesh", vary={"pofk_nmesh": [1024]}, seeds=[1001], phases=["normal"],
            fixed={**PLANCK, "gravity_model": "GR", "ic_use_gravity_model_GR": True})
save_pages(scan_figures(mesh), "scans/planck_mesh.pdf")
