#!/usr/bin/env python3
"""Force mesh f x box size L x particles n, each in {256, 512, 1024}, 10 steps. Writes
figures/resolution_scan/all_spectra.pdf and figures/resolution_scan/P_ij.pdf.

    python3 scripts/run_resolution_scan.py [--dry-run] [--plot-only]
"""
from plot_page import save_pages
from plot_single import save_single
from runs import run, summary
from spectra import mean_curve

spectra = {}
for f in [256, 512, 1024]:
    for box in [256, 512, 1024]:
        for n in [256, 512, 1024]:
            folder = run(f"resolution_scan/f{f}_L{box}_n{n}", force_nmesh=f, simulation_boxsize=box, particle_Npart_1D=n,
                         timestep_nsteps=[10])
            spectra[(f, box, n)] = mean_curve([folder])
summary()

figures = [{"title": "z = 0.0", "curves": {setup: curve for setup, curve in spectra.items() if curve},
            "names": ["f", "L", "n"]}]
save_pages(figures, "resolution_scan/all_spectra.pdf")
save_single(figures, "resolution_scan")
