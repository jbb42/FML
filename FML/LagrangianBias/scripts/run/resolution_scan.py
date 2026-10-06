#!/usr/bin/env python3
"""Force mesh f x box size L x particles n, each in {256, 512, 1024}, in results/resolution_scan/f<f>_L<L>_n<n>/.
Plot with plot_runs.py, e.g. python3 scripts/plot/plot_runs.py resolution_scan/all_spectra.pdf results/resolution_scan/*

Usage: python3 scripts/run/resolution_scan.py [--dry-run]
"""
from simulations import run, summary

for f in [256, 512, 1024]:
    for box in [256, 512, 1024]:
        for n in [256, 512, 1024]:
            run(f"resolution_scan/f{f}_L{box}_n{n}", force_nmesh=f, simulation_boxsize=box, particle_Npart_1D=n)
summary()
