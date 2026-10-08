#!/usr/bin/env python3
"""100 realisations of parameterfile.lua with different seeds (12300 + i), at 512^3 with 10 steps. Writes
figures/ensemble/all_spectra.pdf and figures/ensemble/P_ij.pdf: their mean and standard deviation.

    python3 scripts/run_ensemble_100.py [--dry-run] [--plot-only]
"""
from plot_page import save_pages
from plot_single import save_single
from runs import run, summary
from spectra import mean_curve

folders = [run(f"ensemble_100seeds/output_seed_{i}", ic_random_seed=12300 + i,
               particle_Npart_1D=512, force_nmesh=512, timestep_nsteps=[10]) for i in range(1, 101)]
summary()

curve = mean_curve(folders)
figures = [{"title": "z = 0.0", "curves": [{**curve, "label": "ensemble_100seeds" + curve["info"]}]}]
save_pages(figures, "ensemble/all_spectra.pdf")
save_single(figures, "ensemble")
