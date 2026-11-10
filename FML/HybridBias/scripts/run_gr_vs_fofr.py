#!/usr/bin/env python3
"""GR vs f(R) (F5, from parameterfile.lua) at several resolutions (force mesh f, particles n), 30 steps, one seed with
fixed amplitude as a phase-reversed pair, both models from identical initial conditions. Writes
figures/fofr_boost/all_boost_gr_vs_fofr.pdf and figures/fofr_boost/gr_vs_fofr/P_ij.pdf: the boosts f(R)/GR.

    python3 scripts/run_gr_vs_fofr.py [--dry-run] [--plot-only]
"""
from plot_page import save_pages
from plot_single import save_single
from runs import Run, run, summary
from spectra import boosts

runs = []
for f, n in [(512, 512), (1024, 512), (1024, 1024), (1280, 1280)]:
    for phase in ["normal", "reversed"]:
        for gravity in ["GR", "f(R)"]:
            name = f"{'gr' if gravity == 'GR' else 'fofr'}_f{f}_n{n}_s30_fixed_{phase}_seed1234"
            folder = run(f"gr_vs_fofr/{name}", ntasks=128 if n > 1024 else 64, gravity_model=gravity, force_nmesh=f,
                         particle_Npart_1D=n, timestep_nsteps=[30], ic_use_gravity_model_GR=True, ic_fix_amplitude=True,
                         ic_reverse_phases=phase == "reversed", ic_random_seed=1234)
            runs.append(Run(folder, values=(1024, f, n, 30), model="GR" if gravity == "GR" else "F5", seed=1234,
                            reverse=phase == "reversed"))
summary()

figures = [{"title": "f(R) / GR at z = 0.0", "names": ["model", "L", "f", "n", "s"], "ratio": True, "kmin": 0.1,
            "kmax": 1.0, "curves": {("F5", *values): curve for values, curve in boosts(runs, "F5").items()}}]
save_pages(figures, "fofr_boost/all_boost_gr_vs_fofr.pdf")
save_single(figures, "fofr_boost/gr_vs_fofr")
