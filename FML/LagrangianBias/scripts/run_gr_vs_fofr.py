#!/usr/bin/env python3
"""GR vs f(R) at several resolutions (force mesh f, particles n), 30 steps, one seed with fixed amplitude as a
phase-reversed pair. Both gravity models start from identical initial conditions (ic_use_gravity_model_GR = true).
The f(R) parameters come from parameterfile.lua. Results in results/gr_vs_fofr/. Plot with plot_fofr_boost.py.

Usage: python3 scripts/run_gr_vs_fofr.py [--dry-run]
"""
from simulations import run, summary

for f, n in [(512, 512), (1024, 512), (1024, 1024), (1280, 1280)]:
    for reverse in [False, True]:
        for gravity in ["GR", "f(R)"]:
            name = f"{'gr' if gravity == 'GR' else 'fofr'}_f{f}_n{n}_s30_fixed_{'reversed' if reverse else 'normal'}_seed1234"
            run(f"gr_vs_fofr/{name}", ntasks=128 if n > 1024 else 64, gravity_model=gravity, force_nmesh=f,
                particle_Npart_1D=n, timestep_nsteps=[30], ic_use_gravity_model_GR=True, ic_fix_amplitude=True,
                ic_reverse_phases=reverse, ic_random_seed=1234)
summary()
