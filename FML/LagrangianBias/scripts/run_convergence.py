#!/usr/bin/env python3
"""Convergence runs in a 1024 Mpc/h box: force mesh f x particles n x time steps s, 10 seeds (versions) each.
Results in results/convergence_runs_v2_fixedseeds/f<f>_n<n>_s<s>_v<version>/. Plot with plot_convergence.py.

Usage: python3 scripts/run_convergence.py [--dry-run]
"""
from simulations import run, summary

# (force mesh, particles, time steps): the 512/1024 grid plus the highest resolution
SETUPS = [(f, n, s) for f in [512, 1024] for n in [512, 1024] for s in [10, 20, 40]] + [(1280, 1280, 40)]

for f, n, s in SETUPS:
    for version in range(1, 11):
        run(f"convergence_runs_v2_fixedseeds/f{f}_n{n}_s{s}_v{version}", ntasks=128 if n > 1024 else 64,
            force_nmesh=f, particle_Npart_1D=n, timestep_nsteps=[s], ic_random_seed=1230 + version)
summary()
