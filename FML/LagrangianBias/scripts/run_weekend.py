#!/usr/bin/env python3
"""Weekend campaign: GR vs f(R) (F4, F5, F6) with f = n = 1024, fixed amplitudes and phase-reversed pairs, identical
initial conditions for every model, outputs at z = 2, 1, 0.5, 0. One simulation at a time, in this order:
    ensemble: for 10 seeds, 30 steps in 1024 Mpc/h, plus alternately 2048 Mpc/h (low k, first seed) and 512 Mpc/h
              (high k), for stitching the boxes together
    steps:    20, 30 and 40 steps in 1024 Mpc/h, first seed only
Results in results/weekend/. Plot with plot_fofr_boost.py weekend. Follow the running simulation with
tail -F results/weekend/current_log.txt

Usage: python3 scripts/run_weekend.py [ensemble] [steps] [--dry-run]   (default: both)
"""
import os
import sys

from simulations import run, steps_per_interval, summary

os.environ["OMP_NUM_THREADS"] = "2"     # 64 MPI tasks x 2 threads fill one node
os.environ["I_MPI_PIN_DOMAIN"] = "omp"  # Intel MPI: give each task its own cores
REDSHIFTS = [2.0, 1.0, 0.5, 0.0]
SEEDS = range(1001, 1011)


def run_models(box, nsteps, seed):
    """GR, F4, F5 and F6, each with normal and reversed phases."""
    for model in ["GR", "F4", "F5", "F6"]:
        for reverse in [False, True]:
            if model == "GR":
                gravity = {"gravity_model": "GR"}
            else:
                gravity = {"gravity_model": "f(R)", "gravity_model_fofr_fofr0": 10.0 ** -int(model[1:])}
            name = f"{'gr' if model == 'GR' else 'fofr' + model}_L{box}_f1024_n1024_s{nsteps}_fixed_" \
                   f"{'reversed' if reverse else 'normal'}_seed{seed}"
            run(f"weekend/{name}", simulation_boxsize=box, force_nmesh=1024, particle_Npart_1D=1024,
                ic_use_gravity_model_GR=True, ic_fix_amplitude=True, ic_reverse_phases=reverse, ic_random_seed=seed,
                output_redshifts=REDSHIFTS, timestep_nsteps=steps_per_interval(nsteps, REDSHIFTS), **gravity)


campaigns = [arg for arg in sys.argv[1:] if not arg.startswith("--")] or ["ensemble", "steps"]
if "ensemble" in campaigns:
    for i, seed in enumerate(SEEDS):
        run_models(1024, 30, seed)
        run_models(2048 if i % 2 == 0 else 512, 30, seed)
if "steps" in campaigns:
    for nsteps in [20, 30, 40]:
        run_models(1024, nsteps, SEEDS[0])
summary()
