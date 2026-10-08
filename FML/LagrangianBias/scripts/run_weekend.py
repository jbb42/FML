#!/usr/bin/env python3
"""GR vs f(R) (F4, F5, F6) with f = n = 1024, fixed amplitudes and phase-reversed pairs, identical initial conditions
for every model, outputs at z = 2, 1, 0.5, 0. One simulation at a time, in this order:
    ensemble: 10 seeds with 30 steps in 1024 Mpc/h, each also in 2048 Mpc/h (odd seeds) or 512 Mpc/h (even seeds)
    steps:    20, 30 and 40 steps in 1024 Mpc/h, first seed only
Writes figures/fofr_boost/all_boost_weekend.pdf and figures/fofr_boost/weekend/P_ij.pdf: the boosts f(R)/GR, varying
the model, box size L or steps s.

    python3 scripts/run_weekend.py [ensemble] [steps] [--dry-run] [--plot-only]    (default: both)
"""
import sys

from plot_page import save_pages
from plot_single import save_single
from runs import Run, gravity_parameters, run, steps_per_interval, summary
from spectra import boosts

MODELS = ["GR", "F4", "F5", "F6"]
REDSHIFTS = [2.0, 1.0, 0.5, 0.0]
SEEDS = range(1001, 1011)
runs = []


def run_models(box, nsteps, seed):
    """Every model, with normal and reversed phases."""
    for model in MODELS:
        for phase in ["normal", "reversed"]:
            name = f"{'gr' if model == 'GR' else 'fofr' + model}_L{box}_f1024_n1024_s{nsteps}_fixed_{phase}_seed{seed}"
            folder = run(f"weekend/{name}", simulation_boxsize=box, force_nmesh=1024, particle_Npart_1D=1024,
                         **gravity_parameters(model), ic_use_gravity_model_GR=True, ic_fix_amplitude=True,
                         ic_reverse_phases=phase == "reversed", ic_random_seed=seed, output_redshifts=REDSHIFTS,
                         timestep_nsteps=steps_per_interval(nsteps, REDSHIFTS))
            runs.append(Run(folder, values=(box, 1024, 1024, nsteps), model=model, seed=seed, reverse=phase == "reversed"))


campaigns = [arg for arg in sys.argv[1:] if not arg.startswith("--")] or ["ensemble", "steps"]
if "ensemble" in campaigns:
    for i, seed in enumerate(SEEDS):
        run_models(1024, 30, seed)
        run_models(2048 if i % 2 == 0 else 512, 30, seed)
if "steps" in campaigns:
    for nsteps in [20, 30, 40]:
        run_models(1024, nsteps, SEEDS[0])
summary()

figures = [{"title": f"f(R) / GR at z = {z}", "names": ["model", "L", "f", "n", "s"], "ratio": True, "kmin": 0.1,
            "kmax": 1.0, "curves": {(model, *values): curve for model in MODELS[1:]
                                    for values, curve in boosts(runs, model, z).items()}}
           for z in sorted(REDSHIFTS)]
save_pages(figures, "fofr_boost/all_boost_weekend.pdf")
save_single(figures, "fofr_boost/weekend")
