#!/usr/bin/env python3
"""100 realisations of parameterfile.lua with different seeds (12300 + i), in results/ensemble_100seeds/output_seed_i/.
Plot with: python3 scripts/plot_runs.py ensemble/all_spectra.pdf 'results/ensemble_100seeds/*'

Usage: python3 scripts/run_ensemble_100.py [--dry-run]
"""
from simulations import run, summary

for i in range(1, 101):
    run(f"ensemble_100seeds/output_seed_{i}", ic_random_seed=12300 + i)
summary()
