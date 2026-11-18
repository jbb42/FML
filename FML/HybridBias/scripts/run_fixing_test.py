#!/usr/bin/env python3
"""How large is the error from fixing h and A_s at Planck? A direct test without an emulator: 3 random cosmologies from
the production ranges (scripts/run_production.py, its design code with another random seed), each with GR + 3 random
realistic mu models, run twice: h and A_s at Planck, and h and A_s shifted (random within +-3 sigma, the same shift for
all models of a cosmology). The difference of D = (P_mu - P_GR) / sqrt(P_ii^GR P_jj^GR) between the two versions is
the error of fixing them at that point, with all interactions. Same seed everywhere; no control-variate prediction
(not needed for a difference at the same seed). 24 runs in results/fixing_test/c<i>/<planck|shifted>/<model>/.

    python3 scripts/run_fixing_test.py --design-only     (on a machine with scipy, once)
    python3 scripts/run_fixing_test.py                   (on the compute nodes; several at once share the queue)
"""
import os
import sys
import time

import numpy as np

import run_production as rp
import runs
from scan_planck import PLANCK

NAME = "fixing_test"
BASE = os.path.join(runs.ROOT, "results", NAME)
DESIGN = os.path.join(BASE, "design.txt")
SIGMA_H, SIGMA_LNAS = 0.0054, 0.014  # Planck 2018


def write_design():
    rng = np.random.default_rng(20261119)
    rows = rp.make_design(3, rng)
    shifts = {c: (rng.uniform(-3, 3), rng.uniform(-3, 3)) for c in range(3)}
    os.makedirs(BASE, exist_ok=True)
    with open(DESIGN, "w") as f:
        f.write("# cosmology model Omega_cdm w0 wa A n mu0_equiv h_shifted As_shifted (planck: h 0.6736, As "
                f"{PLANCK['cosmology_As']}); shifts in sigma: {shifts}\n")
        for r in rows:
            xh, xa = shifts[r[0]]
            h = round(PLANCK["cosmology_h"] + xh * SIGMA_H, 6)
            As = float(f"{PLANCK['cosmology_As'] * np.exp(xa * SIGMA_LNAS):.6g}")
            f.write(f"{r[0]:4d} {r[1]:4s} {r[2]:9.6f} {r[3]:9.6f} {r[4]:9.6f} {r[5]:9.6f} {r[6]:8.6f} {r[7]:9.6f} "
                    f"{h:8.6f} {As:.6g}\n")
    print("Wrote", DESIGN)


def jobs():
    out = []
    for line in open(DESIGN):
        if line.startswith("#"):
            continue
        c, model, oc, w0, wa, A, n, _, h, As = line.split()
        for version in ("planck", "shifted"):
            out.append((int(c), model, version, float(oc), float(w0), float(wa), float(A), float(n),
                        float(h) if version == "shifted" else PLANCK["cosmology_h"],
                        float(As) if version == "shifted" else PLANCK["cosmology_As"]))
    return out


def run_job(job):
    c, model, version, oc, w0, wa, A, n, h, As = job
    cosmo = {**PLANCK, "cosmology_OmegaCDM": oc, "cosmology_w0": w0, "cosmology_wa": wa, "cosmology_h": h,
             "cosmology_As": As}
    input_file = os.path.join(BASE, "input", f"pofk_c{c:03d}_{version}.txt")
    runs.write_input_power({**runs.PROJECT, **cosmo}, input_file)
    gravity = {"gravity_model": "GR"} if model == "GR" else {
        "gravity_model": "mu", "gravity_model_mu_parametrization": "powerlaw", "gravity_model_mu_powerlaw_params": [A, n]}
    params = {**cosmo, **gravity, "ic_type_of_input": "powerspectrum", "ic_input_filename": input_file,
              "ic_random_seed": rp.SEED, "ic_fix_amplitude": True, "ic_reverse_phases": False,
              "output_redshifts": [1.0, 0.0], "timestep_nsteps": runs.steps_per_interval(30, [1.0, 0.0]),
              "hybridbias_zeldovich": False, "hybridbias_zeldovich_field": False}
    return runs.run(f"{NAME}/c{c:03d}/{version}/{model}", ntasks=64, **params)


if __name__ == "__main__":
    if "--design-only" in sys.argv:
        write_design()
        sys.exit()
    for job in jobs():
        folder = os.path.join(BASE, f"c{job[0]:03d}", job[2], job[1])
        if runs.is_finished(folder) or not rp.claim(folder):
            continue
        start = time.time()
        run_job(job)
        with open(os.path.join(BASE, "progress.txt"), "a") as f:
            f.write(f"{time.strftime('%F %T')} c{job[0]:03d}/{job[2]}/{job[1]} {time.time() - start:.0f} s "
                    f"{'ok' if runs.is_finished(folder) else 'FAILED'} {os.uname().nodename}\n")
