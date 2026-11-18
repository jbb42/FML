#!/usr/bin/env python3
"""Production runs: hybrid bias spectra for GR and mu(a) = 1 + A a^n in many cosmologies, for an emulator of
D_ij = (P_mu - P_GR) / sqrt(P_ii^GR P_jj^GR).

Design (nested, written once to results/<name>/design.txt and never changed afterwards):
  - N cosmologies, an optimised Latin hypercube (random-cd) in (Omega_cdm, w0, wa) within +-5 sigma of Planck 2018 (wa: -5 / +3 sigma), with
    w0 + wa <= -0.13 (CLASS fails as w0 + wa -> 0; -0.13 is the largest value tested): wa is drawn uniformly between
    its lower limit and min(its upper limit, -0.13 - w0), so the hypercube property holds in (Omega_cdm, w0, u_wa).
    h, A_s, n_s, Omega_b fixed at Planck (results/mu_cosmology: they change D by < 0.5% at k < 0.3, ~1% at k ~ 1 for
    1 sigma).
  - M = 3 mu models per cosmology: one Latin hypercube of 3N points in (n, mu0_equiv), shuffled and dealt out 3 per
    cosmology. n in [1, 4]; mu0_equiv (the DESI mu0 with the same growth, scripts/mu_prior.py) uniform between the 3
    sigma limits [-1.08, 0.71], further limited at each n so that A stays in [-0.5, 0.5]; A from mu0_equiv by root
    finding.
  - Every cosmology also has a GR run; GR and mu share the initial conditions (seed 1001, fixed amplitude, normal
    phase) so that the noise cancels in D.
Settings: scripts/runs.py PROJECT (1024^3 particles on 2048^3 in 1024 Mpc/h, 30 steps, single precision, R = 1 Mpc/h,
kcut 1, hybridbias_nmesh 512, Zel'dovich control variate). particle_allocation_factor: the default 1.25; a run that
overfills a buffer ("exceeds allocated particle capacity") is rerun once with 1.5.

Jobs run in order (cosmology 0: GR, mu 0-2; cosmology 1: ...). Several nodes can run the same command at once: a node
claims a job by creating <folder>/.claim (O_EXCL), so no job runs twice. Finished jobs are skipped; a job that failed
keeps its claim (rerun it by deleting the .claim file). Progress in results/<name>/progress.txt.

    python3 scripts/run_production.py NAME [N] [--design-only]       e.g. production 150; smoke test: production_smoke 1
Write the design first on a machine with scipy (--design-only); the compute nodes only read it.
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# The single-precision native build of FML 2.0 (build_euclid22.sh; regression-tested against the old nbody_float on
# 2026-11-18: Zel'dovich prediction bit-identical, spectra to <= 4e-5 of sqrt(P_ii P_jj), 5% faster)
os.environ.setdefault("NBODY", os.path.join(os.path.dirname(os.path.dirname(HERE)), "build_euclid22_float_v2", "nbody"))
import runs  # noqa: E402
from scan_planck import PLANCK  # noqa: E402

REDSHIFTS = [1.0, 0.0]  # To be decided with Hans; each extra output costs ~8 min (the control-variate prediction)
SEED, M = 1001, 3
SIGMA = {"cosmology_OmegaCDM": 0.0073, "cosmology_w0": 0.080, "cosmology_wa": 0.29}
RANGE = {"cosmology_OmegaCDM": (PLANCK["cosmology_OmegaCDM"] - 5 * 0.0073, PLANCK["cosmology_OmegaCDM"] + 5 * 0.0073),
         "cosmology_w0": (-1.4, -0.6), "cosmology_wa": (-1.45, 0.87)}
W0WA_MAX = -0.13
N_RANGE, MU0_RANGE, A_MAX = (1.0, 4.0), (-1.08, 0.71), 0.5


def lhs(n, d, rng):
    """A Latin hypercube of n points in [0, 1]^d (one point in each of the n strata of every dimension), optimised
    with 'random-cd' (coordinate swaps that lower the centred discrepancy, keeping the Latin hypercube): 6x (n = 150,
    d = 3) to 13x (n = 450, d = 2) more uniform than a plain random one. scipy.stats.qmc.LatinHypercube; FML 2.0 has
    the same algorithm in C++ (LatinHyperCube/LatinHyperCube.h). Needs scipy, like the rest of make_design."""
    from scipy.stats import qmc
    return qmc.LatinHypercube(d, optimization="random-cd", seed=int(rng.integers(2**31))).random(n)


def make_design(N, rng):
    """Needs scipy (via mu_prior), so it is only called where the design is written, not on the compute nodes."""
    import mu_prior
    from scipy.optimize import brentq

    # mu0_equiv <-> growth change, and A at given n (scripts/mu_prior.py)
    mu0s = np.linspace(-1.5, 1.5, 61)
    growth = np.array([mu_prior.summary(mu_prior.mu_de(m))[0] for m in mu0s])

    def mu0_equiv(A, n):
        return float(np.interp(mu_prior.summary(mu_prior.powerlaw(A, n))[0], growth, mu0s))

    def A_of(mu0, n):
        if mu0 == 0:
            return 0.0
        return brentq(lambda A: mu0_equiv(A, n) - mu0, -A_MAX if mu0 < 0 else 0.0, 0.0 if mu0 < 0 else A_MAX,
                      xtol=1e-6)

    rows = []
    u = lhs(N, 3, rng)
    for c in range(N):
        oc = RANGE["cosmology_OmegaCDM"][0] + u[c, 0] * np.diff(RANGE["cosmology_OmegaCDM"])[0]
        w0 = RANGE["cosmology_w0"][0] + u[c, 1] * np.diff(RANGE["cosmology_w0"])[0]
        wa_lo, wa_hi = RANGE["cosmology_wa"][0], min(RANGE["cosmology_wa"][1], W0WA_MAX - w0)
        wa = wa_lo + u[c, 2] * (wa_hi - wa_lo)
        rows.append((c, "GR", round(oc, 6), round(w0, 6), round(wa, 6), 0.0, 0.0, 0.0))
    v = lhs(M * N, 2, rng)
    v = v[rng.permutation(M * N)]
    for i, (un, um) in enumerate(v):
        c = i // M
        n = N_RANGE[0] + un * (N_RANGE[1] - N_RANGE[0])
        lo, hi = max(MU0_RANGE[0], mu0_equiv(-A_MAX, n)), min(MU0_RANGE[1], mu0_equiv(A_MAX, n))
        mu0 = lo + um * (hi - lo)
        A = A_of(mu0, n)
        _, oc, w0, wa = rows[c][1:5]
        rows.append((c, f"mu{i % M}", oc, w0, wa, round(A, 6), round(n, 6), round(mu0, 6)))
    return sorted(rows, key=lambda r: (r[0], r[1] != "GR", r[1]))


def read_design(path):
    rows = []
    for line in open(path):
        if line.startswith("#"):
            continue
        c, model, *vals = line.split()
        rows.append((int(c), model, *map(float, vals)))
    return rows


def finished(folder):
    return all(os.path.exists(os.path.join(folder, f"snapshot_TestSim_z{z:.3f}", f"hybridbias_cv_z{z:.3f}.txt"))
               for z in REDSHIFTS)


def claim(folder):
    os.makedirs(folder, exist_ok=True)
    try:
        os.close(os.open(os.path.join(folder, ".claim"), os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        return True
    except FileExistsError:
        return False


def run_job(name, row, base):
    c, model, oc, w0, wa, A, n, _ = row
    cosmo = {**PLANCK, "cosmology_OmegaCDM": oc, "cosmology_w0": w0, "cosmology_wa": wa}
    input_file = os.path.join(base, "input", f"pofk_c{c:03d}.txt")
    runs.write_input_power({**runs.PROJECT, **cosmo}, input_file)
    gravity = {"gravity_model": "GR"} if model == "GR" else {
        "gravity_model": "mu", "gravity_model_mu_parametrization": "powerlaw", "gravity_model_mu_powerlaw_params": [A, n]}
    params = {**cosmo, **gravity, "ic_type_of_input": "powerspectrum", "ic_input_filename": input_file,
              "ic_random_seed": SEED, "ic_fix_amplitude": True, "ic_reverse_phases": False,
              "output_redshifts": REDSHIFTS, "timestep_nsteps": runs.steps_per_interval(30, REDSHIFTS)}
    job = f"{name}/c{c:03d}/{model}"
    folder = runs.run(job, ntasks=64, **params)
    if not finished(folder) and "exceeds allocated particle capacity" in open(os.path.join(folder, "log.txt")).read():
        print(time.strftime("%T"), job, "overfilled a particle buffer: rerun with particle_allocation_factor 1.5",
              flush=True)
        folder = runs.run(job, ntasks=64, **params, particle_allocation_factor=1.5)
    return finished(folder)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    name = args[0]
    base = os.path.join(runs.ROOT, "results", name)
    design = os.path.join(base, "design.txt")
    if not os.path.exists(design):
        N = int(args[1])
        rows = make_design(N, np.random.default_rng(20261118))
        os.makedirs(base, exist_ok=True)
        with open(design, "w") as f:
            f.write(f"# {N} cosmologies x (GR + {M} mu), seed {SEED}, redshifts {REDSHIFTS}. Written once; do not edit.\n"
                    "# cosmology model Omega_cdm w0 wa A n mu0_equiv (h, A_s, n_s, Omega_b: Planck 2018, scan_planck.py)\n")
            for r in rows:
                f.write(f"{r[0]:4d} {r[1]:4s} {r[2]:9.6f} {r[3]:9.6f} {r[4]:9.6f} {r[5]:9.6f} {r[6]:8.6f} {r[7]:9.6f}\n")
        print("Wrote", design)
    if "--design-only" in sys.argv:
        sys.exit()
    for row in read_design(design):
        folder = os.path.join(base, f"c{row[0]:03d}", row[1])
        if finished(folder) or not claim(folder):
            continue
        start = time.time()
        ok = run_job(name, row, base)
        with open(os.path.join(base, "progress.txt"), "a") as f:
            f.write(f"{time.strftime('%F %T')} c{row[0]:03d}/{row[1]} {time.time() - start:.0f} s "
                    f"{'ok' if ok else 'FAILED'} {os.uname().nodename}\n")
