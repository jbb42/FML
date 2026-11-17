"""How much does the cosmology change the mu/GR boost, for different mu(a) = 1 + A a^n? Full scale, one run at a time
per node (~320 GB, ~22 min), the jobs in priority order and split over nodes:
    python3 run_all.py PART NPARTS [--dry-run]     (e.g. 0 2 on euclid24 and 1 2 on euclid25)

mu models: A0.3_n3 (fiducial) and realistic corners: A in {-0.35, 0.2} at n = 1, {-0.5, 0.6} at n = 2, {-0.5, 1} at
n = 4; and GR for the boosts.
Cosmologies: Planck 2018; Omega_cdm, h, A_s, n_s, w0, wa one at a time at -5/+5 sigma (wa: -5/+3 sigma); four joint
extremes of all six (hi, lo and two mixed) and the w0-wa combinations allowed by w0 + wa < 0.
Setup (production): 1024^3 particles on a 2048^3 force grid in 1024 Mpc/h, 30 steps, outputs z = 1, 0, R = 1 Mpc/h,
cut 1.0, bias spectra on 512^3, one phase (fixed amplitude, seed 1001), single precision. Finished runs are skipped.
Order: the fiducial cosmology for all models, then one parameter at a time, then the joint extremes. Runs go to results/mu_cosmology/<cosmology>/<model>/. Progress in progress.txt.
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from runs import steps_per_interval, write_input_power  # noqa: E402
from scan_planck import PLANCK  # noqa: E402

NBODY = os.path.join(os.path.dirname(ROOT), "SimulationCode", "nbody_float")
DRY = "--dry-run" in sys.argv
Z = [1.0, 0.0]
COMMON = dict(simulation_boxsize=1024.0, particle_Npart_1D=1024, force_nmesh=2048, pofk_nmesh=512, hybridbias_nmesh=512,
              hybridbias_smoothing_scale=1.0, hybridbias_kcut_over_nyquist=1.0, ic_use_gravity_model_GR=True,
              ic_fix_amplitude=True, ic_reverse_phases=False, ic_random_seed=1001, output_redshifts=Z,
              timestep_nsteps=steps_per_interval(30, Z), particle_allocation_factor=1.5)
MODELS = {"GR": dict(gravity_model="GR")}
# Realistic models only (scripts/mu_prior.py: equivalent DESI mu0 within ~3 sigma). Replaced 2026-11-17: A = -0.5 and
# A = 1 at n = 0.5 (growth -36% / +125%, excluded); their finished runs stay in the folders but are not in the list.
# A0.2_n1 and A0.6_n2 have the same growth change (+7.2%), to separate the effect of n.
for A, n in [(0.3, 3.0), (-0.35, 1.0), (-0.5, 4.0), (0.2, 1.0), (1.0, 4.0), (0.6, 2.0), (-0.5, 2.0)]:
    MODELS[f"A{A:g}_n{n:g}"] = dict(gravity_model="mu", gravity_model_mu_parametrization="powerlaw",
                                     gravity_model_mu_powerlaw_params=[A, n])

# Planck 2018 sigmas (as in scan_planck.py; n_s: 0.0042) and the +-5 sigma values (wa: +3 sigma, see scan_planck.py)
P = PLANCK
SIG = {"cosmology_OmegaCDM": 0.0073, "cosmology_h": 0.0054, "cosmology_ns": 0.0042}
lo = {k: round(P[k] - 5 * s, 6) for k, s in SIG.items()}
hi = {k: round(P[k] + 5 * s, 6) for k, s in SIG.items()}
lo["cosmology_As"], hi["cosmology_As"] = 1.957e-09, 2.25109e-09
W0 = {"lo": -1.4, "hi": -0.6}
WA = {"lo": -1.45, "hi": 0.87}
COSMOLOGIES = {"fiducial": {}}
for key in ["cosmology_OmegaCDM", "cosmology_h", "cosmology_As", "cosmology_ns"]:
    short = key.removeprefix("cosmology_")
    COSMOLOGIES[f"{short}_lo"] = {key: lo[key]}
    COSMOLOGIES[f"{short}_hi"] = {key: hi[key]}
for tag in ("lo", "hi"):
    COSMOLOGIES[f"w0_{tag}"] = {"cosmology_w0": W0[tag]}
    COSMOLOGIES[f"wa_{tag}"] = {"cosmology_wa": WA[tag]}
# Joint extremes: all up (more clustering / later dark energy), all down, and two mixed; w0 + wa < 0 in all
COSMOLOGIES["all_hi"] = {**hi, "cosmology_w0": W0["hi"], "cosmology_wa": WA["lo"]}
COSMOLOGIES["all_lo"] = {**lo, "cosmology_w0": W0["lo"], "cosmology_wa": WA["hi"]}
COSMOLOGIES["mix_a"] = {"cosmology_OmegaCDM": hi["cosmology_OmegaCDM"], "cosmology_h": lo["cosmology_h"],
                        "cosmology_As": lo["cosmology_As"], "cosmology_ns": hi["cosmology_ns"],
                        "cosmology_w0": W0["hi"], "cosmology_wa": WA["lo"]}
COSMOLOGIES["mix_b"] = {"cosmology_OmegaCDM": lo["cosmology_OmegaCDM"], "cosmology_h": hi["cosmology_h"],
                        "cosmology_As": hi["cosmology_As"], "cosmology_ns": lo["cosmology_ns"],
                        "cosmology_w0": W0["lo"], "cosmology_wa": WA["hi"]}
COSMOLOGIES["w0lo_walo"] = {"cosmology_w0": W0["lo"], "cosmology_wa": WA["lo"]}  # (hi, hi) would have w0 + wa > 0
COSMOLOGIES["w0hi_walo"] = {"cosmology_w0": W0["hi"], "cosmology_wa": WA["lo"]}
COSMOLOGIES["w0lo_wahi"] = {"cosmology_w0": W0["lo"], "cosmology_wa": WA["hi"]}


def arg(key, value):
    if isinstance(value, (list, tuple)):
        return f"{key}=[{', '.join(str(v) for v in value)}]"
    return f"{key}={value}"


def run(job):
    cosmo_label, model = job
    cosmo = {**P, **COSMOLOGIES[cosmo_label]}
    assert cosmo["cosmology_w0"] + cosmo["cosmology_wa"] < 0, cosmo_label
    path = os.path.join(HERE, "input", f"pofk_{cosmo_label}.txt")
    params = dict(cosmo, **COMMON, **MODELS[model], ic_type_of_input="powerspectrum", ic_input_filename=path)
    env = dict(os.environ, NBODY=NBODY)  # One run per node: runs.py gives each of the 64 tasks 2 threads
    cmd = ["python3", "scripts/runs.py", f"mu_cosmology/{cosmo_label}/{model}"] + \
          [arg(k, v) for k, v in params.items()] + ["ntasks=64"] + (["--dry-run"] if DRY else [])
    start = time.time()
    with open(os.path.join(HERE, "runs_stdout.txt"), "a") as out:
        subprocess.run(cmd, cwd=ROOT, env=env, stdout=out, stderr=subprocess.STDOUT)
    with open(os.path.join(HERE, "progress.txt"), "a") as f:
        f.write(f"{time.strftime('%T')} {cosmo_label}/{model} {time.time() - start:.0f} s\n")


if __name__ == "__main__":
    for label in COSMOLOGIES:  # The input P(k) once per cosmology, before the runs start in parallel
        cosmo = {**P, **COSMOLOGIES[label]}
        write_input_power({**COMMON, **cosmo}, os.path.join(HERE, "input", f"pofk_{label}.txt"))
    jobs = [(c, m) for c in COSMOLOGIES for m in MODELS]  # COSMOLOGIES is in priority order
    part, nparts = (int(a) for a in [a for a in sys.argv[1:] if not a.startswith("--")][:2])
    mine = jobs[part::nparts]
    print(len(jobs), "runs in all,", len(mine), "in part", part, "of", nparts)
    for job in mine:
        run(job)
    with open(os.path.join(HERE, "progress.txt"), "a") as f:
        f.write(f"ALLDONE part {part}\n")
