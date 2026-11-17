"""One-loop (and Zel'dovich) predictions for the fiducial runs of run_mu_convergence.py, without time stepping
(hybridbias_zeldovich_field_only), with the same settings: Planck 2018, 1024^3 lattice in 1024 Mpc/h, bias spectra on
512^3, R = 1 Mpc/h, cut 1.0, z = 0 (the Zel'dovich prediction is in the runs). Run on euclid26 (native single-precision build):
    python3 run_theory.py KERNELS NPARALLEL MODEL [MODEL ...]     (KERNELS: EdS or model; MODEL: GR, A0.3_n3, ...)
Folders: results/mu_convergence/theory/<model>_<kernels>/
Needs the one-loop prediction in the field-only mode: the last argument of write_hybrid_bias_predictions in the
hybridbias_zeldovich_field_only branch of SimulationCode/src/SimulationStep.h set to true (used on 2026-11-17 for
figures/mu_boost/theory_vs_sim.pdf, then reverted since production uses the Zel'dovich prediction only).
"""
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from scan_planck import PLANCK  # noqa: E402

NBODY = os.path.join(os.path.dirname(ROOT), "build_euclid22_float", "nbody")
COMMON = dict(simulation_boxsize=1024.0, particle_Npart_1D=1024, force_nmesh=2048, hybridbias_nmesh=512,
              hybridbias_smoothing_scale=1.0, hybridbias_kcut_over_nyquist=1.0, ic_use_gravity_model_GR=True,
              ic_fix_amplitude=True, ic_random_seed=1001, output_redshifts=[0.0],
              hybridbias_zeldovich_field=True, hybridbias_zeldovich_field_only=True, hybridbias_zeldovich=False,
              hybridbias_oneloop=True, ic_type_of_input="powerspectrum",
              ic_input_filename=os.path.join(ROOT, "results", "mu_convergence", "input", "pofk_fiducial.txt"))


def model_params(m):
    if m == "GR":
        return dict(gravity_model="GR")
    A, n = (float(x[1:]) for x in m.split("_"))
    return dict(gravity_model="mu", gravity_model_mu_parametrization="powerlaw", gravity_model_mu_powerlaw_params=[A, n])


def arg(key, value):
    return f"{key}=[{', '.join(str(v) for v in value)}]" if isinstance(value, list) else f"{key}={value}"


def run(model, kernels, ntasks):
    params = dict(PLANCK, **COMMON, **model_params(model), hybridbias_oneloop_kernels=kernels)
    env = dict(os.environ, NBODY=NBODY, OMP_NUM_THREADS=str(max(1, 128 // (ntasks * NPAR))), I_MPI_PIN="off")
    start = time.time()
    with open(os.path.join(HERE, "runs_stdout.txt"), "a") as out:
        subprocess.run(["python3", "scripts/runs.py", f"mu_convergence/theory/{model}_{kernels}"]
                       + [arg(k, v) for k, v in params.items()] + [f"ntasks={ntasks}"],
                       cwd=ROOT, env=env, stdout=out, stderr=subprocess.STDOUT)
    with open(os.path.join(HERE, "progress.txt"), "a") as f:
        f.write(f"{time.strftime('%T')} {model}_{kernels} {time.time() - start:.0f} s\n")


if __name__ == "__main__":
    kernels, NPAR, models = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
    with ThreadPoolExecutor(max_workers=NPAR) as pool:
        list(pool.map(lambda m: run(m, kernels, 64 // NPAR if NPAR > 1 else 64), models))
