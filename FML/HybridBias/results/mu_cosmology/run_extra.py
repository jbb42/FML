"""Extra jobs for the cosmology test, run by every node after its run_all.py / run_helper.py has finished (added
2026-11-17, nodes free until 09:00 on 2026-11-18):
  1. Omega_b at -5 / +5 sigma (Planck: omega_b = 0.02237 +- 0.00015, so Omega_b +- 0.00033 per sigma; Omega_cdm fixed),
     all models: Ob_lo, Ob_hi
  2. two more seeds of the fiducial cosmology, all models: fiducial_s1002, fiducial_s1003 (the realisation noise of
     the boost, to tell a real cosmology dependence from noise)
Same settings and models as run_all.py. A node claims a job by creating <cosmology>/<model>/.claim (O_EXCL), so no job
runs twice; finished jobs are skipped. Usage on each node: python3 run_extra.py
"""
import os
import subprocess
import sys
import time

import run_all
from runs import write_input_power  # noqa: E402 (run_all puts scripts/ on the path)

HERE = run_all.HERE
OB_SIGMA = 0.00015 / run_all.P["cosmology_h"] ** 2
EXTRA = {
    "Ob_lo": ({"cosmology_Omegab": round(run_all.P["cosmology_Omegab"] - 5 * OB_SIGMA, 6)}, {}),
    "Ob_hi": ({"cosmology_Omegab": round(run_all.P["cosmology_Omegab"] + 5 * OB_SIGMA, 6)}, {}),
    "fiducial_s1002": ({}, {"ic_random_seed": 1002}),
    "fiducial_s1003": ({}, {"ic_random_seed": 1003}),
}


def input_file(label):
    return os.path.join(HERE, "input", f"pofk_{'fiducial' if label.startswith('fiducial') else label}.txt")


def finished(folder):
    import glob
    return bool(glob.glob(os.path.join(folder, "snapshot_*_z0.000", "hybridbias_z0.000.txt")))


def claim(folder):
    os.makedirs(folder, exist_ok=True)
    try:
        os.close(os.open(os.path.join(folder, ".claim"), os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        return True
    except FileExistsError:
        return False


def run(label, model):
    cosmo_change, sim_change = EXTRA[label]
    cosmo = {**run_all.P, **cosmo_change}
    params = {**cosmo, **run_all.COMMON, **sim_change, **run_all.MODELS[model], "ic_type_of_input": "powerspectrum",
              "ic_input_filename": input_file(label)}
    cmd = ["python3", "scripts/runs.py", f"mu_cosmology/{label}/{model}"] + \
          [run_all.arg(k, v) for k, v in params.items()] + ["ntasks=64"]
    start = time.time()
    with open(os.path.join(HERE, "runs_stdout.txt"), "a") as out:
        subprocess.run(cmd, cwd=run_all.ROOT, env=dict(os.environ, NBODY=run_all.NBODY), stdout=out,
                       stderr=subprocess.STDOUT)
    with open(os.path.join(HERE, "progress.txt"), "a") as f:
        f.write(f"{time.strftime('%T')} {label}/{model} {time.time() - start:.0f} s (extra)\n")


if __name__ == "__main__":
    if "--inputs" in sys.argv:  # Once, before the nodes start
        for label, (change, _) in EXTRA.items():
            if not label.startswith("fiducial"):
                write_input_power({**run_all.COMMON, **run_all.P, **change}, input_file(label))
        sys.exit()
    for label in EXTRA:  # Priority order
        for model in run_all.MODELS:
            folder = os.path.join(HERE, label, model)
            if finished(folder) or not claim(folder):
                continue
            print(time.strftime("%T"), "extra:", label, model, flush=True)
            run(label, model)
