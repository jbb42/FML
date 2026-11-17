"""Pre-production tests on euclid25 for the mu(a) = 1 + A a^n campaign. Run there: python3 run_all.py [A] [C] [B] [D]

A: is the control variate unbiased for mu models and w0wa? Zel'dovich fields without time stepping (32 seeds, both
   phases, 1/64-scale box) against the prediction, for A0.3_n3 (Planck) and A1.0_n1 on w0 = -0.6, wa = -1.5
B: corners of the hypercube incl. negative A, 512^3 on 1024^3 in 1024 Mpc/h, 30 steps, default particle buffer
C: cost of the Zel'dovich prediction for hybridbias_nmesh 256, 384, 512 at full size (no time stepping)
D: n_s +-5 sigma for GR and A0.3_n3, the settings of run_mu_convergence.py (results go to results/mu_convergence)
One line per finished run in progress.txt; full-size runs under watchdog.sh.
"""
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from runs import steps_per_interval, write_input_power  # noqa: E402
from scan_planck import PLANCK  # noqa: E402

NBODY = os.path.join(os.path.dirname(ROOT), "SimulationCode", "nbody_float")
REL = os.path.relpath(HERE, os.path.join(ROOT, "results"))
OP = dict(hybridbias_smoothing_scale=1.0, hybridbias_kcut_over_nyquist=1.0, ic_use_gravity_model_GR=True,
          ic_fix_amplitude=True)


def mu(A, n):
    return dict(gravity_model="mu", gravity_model_mu_parametrization="powerlaw", gravity_model_mu_powerlaw_params=[A, n])


def arg(key, value):
    if isinstance(value, (list, tuple)):
        return f"{key}=[{', '.join(str(v) for v in value)}]"
    return f"{key}={value}"


def run(name, params, ntasks, parallel=True, watch=False, rel=REL):
    """parallel (always, as two streams of jobs share the node): one thread per task and no pinning. The memory is
    watched for the whole node by main()."""
    env = dict(os.environ, NBODY=NBODY)
    if parallel:
        env.update(OMP_NUM_THREADS="1", I_MPI_PIN="off")
    start = time.time()
    with open(os.path.join(HERE, "runs_stdout.txt"), "a") as out:
        subprocess.run(["python3", "scripts/runs.py", f"{rel}/{name}" + ("/r" if rel == REL else "")]
                       + [arg(k, v) for k, v in params.items()] + [f"ntasks={ntasks}"],
                       cwd=ROOT, env=env, stdout=out, stderr=subprocess.STDOUT)
    with open(os.path.join(HERE, "progress.txt"), "a") as f:
        f.write(f"{time.strftime('%T')} {name} {time.time() - start:.0f} s\n")


def input_power(label, cosmo, folder):
    path = os.path.join(folder, f"pofk_{label}.txt")
    write_input_power({"simulation_boxsize": 1024.0, "particle_Npart_1D": 1024, **cosmo}, path)
    return dict(ic_type_of_input="powerspectrum", ic_input_filename=path)


def test_A():
    inp = os.path.join(HERE, "input")
    corner = {**PLANCK, "cosmology_model": "w0waCDM", "cosmology_w0": -0.6, "cosmology_wa": -1.5}
    cases = {"A0.3_n3_planck": ({**PLANCK, **input_power("planck", PLANCK, inp)}, mu(0.3, 3.0)),
             "A1.0_n1_w0wa": ({**corner, **input_power("w0wa_corner", corner, inp)}, mu(1.0, 1.0))}
    jobs = []
    for tag, (cosmo, model) in cases.items():
        base = dict(cosmo, **model, **OP, simulation_boxsize=256.0, particle_Npart_1D=256, force_nmesh=256,
                    pofk_nmesh=256, hybridbias_nmesh=256, hybridbias_zeldovich_field=True,
                    hybridbias_zeldovich_field_only=True)
        jobs.append((f"A/{tag}_pred", dict(base, hybridbias_zeldovich=True, ic_random_seed=1)))
        for seed in range(1, 33):
            for rev in (False, True):
                jobs.append((f"A/{tag}_s{seed}_{'rev' if rev else 'nor'}",
                             dict(base, ic_random_seed=seed, ic_reverse_phases=rev)))
    with ThreadPoolExecutor(max_workers=4) as pool:  # 4 x 16 tasks: half the node (D runs alongside)
        list(pool.map(lambda j: run(j[0], j[1], 16), jobs))


def test_C():
    fid = {**PLANCK, **input_power("planck", PLANCK, os.path.join(HERE, "input"))}
    for nm in (256, 384, 512):
        run(f"C/nmesh{nm}", dict(fid, **OP, pofk_nmesh=nm, hybridbias_nmesh=nm, hybridbias_zeldovich_field=True,
                                 hybridbias_zeldovich_field_only=True, hybridbias_zeldovich=True,
                                 output_redshifts=[0.0]), 64)


def test_B():
    sig = dict(OmegaCDM=0.0073, h=0.0054, ns=0.0042)
    lo = {"cosmology_OmegaCDM": round(PLANCK["cosmology_OmegaCDM"] - 5 * sig["OmegaCDM"], 6),
          "cosmology_h": round(PLANCK["cosmology_h"] - 5 * sig["h"], 6), "cosmology_As": 1.957e-09,
          "cosmology_ns": round(PLANCK["cosmology_ns"] - 5 * sig["ns"], 6)}
    hi = {"cosmology_OmegaCDM": round(PLANCK["cosmology_OmegaCDM"] + 5 * sig["OmegaCDM"], 6),
          "cosmology_h": round(PLANCK["cosmology_h"] + 5 * sig["h"], 6), "cosmology_As": 2.25109e-09,
          "cosmology_ns": round(PLANCK["cosmology_ns"] + 5 * sig["ns"], 6)}
    w = lambda w0, wa: {"cosmology_model": "w0waCDM", "cosmology_w0": w0, "cosmology_wa": wa}
    corners = {
        "Am0.5_n0.5_w0-1.4_wa0.5_low": ({**PLANCK, **lo, **w(-1.4, 0.5)}, mu(-0.5, 0.5)),
        "Am0.5_n4_w0-0.5_wa-2": ({**PLANCK, **w(-0.5, -2.0)}, mu(-0.5, 4.0)),
        "A1_n0.5_w0-0.5_wa-2_high": ({**PLANCK, **hi, **w(-0.5, -2.0)}, mu(1.0, 0.5)),
        "A1_n0.5_w0-1.4_wa0.5": ({**PLANCK, **w(-1.4, 0.5)}, mu(1.0, 0.5)),
        "A1_n4_w0-0.5_wa-2_high": ({**PLANCK, **hi, **w(-0.5, -2.0)}, mu(1.0, 4.0)),
    }
    for tag, (cosmo, model) in corners.items():
        run(f"B/{tag}", dict(cosmo, **model, **OP, particle_Npart_1D=512, force_nmesh=1024, pofk_nmesh=512,
                             hybridbias_nmesh=512, output_redshifts=[0.0], timestep_nsteps=[30],
                             ic_random_seed=1001), 64)


def test_D():
    """n_s +-5 sigma with the settings and folder names of run_mu_convergence.py, into results/mu_convergence."""
    common = dict(force_nmesh=2048, particle_Npart_1D=1024, hybridbias_nmesh=512, hybridbias_smoothing_scale=1.0,
                  hybridbias_kcut_over_nyquist=1.0, hybridbias_zeldovich_field=True, ic_use_gravity_model_GR=True,
                  ic_fix_amplitude=True, ic_random_seed=1001, simulation_boxsize=1024, ic_reverse_phases=False,
                  output_redshifts=[0.0], timestep_nsteps=steps_per_interval(30, [0.0]), hybridbias_zeldovich=True)
    inp = os.path.join(ROOT, "results", "mu_convergence", "input")
    for v in (round(PLANCK["cosmology_ns"] - 5 * 0.0042, 6), round(PLANCK["cosmology_ns"] + 5 * 0.0042, 6)):
        cosmo = {**PLANCK, "cosmology_ns": v}
        ic = input_power(f"ns{v}", cosmo, inp)
        for label, model in (("gr", dict(gravity_model="GR")), ("mu_A0.3_n3", mu(0.3, 3.0))):
            run(f"{label}_ns{v}_L1024_f2048_n1024_s30_normal", dict(cosmo, **ic, **model, **common), 64,
                rel="mu_convergence")


def stream(tests):
    for t in tests:
        globals()["test_" + t]()


if __name__ == "__main__":
    # Two streams side by side, 64 cores each: D (320 GB per run, one at a time) then C, and A then B (small)
    watch = subprocess.Popen(["bash", os.path.join(HERE, "watchdog.sh"), os.path.join(HERE, "watch_node.txt")])
    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            list(pool.map(stream, [["D", "C"], ["A", "B"]]))
    finally:
        watch.terminate()
    with open(os.path.join(HERE, "progress.txt"), "a") as f:
        f.write("ALLDONE\n")
