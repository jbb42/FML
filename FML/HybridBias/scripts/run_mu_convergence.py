#!/usr/bin/env python3
"""GR vs mu(a) = 1 + A a^n (powerlaw) convergence & cosmology test in single precision (nbody_float), stripped
down from run_weekend.py + scan_planck.py to finish within ~18 hours at high resolution (1024^3 particles on
2048^3 force mesh, ~320 GB peak RAM on euclid22):

    Models:      GR and six (A, n) pairs: A in {0.1, 0.3, 1.0} x n in {3, 1} (A1.0_n1 excluded from plots)
    fiducial:    Planck 2018 fiducial at L = 1024 Mpc/h, 30 steps, normal + reversed phases, all 7 models (14 runs)
    steps:       40 steps at L = 1024 Mpc/h, normal phase, for GR, A0.3_n3, A1.0_n3, A1.0_n1 (4 runs)
    box:         L = 512 Mpc/h, 30 steps, normal phase for GR and A0.3_n3 (2 runs)
    cosmologies: Planck 2018 +-5 sigma (wa: -5, +3 sigma) in Omega_m, H_0, A_s, w0, wa for GR and A0.3_n3,
                 normal phase, z = 0, with Zel'dovich prediction + field for CV (20 new runs; fiducial shared)

Writes figures/mu_boost/all_boost_mu.pdf (all 15 spectra per page, organized by parameter test; no single-spectrum files).

    python3 scripts/run_mu_convergence.py [fiducial] [steps] [box] [cosmologies] [--dry-run] [--plot-only]
"""
import glob
import os
import sys
import warnings

import numpy as np

import runs as runs_mod
from plot_page import save_pages
from runs import PLOT_ONLY, Run, run, steps_per_interval, summary, write_input_power
from scan_planck import PLANCK
from spectra import PAIRS, load_run

# Default to single-precision nbody_float (fits 1024^3 on 2048^3 in ~320 GB and is ~30% faster)
if "NBODY" not in os.environ:
    runs_mod.NBODY = os.path.join(os.path.dirname(runs_mod.ROOT), "SimulationCode", "nbody_float")

MODELS = [
    ("GR", None),
    ("A0.1_n3", [0.1, 3.0]),
    ("A0.3_n3", [0.3, 3.0]),
    ("A1.0_n3", [1.0, 3.0]),
    ("A0.1_n1", [0.1, 1.0]),
    ("A0.3_n1", [0.3, 1.0]),
    ("A1.0_n1", [1.0, 1.0]),
]
PLOT_MODELS = ["A0.1_n3", "A0.3_n3", "A1.0_n3", "A0.1_n1", "A0.3_n1"]  # A1.0_n1 excluded (too large a boost)
REDSHIFTS = [1.0, 0.0]
SEED = 1001
COMMON = {
    "force_nmesh": 2048,
    "particle_Npart_1D": 1024,
    "hybridbias_nmesh": 512,
    "hybridbias_smoothing_scale": 1.0,
    "hybridbias_kcut_over_nyquist": 1.0,
    "hybridbias_zeldovich_field": True,
    "ic_use_gravity_model_GR": True,
    "ic_fix_amplitude": True,
    "ic_random_seed": SEED,
}

# Cosmologies from scan_planck.py (extreme +-5 sigma points; wa: -5, +3 sigma)
h = PLANCK["cosmology_h"]
Omega_b = PLANCK["cosmology_Omegab"]
Omega_m = PLANCK["cosmology_OmegaCDM"] + Omega_b
ln_As = 3.044
COSMO_PARAMS = [
    ("OmegaCDM", r"\Omega_{\rm cdm}", PLANCK["cosmology_OmegaCDM"],
     [round(Omega_m + s * 0.0073 - Omega_b, 6) for s in (-5, 5)], ("-5σ", "+5σ")),
    ("h", r"h", h, [round(h + s * 0.0054, 6) for s in (-5, 5)], ("-5σ", "+5σ")),
    ("As", r"A_s", PLANCK["cosmology_As"],
     [float(f"{np.exp(ln_As + s * 0.014) * 1e-10:.6g}") for s in (-5, 5)], ("-5σ", "+5σ")),
    ("w0", r"w_0", -1.0, [round(-1.0 + s * 0.080, 6) for s in (-5, 5)], ("-5σ", "+5σ")),
    ("wa", r"w_a", 0.0, [round(s * 0.29, 6) for s in (-5, 3)], ("-5σ", "+3σ")),
]
EXTREME_COSMOLOGIES = []
for param, _, _, values, _ in COSMO_PARAMS:
    for v in values:
        EXTREME_COSMOLOGIES.append((f"{param}{v}", {**PLANCK, f"cosmology_{param}": v}))

runs = []


def model_parameters(label, params):
    if label == "GR":
        return {"gravity_model": "GR"}
    return {
        "gravity_model": "mu",
        "gravity_model_mu_parametrization": "powerlaw",
        "gravity_model_mu_powerlaw_params": params,
    }


def input_power_settings(cosmo_label, cosmo_params):
    """Compute CLASS input P(k) once per cosmology and return ic_* settings."""
    input_file = os.path.join(runs_mod.ROOT, "results", "mu_convergence", "input", f"pofk_{cosmo_label}.txt")
    if not PLOT_ONLY:
        write_input_power({**COMMON, "simulation_boxsize": 1024.0, **cosmo_params}, input_file)
    return {"ic_type_of_input": "powerspectrum", "ic_input_filename": input_file}


def run_case(cosmo_label, cosmo_params, ic_settings, models, box, nsteps, phases=("normal",), redshifts=REDSHIFTS,
             with_zeldovich=False):
    for label, params in models:
        for phase in phases:
            tag = f"{'gr' if label == 'GR' else 'mu_' + label}_{cosmo_label}_L{box}_f2048_n1024_s{nsteps}_{phase}"
            folder = run(
                f"mu_convergence/{tag}",
                **COMMON,
                **cosmo_params,
                **ic_settings,
                **model_parameters(label, params),
                simulation_boxsize=box,
                ic_reverse_phases=phase == "reversed",
                output_redshifts=redshifts,
                timestep_nsteps=steps_per_interval(nsteps, redshifts),
                hybridbias_zeldovich=with_zeldovich and phase == "normal",
            )
            runs.append(Run(folder, values=(cosmo_label, box, nsteps), model=label, seed=SEED,
                            reverse=phase == "reversed"))


fid_ic = input_power_settings("fiducial", PLANCK)
campaigns = [arg for arg in sys.argv[1:] if not arg.startswith("--")] or ["fiducial", "steps", "box", "cosmologies"]

if "fiducial" in campaigns:
    run_case("fiducial", PLANCK, fid_ic, MODELS, 1024, 30, phases=("normal", "reversed"), with_zeldovich=True)
if "steps" in campaigns:
    step_models = [MODELS[0], MODELS[2], MODELS[3], MODELS[6]]  # GR, A0.3_n3, A1.0_n3, A1.0_n1
    run_case("fiducial", PLANCK, fid_ic, step_models, 1024, 40)
if "box" in campaigns:
    run_case("fiducial", PLANCK, fid_ic, [MODELS[0], MODELS[2]], 512, 30)
if "cosmologies" in campaigns:
    cosmo_models = [MODELS[0], MODELS[2]]  # GR and A0.3_n3
    for cosmo_label, cosmo_params in EXTREME_COSMOLOGIES:
        ic = input_power_settings(cosmo_label, cosmo_params)
        run_case(cosmo_label, cosmo_params, ic, cosmo_models, 1024, 30, redshifts=[0.0], with_zeldovich=True)

summary()


# =====================================================================
# Plotting helpers: raw and control-variate (CV) spectra & boosts
# =====================================================================
def folder_of(model, cosmo="fiducial", box=1024, nsteps=30, phase="normal"):
    prefix = "gr" if model == "GR" else f"mu_{model}"
    return os.path.join(runs_mod.ROOT, "results", "mu_convergence",
                        f"{prefix}_{cosmo}_L{box}_f2048_n1024_s{nsteps}_{phase}")


def load_spectra(model, cosmo="fiducial", box=1024, nsteps=30, phase="normal", z=0.0, cv=True):
    """Return (k, P[15, len(k)]). If cv=True, uses hybridbias_cv_z*.txt if present, or constructs
    P - beta(k) (P_zfield - P_zel) using the matching (model, cosmo, L=1024, s=30, normal) prediction."""
    folder = folder_of(model, cosmo, box, nsteps, phase)
    snaps = glob.glob(os.path.join(folder, f"snapshot_*_z{z:.3f}"))
    if not snaps:
        return None, None
    snap = snaps[0]
    if cv:
        cv_file = os.path.join(snap, f"hybridbias_cv_z{z:.3f}.txt")
        if os.path.exists(cv_file):
            d = np.loadtxt(cv_file)
            return d[:, 0], d[:, 2:2 + len(PAIRS)].T
        # Construct CV from this run's zeldovichfield and the L=1024, s=30 normal run's zeldovich prediction
        zf_file = os.path.join(snap, f"hybridbias_zeldovichfield_z{z:.3f}.txt")
        ref_snap = glob.glob(os.path.join(folder_of(model, cosmo, 1024, 30, "normal"), f"snapshot_*_z{z:.3f}"))
        if os.path.exists(zf_file) and ref_snap:
            zp_file = os.path.join(ref_snap[0], f"hybridbias_zeldovich_z{z:.3f}.txt")
            if os.path.exists(zp_file):
                raw = np.loadtxt(os.path.join(snap, f"hybridbias_z{z:.3f}.txt"))
                zf = np.loadtxt(zf_file)
                zp = np.loadtxt(zp_file)
                k = raw[:, 0]
                zp_cols = zp[:, 2:17] if len(k) == len(zp) and np.allclose(k, zp[:, 0]) else \
                    np.column_stack([np.interp(k, zp[:, 0], zp[:, c]) for c in range(2, 17)])
                beta = 0.5 * (1.0 - np.tanh((k - 0.618) / 0.167))
                P_cv = raw[:, 2:17] - beta[:, None] * (zf[:, 2:17] - zp_cols)
                return k, P_cv.T
    return load_run(folder, z)


def safe_ratio(num, den):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.where(np.sign(num) == np.sign(den), num / den, np.nan)


def get_spectrum(model, cosmo="fiducial", box=1024, nsteps=30, phases=("normal",), z=0.0, cv=True):
    """Mean spectrum over phases."""
    arrs = []
    k_out = None
    for ph in phases:
        k, P = load_spectra(model, cosmo, box, nsteps, ph, z, cv=cv)
        if k is not None:
            k_out = k
            arrs.append(P)
    if not arrs:
        return None, None
    return k_out, np.mean(arrs, axis=0)


def get_boost(model, cosmo="fiducial", box=1024, nsteps=30, phases=("normal",), z=0.0, cv=True):
    """Ratio of phase-averaged P_mu / P_GR."""
    k, P_mu = get_spectrum(model, cosmo, box, nsteps, phases, z, cv=cv)
    _, P_gr = get_spectrum("GR", cosmo, box, nsteps, phases, z, cv=cv)
    if k is None or P_gr is None:
        return None, None
    return k, safe_ratio(P_mu, P_gr)


def curve(k, P, label, **style):
    return {"k": k, "P": P, "label": label, **style}


KRANGE = {"kmin": 0.01, "kmax": 1.0}
figures = []

# =====================================================================
# 1. Gravity model parameters (A, n): how A and n change the boost
# =====================================================================
MODEL_COLORS = {"A0.1": "C0", "A0.3": "C1", "A1.0": "C3"}
for z in [0.0, 1.0]:
    c_all = []
    for m in PLOT_MODELS:
        A_tag, n_tag = m.split("_")
        k, B = get_boost(m, phases=("normal", "reversed"), z=z, cv=True)
        if k is not None:
            c_all.append(curve(k, B, f"{A_tag}, {n_tag} (paired+CV)",
                               color=MODEL_COLORS[A_tag], ls="-" if n_tag == "n3" else "--"))
    figures.append({
        "title": f"1. Gravity models mu(a) = 1 + A a^n: boost mu / GR at z = {z} (solid: n=3, dashed: n=1)",
        "ratio": True, "ylabel": r"$P_{\mu} / P_{\rm GR}$", **KRANGE, "curves": c_all,
    })

    c_n3 = []
    for m in ["A0.1_n3", "A0.3_n3", "A1.0_n3"]:
        k, B = get_boost(m, phases=("normal", "reversed"), z=z, cv=True)
        if k is not None:
            c_n3.append(curve(k, B, m.replace("_", ", "), color=MODEL_COLORS[m.split("_")[0]]))
    figures.append({
        "title": f"1b. Varying amplitude A at fixed n = 3: boost mu / GR at z = {z} (paired+CV)",
        "ratio": True, "ylabel": r"$P_{\mu} / P_{\rm GR}$", **KRANGE, "curves": c_n3,
    })

    c_n = []
    for A_tag in ["A0.1", "A0.3"]:
        for n_tag, ls in [("n3", "-"), ("n1", "--")]:
            k, B = get_boost(f"{A_tag}_{n_tag}", phases=("normal", "reversed"), z=z, cv=True)
            if k is not None:
                c_n.append(curve(k, B, f"{A_tag}, {n_tag}", color=MODEL_COLORS[A_tag], ls=ls))
    figures.append({
        "title": f"1c. Varying exponent n (n=3 solid vs n=1 dashed) at A = 0.1, 0.3: boost mu / GR at z = {z}",
        "ratio": True, "ylabel": r"$P_{\mu} / P_{\rm GR}$", **KRANGE, "curves": c_n,
    })

# =====================================================================
# 2. Variance reduction: Normal vs Reversed vs Paired vs Control Variate
# =====================================================================
for z in [0.0, 1.0]:
    k, B_ref = get_boost("A0.3_n3", phases=("normal", "reversed"), z=z, cv=True)
    _, B_norm_raw = get_boost("A0.3_n3", phases=("normal",), z=z, cv=False)
    _, B_rev_raw = get_boost("A0.3_n3", phases=("reversed",), z=z, cv=False)
    _, B_pair_raw = get_boost("A0.3_n3", phases=("normal", "reversed"), z=z, cv=False)
    _, B_norm_cv = get_boost("A0.3_n3", phases=("normal",), z=z, cv=True)

    figures.append({
        "title": f"2a. Phases & Control Variate (A0.3_n3, z = {z}): boost mu / GR",
        "ratio": True, "ylabel": r"$P_{\mu} / P_{\rm GR}$", **KRANGE,
        "curves": [
            curve(k, B_norm_raw, "normal (raw)", color="C0", ls=":"),
            curve(k, B_rev_raw, "reversed (raw)", color="C1", ls=":"),
            curve(k, B_pair_raw, "paired normal+reversed (raw)", color="C2", ls="--"),
            curve(k, B_norm_cv, "normal (with CV)", color="C0", ls="-"),
            curve(k, B_ref, "paired normal+reversed (with CV)", color="k", ls="-"),
        ],
    })
    figures.append({
        "title": f"2b. Phases & Control Variate (A0.3_n3, z = {z}): boost relative to paired+CV (grey band: +-0.5%)",
        "ratio": True, "ylabel": r"$B / B_{\rm paired+CV}$", "ylim": (0.985, 1.015), "band": 0.005, **KRANGE,
        "curves": [
            curve(k, safe_ratio(B_norm_raw, B_ref), "normal (raw) / paired+CV", color="C0", ls=":"),
            curve(k, safe_ratio(B_rev_raw, B_ref), "reversed (raw) / paired+CV", color="C1", ls=":"),
            curve(k, safe_ratio(B_pair_raw, B_ref), "paired (raw) / paired+CV", color="C2", ls="--"),
            curve(k, safe_ratio(B_norm_cv, B_ref), "normal (with CV) / paired+CV", color="C3", ls="-"),
        ],
    })

# =====================================================================
# 3. Time-step convergence: 30 steps vs 40 steps (normal phase, L=1024)
# =====================================================================
for z in [0.0, 1.0]:
    c_step_b, c_step_bratio, c_step_pratio = [], [], []
    for m, col in [("A0.3_n3", "C1"), ("A1.0_n3", "C3")]:
        k, B30 = get_boost(m, nsteps=30, phases=("normal",), z=z, cv=True)
        _, B40 = get_boost(m, nsteps=40, phases=("normal",), z=z, cv=True)
        c_step_b.append(curve(k, B40, f"{m}, 40 steps", color=col, ls="-"))
        c_step_b.append(curve(k, B30, f"{m}, 30 steps", color=col, ls="--"))
        c_step_bratio.append(curve(k, safe_ratio(B30, B40), f"{m}: B(30 steps) / B(40 steps)", color=col, ls="-"))

    for m, col in [("GR", "k"), ("A0.3_n3", "C1"), ("A1.0_n3", "C3")]:
        k, P30 = get_spectrum(m, nsteps=30, phases=("normal",), z=z, cv=False)
        _, P40 = get_spectrum(m, nsteps=40, phases=("normal",), z=z, cv=False)
        c_step_pratio.append(curve(k, safe_ratio(P30, P40), f"{m}: P(30 steps) / P(40 steps)", color=col, ls="-"))

    figures.append({
        "title": f"3a. Time steps (30 dashed vs 40 solid, z = {z}): boost mu / GR",
        "ratio": True, "ylabel": r"$P_{\mu} / P_{\rm GR}$", **KRANGE, "curves": c_step_b,
    })
    figures.append({
        "title": f"3b. Time-step convergence of the BOOST at z = {z}: B(30 steps) / B(40 steps) (grey band: +-0.5%)",
        "ratio": True, "ylabel": r"$B(30) / B(40)$", "ylim": (0.98, 1.02), "band": 0.005, **KRANGE,
        "curves": c_step_bratio,
    })
    figures.append({
        "title": f"3c. Time-step convergence of the SPECTRA at z = {z}: P(30 steps) / P(40 steps) (grey band: +-1%)",
        "ratio": True, "ylabel": r"$P(30) / P(40)$", "ylim": (0.96, 1.04), "band": 0.01, **KRANGE,
        "curves": c_step_pratio,
    })

# =====================================================================
# 4. Box size & spatial resolution: L = 512 vs L = 1024 (s=30, normal)
# =====================================================================
for z in [0.0, 1.0]:
    k1024, B1024_cv = get_boost("A0.3_n3", box=1024, nsteps=30, phases=("normal",), z=z, cv=True)
    _, B1024_raw = get_boost("A0.3_n3", box=1024, nsteps=30, phases=("normal",), z=z, cv=False)
    k512, B512_cv = get_boost("A0.3_n3", box=512, nsteps=30, phases=("normal",), z=z, cv=True)
    _, B512_raw = get_boost("A0.3_n3", box=512, nsteps=30, phases=("normal",), z=z, cv=False)
    B512_cv_interp = np.array([np.interp(k1024, k512, b, left=np.nan, right=np.nan) for b in B512_cv])
    B512_raw_interp = np.array([np.interp(k1024, k512, b, left=np.nan, right=np.nan) for b in B512_raw])

    figures.append({
        "title": f"4a. Box size / resolution (L = 512 vs L = 1024 Mpc/h, A0.3_n3, z = {z}): boost mu / GR",
        "ratio": True, "ylabel": r"$P_{\mu} / P_{\rm GR}$", **KRANGE,
        "curves": [
            curve(k1024, B1024_cv, "L = 1024 Mpc/h (normal, with CV)", color="k", ls="-"),
            curve(k512, B512_cv, "L = 512 Mpc/h (2x finer res, normal, with CV)", color="C3", ls="-"),
            curve(k512, B512_raw, "L = 512 Mpc/h (2x finer res, normal, raw)", color="C0", ls=":"),
        ],
    })
    figures.append({
        "title": f"4b. Box size / resolution convergence of BOOST at z = {z}: B(L=512) / B(L=1024) (grey band: +-1%)",
        "ratio": True, "ylabel": r"$B(L=512) / B(L=1024)$", "ylim": (0.95, 1.05), "band": 0.01, **KRANGE,
        "curves": [
            curve(k1024, safe_ratio(B512_cv_interp, B1024_cv), "B(L=512 CV) / B(L=1024 CV)", color="C3", ls="-"),
            curve(k1024, safe_ratio(B512_raw_interp, B1024_raw), "B(L=512 raw) / B(L=1024 raw)", color="C0", ls=":"),
        ],
    })

# =====================================================================
# 5. Cosmological parameters: OmegaCDM, h, As, w0, wa at z = 0.0 (A0.3_n3)
# =====================================================================
k_fid, B_fid_cv = get_boost("A0.3_n3", cosmo="fiducial", phases=("normal",), z=0.0, cv=True)
_, P_fid_gr = get_spectrum("GR", cosmo="fiducial", phases=("normal",), z=0.0, cv=True)

summary_plus, summary_minus = [], []
for idx, (param, tex, fid_val, values, sig_labels) in enumerate(COSMO_PARAMS):
    c_boost, c_bratio, c_pratio = [curve(k_fid, B_fid_cv, f"Planck fiducial ({param}={fid_val})", color="k", ls="-")], [], []
    for v, s_lbl, col, ls, target_list in zip(values, sig_labels, ("C0", "C3"), ("--", "-"), (summary_minus, summary_plus)):
        clabel = f"{param}{v}"
        k, B_cv = get_boost("A0.3_n3", cosmo=clabel, phases=("normal",), z=0.0, cv=True)
        _, B_raw = get_boost("A0.3_n3", cosmo=clabel, phases=("normal",), z=0.0, cv=False)
        _, B_fid_raw = get_boost("A0.3_n3", cosmo="fiducial", phases=("normal",), z=0.0, cv=False)
        _, P_gr = get_spectrum("GR", cosmo=clabel, phases=("normal",), z=0.0, cv=True)

        c_boost.append(curve(k, B_cv, f"{param} = {v} ({s_lbl}, CV)", color=col, ls="-"))
        c_bratio.append(curve(k, safe_ratio(B_cv, B_fid_cv), f"{param} = {v} ({s_lbl}, CV)", color=col, ls="-"))
        c_bratio.append(curve(k, safe_ratio(B_raw, B_fid_raw), f"{param} = {v} ({s_lbl}, raw)", color=col, ls=":"))
        c_pratio.append(curve(k, safe_ratio(P_gr, P_fid_gr), f"{param} = {v} ({s_lbl})", color=col, ls="-"))
        target_list.append(curve(k, safe_ratio(B_cv, B_fid_cv), f"{param} {s_lbl} ({v})", color=f"C{idx}", ls="-"))

    figures.append({
        "title": f"5.{idx+1}a Varying {param} (A0.3_n3, z = 0): boost mu / GR (with CV)",
        "ratio": True, "ylabel": r"$P_{\mu} / P_{\rm GR}$", **KRANGE, "curves": c_boost,
    })
    figures.append({
        "title": f"5.{idx+1}b Varying {param} (A0.3_n3, z = 0): boost relative to Planck fiducial (grey band: +-0.5%)",
        "ratio": True, "ylabel": r"$B(\theta) / B_{\rm fid}$", "ylim": (0.97, 1.03), "band": 0.005, **KRANGE,
        "curves": c_bratio,
    })
    figures.append({
        "title": f"5.{idx+1}c Varying {param} (GR, z = 0): spectra relative to Planck fiducial (grey band: +-5%)",
        "ratio": True, "ylabel": r"$P_{\rm GR}(\theta) / P_{\rm GR, fid}$", "ylim": (0.5, 1.5), "band": 0.05, **KRANGE,
        "curves": c_pratio,
    })

figures.append({
    "title": "5.6a Summary (+5σ / +3σ): cosmological parameters' effect on boost B(theta) / B_fid at z = 0 (band: +-0.5%)",
    "ratio": True, "ylabel": r"$B(\theta) / B_{\rm fid}$", "ylim": (0.97, 1.03), "band": 0.005, **KRANGE,
    "curves": summary_plus,
})
figures.append({
    "title": "5.6b Summary (-5σ): cosmological parameters' effect on boost B(theta) / B_fid at z = 0 (band: +-0.5%)",
    "ratio": True, "ylabel": r"$B(\theta) / B_{\rm fid}$", "ylim": (0.97, 1.03), "band": 0.005, **KRANGE,
    "curves": summary_minus,
})

save_pages(figures, "mu_boost/all_boost_mu.pdf")
