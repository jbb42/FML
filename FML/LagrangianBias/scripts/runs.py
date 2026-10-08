"""Run COLASolver with parameterfile.lua, changing any of its parameters.

From the command line, one run in results/<folder>/ (lists as Python lists, true/false as True/False):
    python3 scripts/runs.py tests/mesh1024 pofk_nmesh=1024 cosmology_h=0.7 [ntasks=64] [--dry-run]

From a run script:
    run("tests/mesh1024", pofk_nmesh=1024)     # one run, returns its folder
    planck = scan("planck", fixed={...}, one_at_a_time={...}, gravity=["F5"], seeds=[1001])

Finished runs are skipped, so a script can be restarted or extended. Options for every script:
    --dry-run     only write the parameter files and print the commands, no figures
    --plot-only   run nothing, only make the figures from the runs that have finished
Changing cosmology_* parameters also computes the input P(k) with CLASS (write_input_power below).
"""
import ast
import itertools
import os
import re
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass

import numpy as np


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # The LagrangianBias folder
NBODY = os.path.join(os.path.dirname(ROOT), "COLASolver", "nbody")
BASE_PARAMETER_FILE = os.path.join(ROOT, "parameterfile.lua")
DRY_RUN, PLOT_ONLY = "--dry-run" in sys.argv, "--plot-only" in sys.argv
counts = {"run": 0, "already finished": 0, "failed": 0}


def base_parameter(name):
    """The value of a parameter in parameterfile.lua, as a string without quotes and trailing -- comment."""
    value = re.search(rf"^\s*{name}\s*=\s*(.*?)\s*(?:--.*)?$", open(BASE_PARAMETER_FILE).read(), re.M).group(1)
    return value.strip('"')


def fiducial(name):
    """The value of a parameter in parameterfile.lua: a number, or a string such as "GR"."""
    try:
        return float(base_parameter(name))
    except ValueError:
        return base_parameter(name)


def around_fiducial(name, fractions):
    """The fiducial value changed by these fractions: [-0.5, 0.1] gives fiducial x 0.5 and x 1.1."""
    return [float(f"{fiducial(name) * (1 + fraction):.6g}") for fraction in fractions]


def lua(value):
    """A Python value as Lua: True -> true, "GR" -> "GR" (with quotes), [10, 20] -> {10, 20}."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return f'"{value}"'
    if isinstance(value, (list, tuple)):
        return "{" + ", ".join(lua(v) for v in value) + "}"
    return str(value)


# The CLASS executable (https://github.com/lesgourg/class_public), or set the environment variable CLASS
CLASS = os.environ.get("CLASS", "/mn/stornext/u3/jonasbbe/pc/Dokumenter/simsilun_backreaction/class_public/class")

# The parameters of parameterfile.lua and the CLASS parameter each one sets
LUA_TO_CLASS = {
    "cosmology_h": "h",
    "cosmology_Omegab": "Omega_b",
    "cosmology_OmegaCDM": "Omega_cdm",
    "cosmology_OmegaK": "Omega_k",
    "cosmology_TCMB_kelvin": "T_cmb",
    "cosmology_As": "A_s",
    "cosmology_ns": "n_s",
    "cosmology_kpivot_mpc": "k_pivot",
    "ic_input_redshift": "z_pk",
}
# The same for cosmology_model = "w0waCDM", where a fluid with w(a) = w0 + wa (1 - a) replaces the cosmological constant
LUA_TO_CLASS_W0WA = {**LUA_TO_CLASS, "cosmology_w0": "w0_fld", "cosmology_wa": "wa_fld"}


def neutrinos(Neffective, OmegaMNu):
    """CLASS neutrino parameters for COLASolver's neutrinos: 3 species with temperature
    T_CMB (N_eff / 3)^(1/4) (4/11)^(1/3), sharing the mass that gives OmegaMNu."""
    if OmegaMNu == 0:
        return {"N_ur": Neffective}
    return {"N_ur": 0, "N_ncdm": 1, "deg_ncdm": 3, "Omega_ncdm": OmegaMNu,
            "T_ncdm": (Neffective / 3) ** 0.25 * (4 / 11) ** (1 / 3)}


def write_input_power(parameters, path):
    """Write the linear CDM+baryon P(k) at ic_input_redshift for the cosmology in parameters (the rest from
    parameterfile.lua) to path, as columns k (h/Mpc) and P(k) (Mpc/h)^3. The CLASS input is kept next to it, in
    <path>_class.ini, and nothing is done if that is unchanged."""
    def value(name):
        return parameters.get(name, base_parameter(name))

    # FML needs P(k) up to the corner of the IC grid (ic_nmesh = particle_Npart_1D), sqrt(3) times its Nyquist frequency
    kmax = 3**0.5 * np.pi * float(value("particle_Npart_1D")) / float(value("simulation_boxsize"))

    model = value("cosmology_model")
    mapping = {"LCDM": LUA_TO_CLASS, "w0waCDM": LUA_TO_CLASS_W0WA}.get(model)
    unknown = [p for p in parameters if p.startswith("cosmology_") and p not in (mapping or {})
               and p not in ("cosmology_model", "cosmology_Neffective", "cosmology_OmegaMNu")]
    if mapping is None or unknown:
        raise SystemExit(f"The input P(k) is set up for cosmology_model = LCDM or w0waCDM, not {model}, "
                         f"and for the parameters in LUA_TO_CLASS, not {unknown}")

    class_parameters = {class_name: value(lua_name) for lua_name, class_name in mapping.items()}
    if model == "w0waCDM":
        class_parameters["Omega_Lambda"] = 0  # The fluid then fills the rest of the energy budget
    class_parameters.update(neutrinos(float(value("cosmology_Neffective")), float(value("cosmology_OmegaMNu"))))
    ini = "output = mPk\n" + f"P_k_max_h/Mpc = {max(100.0, 1.1 * kmax)}\nk_per_decade_for_pk = 32\n"
    ini += "".join(f"{name} = {v}\n" for name, v in class_parameters.items())
    ini_file = path.removesuffix(".txt") + "_class.ini"
    if os.path.exists(path) and os.path.exists(ini_file) and open(ini_file).read() == ini:
        return  # Already computed with exactly this CLASS input
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:  # CLASS only accepts short output paths
        with open(os.path.join(tmp, "input.ini"), "w") as f:
            f.write(ini + f"root = {tmp}/\n")
        result = subprocess.run([CLASS, os.path.join(tmp, "input.ini")], capture_output=True, text=True)
        if result.returncode != 0:
            raise SystemExit(f"CLASS failed for {path}:\n{ini}{result.stdout[-2000:]}{result.stderr[-2000:]}")
        massive_neutrinos = float(value("cosmology_OmegaMNu")) > 0
        k, P = np.loadtxt(os.path.join(tmp, "00_pk_cb.dat" if massive_neutrinos else "00_pk.dat"), unpack=True)  # CDM + baryons
    np.savetxt(path, np.column_stack([k, P]), header="k (h/Mpc)   P(k) (Mpc/h)^3, linear CDM+baryon from CLASS")
    with open(ini_file, "w") as f:  # The CLASS input, as a record and to know when to recompute
        f.write(ini)
    print(f"Computed the input P(k) with CLASS: {path}")


def is_finished(folder):
    return os.path.exists(os.path.join(folder, "snapshot_TestSim_z0.000", "pofk_44.txt"))


def run(name, ntasks=64, **parameters):
    """Run parameterfile.lua with these parameters changed in results/<name>/, unless it has finished. Returns the
    folder."""
    folder = os.path.join(ROOT, "results", name)
    if PLOT_ONLY:
        return folder
    if is_finished(folder):
        counts["already finished"] += 1
        return folder
    os.makedirs(folder, exist_ok=True)
    if any(p.startswith("cosmology_") for p in parameters) and "ic_input_filename" not in parameters:
        input_file = os.path.join(folder, "input_power.txt")
        write_input_power(parameters, input_file)
        parameters.update(ic_type_of_input="powerspectrum", ic_input_filename=input_file)

    text = open(BASE_PARAMETER_FILE).read()
    for key, value in {**parameters, "output_folder": folder}.items():
        text = re.sub(rf"^(\s*){key}\s*=.*$", rf"\g<1>{key} = {lua(value)}", text, flags=re.M)
    parameter_file = os.path.join(folder, "parameterfile.lua")
    open(parameter_file, "w").write(text)

    command = ["mpirun", "-np", str(ntasks), NBODY, parameter_file]
    counts["run"] += 1
    if DRY_RUN:
        print("DRY_RUN:", " ".join(command))
        return folder
    print(f"{time.strftime('%F %T')} Running {name} with {ntasks} tasks")
    current_log = os.path.join(os.path.dirname(folder), "current_log.txt")  # Follow with tail -F
    if os.path.lexists(current_log):
        os.remove(current_log)
    os.symlink(os.path.join(folder, "log.txt"), current_log)
    # Each MPI task gets its share of the hardware threads, and Intel MPI pins it to its own cores
    env = {"OMP_NUM_THREADS": str(max(1, os.cpu_count() // ntasks)), "I_MPI_PIN_DOMAIN": "omp", **os.environ}
    start = time.time()
    with open(os.path.join(folder, "log.txt"), "w") as log:  # Input paths in the parameter file are relative to ROOT
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT, env=env)
    if is_finished(folder):
        print(f"{time.strftime('%F %T')} Finished {name} in {(time.time() - start) / 60:.0f} min")
    else:
        print(f"{time.strftime('%F %T')} ERROR: {name} failed, see {folder}/log.txt")
        counts["failed"] += 1
    return folder


def summary():
    """Print how many runs were done, skipped and failed; with --dry-run, stop the script here (no figures)."""
    if not PLOT_ONLY:
        print(", ".join(f"{n} {what}" for what, n in counts.items()))
    if DRY_RUN:
        sys.exit()


def gravity_parameters(model):
    """The parameters of a gravity model: "GR", or "F5" for f(R) with fofr0 = 1e-5 (F4: 1e-4, ...)."""
    if model == "GR":
        return {"gravity_model": "GR"}
    return {"gravity_model": "f(R)", "gravity_model_fofr_fofr0": 10.0 ** -int(model[1:])}


def steps_per_interval(nsteps, redshifts):
    """nsteps time steps split over the intervals between outputs at these redshifts, in proportion to the change in
    a, so the steps are as long as in a run with only the last output. Rounded down, then the intervals with the
    largest remainders get one more step each, until they add up to nsteps."""
    a = [1 / (1 + float(base_parameter("ic_initial_redshift")))] + [1 / (1 + z) for z in redshifts]
    exact = [nsteps * (a[i + 1] - a[i]) / (1 - a[0]) for i in range(len(redshifts))]
    steps = [int(x) for x in exact]
    for i in sorted(range(len(exact)), key=lambda i: steps[i] - exact[i])[:nsteps - sum(steps)]:
        steps[i] += 1
    return steps


@dataclass
class Run:
    folder: str
    values: tuple  # Of the parameters that are varied
    model: str     # Gravity model, "GR", "F5", ..., or None
    seed: int
    reverse: bool  # Reversed phases


@dataclass
class Scan:
    name: str
    parameters: list  # Names of the varied parameters, without "cosmology_"
    reference: tuple  # Their fiducial values, or None
    models: list
    runs: list


def short(name):
    return name.removeprefix("cosmology_")


def scan(name, vary=None, one_at_a_time=None, fixed={}, gravity=None, seeds=None, phases=("normal", "reversed"),
         ntasks=64):
    """Runs in results/scans/<name>/<values>[_<gravity>]_seed<seed>_<phase>/, where <values> are the parameters that
    differ from the fiducial (or "fiducial"). Give either
        vary:          {parameter: [values]}, to run every combination, or
        one_at_a_time: {parameter: [values]}, to change one parameter at a time, the others fiducial (from fixed, or
                       else parameterfile.lua), plus one fiducial run.
    fixed:   {parameter: value} for all runs
    gravity: f(R) models, e.g. ["F5"]: every point is then run in GR and these, from the same initial conditions
    seeds:   default the one in parameterfile.lua
    Returns the Scan, for scan_figures() in spectra.py."""
    parameters = list(vary or one_at_a_time)
    if vary:
        reference, points = None, list(itertools.product(*vary.values()))
    else:
        reference = tuple(fixed.get(p, fiducial(p)) for p in parameters)
        points = [reference] + [reference[:i] + (value,) + reference[i + 1:]
                                for i, values in enumerate(one_at_a_time.values()) for value in values if value != reference[i]]
    models = ["GR"] + gravity if gravity else [None]
    runs = []
    for values in points:
        changed = {p: v for p, v, ref in zip(parameters, values, reference or [None] * len(values)) if v != ref}
        settings = {**fixed, **dict(zip(parameters, values))}
        if any(p.startswith("cosmology_") for p in settings):  # The input P(k), computed once per cosmology
            cosmology = "_".join(f"{short(p)}{v}" for p, v in changed.items() if p.startswith("cosmology_"))
            input_file = os.path.join(ROOT, "results", "scans", name, "input", f"pofk_{cosmology or 'fiducial'}.txt")
            if not PLOT_ONLY:
                write_input_power(settings, input_file)
            settings.update(ic_type_of_input="powerspectrum", ic_input_filename=input_file)
        label = "_".join(f"{short(p)}{v}" for p, v in changed.items()) or "fiducial"
        for model in models:
            model_settings = {**gravity_parameters(model), "ic_use_gravity_model_GR": True} if model else {}
            for seed in seeds or [int(base_parameter("ic_random_seed"))]:
                for phase in phases:
                    folder = run(f"scans/{name}/{label}{'_' + model if model else ''}_seed{seed}_{phase}", ntasks,
                                 **settings, **model_settings, ic_random_seed=seed, ic_reverse_phases=phase == "reversed")
                    runs.append(Run(folder, values, model, seed, phase == "reversed"))
    summary()
    return Scan(name, [short(p) for p in parameters], reference, models, runs)


if __name__ == "__main__":
    def parse(value):
        """1024 -> 1024, [1.0, 0.0] -> list, True/true -> True, f(R) -> "f(R)"."""
        if value.lower() in ("true", "false"):
            return value.lower() == "true"
        try:
            return ast.literal_eval(value)
        except (ValueError, SyntaxError):
            return value

    name, *assignments = [arg for arg in sys.argv[1:] if not arg.startswith("--")]
    parameters = {key: parse(value) for key, value in (assignment.split("=", 1) for assignment in assignments)}
    run(name, **parameters)
    summary()
