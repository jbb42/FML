"""Run COLASolver simulations for the campaign scripts in this folder.

run("campaign/name", ...) runs parameterfile.lua with some parameters changed, in results/campaign/name/, which then
holds the parameter file, log.txt and the output. Finished runs are skipped, so a campaign can be restarted.
Run a campaign script with --dry-run to only write the parameter files and print the commands.
"""
import os
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # The LagrangianBias folder
NBODY = os.path.join(os.path.dirname(ROOT), "COLASolver", "nbody")
BASE_PARAMETER_FILE = os.path.join(ROOT, "parameterfile.lua")
DRY_RUN = "--dry-run" in sys.argv
counts = {"run": 0, "already finished": 0, "failed": 0}


def lua(value):
    """A Python value as Lua: True -> true, "GR" -> "GR" (with quotes), [10, 20] -> {10, 20}."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return f'"{value}"'
    if isinstance(value, (list, tuple)):
        return "{" + ", ".join(lua(v) for v in value) + "}"
    return str(value)


def base_parameter(name):
    """The value of a parameter in parameterfile.lua, as a string without quotes and trailing -- comment."""
    value = re.search(rf"^\s*{name}\s*=\s*(.*?)\s*(?:--.*)?$", open(BASE_PARAMETER_FILE).read(), re.M).group(1)
    return value.strip('"')


def write_parameter_file(path, parameters):
    """parameterfile.lua with the given parameters replaced."""
    text = open(BASE_PARAMETER_FILE).read()
    for name, value in parameters.items():
        text = re.sub(rf"^(\s*){name}\s*=.*$", rf"\g<1>{name} = {lua(value)}", text, flags=re.M)
    with open(path, "w") as f:
        f.write(text)


def run(name, ntasks=64, **parameters):
    """Run parameterfile.lua with these parameters in results/<name>, unless it has already finished."""
    out_dir = os.path.join(ROOT, "results", name)
    if os.path.exists(os.path.join(out_dir, "snapshot_TestSim_z0.000", "pofk_bias_info.txt")):
        counts["already finished"] += 1
        return
    os.makedirs(out_dir, exist_ok=True)
    parameter_file = os.path.join(out_dir, "parameterfile.lua")
    write_parameter_file(parameter_file, {**parameters, "output_folder": out_dir})
    command = ["mpirun", "-np", str(ntasks), NBODY, parameter_file]
    counts["run"] += 1
    if DRY_RUN:
        print("DRY_RUN:", " ".join(command))
        return

    print(f"{time.strftime('%F %T')} Running {name} with {ntasks} tasks")
    current_log = os.path.join(os.path.dirname(out_dir), "current_log.txt")  # Follow with tail -F
    if os.path.lexists(current_log):
        os.remove(current_log)
    os.symlink(os.path.join(out_dir, "log.txt"), current_log)
    start = time.time()
    with open(os.path.join(out_dir, "log.txt"), "w") as log:
        status = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT).returncode  # Input paths are relative to ROOT
    if status == 0 and os.path.exists(os.path.join(out_dir, "snapshot_TestSim_z0.000", "pofk_bias_info.txt")):
        print(f"{time.strftime('%F %T')} Finished {name} in {(time.time() - start) / 60:.0f} min")
    else:
        print(f"{time.strftime('%F %T')} ERROR: {name} failed, see {out_dir}/log.txt")
        counts["failed"] += 1


def summary():
    print(", ".join(f"{n} {what}" for what, n in counts.items()))


def steps_per_interval(nsteps, redshifts):
    """Split nsteps time steps over the intervals between the outputs at these redshifts, in proportion to the change
    in a, so the step size matches a run with a single output. Rounded so they add up to nsteps."""
    a_ini = 1.0 / (1.0 + float(base_parameter("ic_initial_redshift")))
    a = [a_ini] + [1.0 / (1.0 + z) for z in redshifts]
    exact = [nsteps * (a[i + 1] - a[i]) / (1.0 - a_ini) for i in range(len(redshifts))]
    steps = [int(x) for x in exact]
    largest_remainders = sorted(range(len(exact)), key=lambda i: steps[i] - exact[i])
    for i in largest_remainders[: nsteps - sum(steps)]:
        steps[i] += 1
    return steps
