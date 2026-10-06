"""The linear input P(k) for scans that change cosmology_* parameters, computed by running the CLASS executable."""
import glob
import os
import subprocess

import numpy as np

from simulations import base_parameter

# The CLASS executable (https://github.com/lesgourg/class_public), or set the environment variable CLASS
CLASS = os.environ.get("CLASS", "/mn/stornext/u3/jonasbbe/pc/Dokumenter/simsilun_backreaction/class_public/class")


def write_input_power(parameters, path):
    """Write the linear CDM+baryon P(k) at ic_input_redshift for the cosmology in parameters (the rest from
    parameterfile.lua) to path, as columns k (h/Mpc) and P(k) (Mpc/h)^3. Nothing is done if path exists. The CLASS
    input file is kept next to it as a record."""
    if os.path.exists(path):
        return
    if parameters.get("cosmology_model", base_parameter("cosmology_model")) != "LCDM":
        raise SystemExit("The input P(k) is only set up for cosmology_model = LCDM")

    def value(name):
        return float(parameters.get(name, base_parameter(name)))

    h = value("cosmology_h")
    m_nu = 93.14 * value("cosmology_OmegaMNu") * h**2  # Sum of neutrino masses in eV
    if m_nu > 0:  # One massive neutrino, the rest of N_eff massless
        neutrinos = f"N_ncdm = 1\nm_ncdm = {m_nu}\nN_ur = {value('cosmology_Neffective') - 1.0132}"
    else:
        neutrinos = f"N_ur = {value('cosmology_Neffective')}"
    root = path.removesuffix(".txt") + "_class_"
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(root + "input.ini", "w") as f:
        f.write(f"""output = mPk
z_pk = {value('ic_input_redshift')}
P_k_max_h/Mpc = 100
k_per_decade_for_pk = 32
h = {h}
omega_b = {value('cosmology_Omegab') * h**2}
omega_cdm = {value('cosmology_OmegaCDM') * h**2}
Omega_k = {value('cosmology_OmegaK')}
T_cmb = {value('cosmology_TCMB_kelvin')}
A_s = {value('cosmology_As')}
n_s = {value('cosmology_ns')}
k_pivot = {value('cosmology_kpivot_mpc')}
{neutrinos}
root = {root}
""")
    subprocess.run([CLASS, root + "input.ini"], check=True, stdout=subprocess.DEVNULL)
    k, P = np.loadtxt(root + ("00_pk_cb.dat" if m_nu > 0 else "00_pk.dat"), unpack=True)  # CDM + baryons
    np.savetxt(path, np.column_stack([k, P]), header="k (h/Mpc)   P(k) (Mpc/h)^3, linear CDM+baryon from CLASS")
    for output in set(glob.glob(root + "*")) - {root + "input.ini"}:
        os.remove(output)
    print(f"Computed the input P(k) with CLASS: {path}")
