"""The linear input P(k) for scans that change cosmology_* parameters, computed by running the CLASS executable."""
import glob
import os
import subprocess

import numpy as np

from simulations import base_parameter

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


def neutrinos(Neffective, OmegaMNu):
    """CLASS neutrino parameters for COLASolver's neutrinos: 3 species with temperature
    T_CMB (N_eff / 3)^(1/4) (4/11)^(1/3), sharing the mass that gives OmegaMNu."""
    if OmegaMNu == 0:
        return {"N_ur": Neffective}
    return {"N_ur": 0, "N_ncdm": 1, "deg_ncdm": 3, "Omega_ncdm": OmegaMNu,
            "T_ncdm": (Neffective / 3) ** 0.25 * (4 / 11) ** (1 / 3)}


def write_input_power(parameters, path):
    """Write the linear CDM+baryon P(k) at ic_input_redshift for the cosmology in parameters (the rest from
    parameterfile.lua) to path, as columns k (h/Mpc) and P(k) (Mpc/h)^3. Nothing is done if path exists. The CLASS
    input file is kept next to it as a record."""
    if os.path.exists(path):
        return

    def value(name):
        return parameters.get(name, base_parameter(name))

    unknown = [p for p in parameters if p.startswith("cosmology_") and p not in LUA_TO_CLASS
               and p not in ("cosmology_Neffective", "cosmology_OmegaMNu")]
    if value("cosmology_model") != "LCDM" or unknown:
        raise SystemExit(f"The input P(k) is only set up for cosmology_model = LCDM and {list(LUA_TO_CLASS)}, "
                         f"cosmology_Neffective and cosmology_OmegaMNu, not {unknown or value('cosmology_model')}")

    class_parameters = {class_name: value(lua_name) for lua_name, class_name in LUA_TO_CLASS.items()}
    class_parameters.update(neutrinos(float(value("cosmology_Neffective")), float(value("cosmology_OmegaMNu"))))
    root = path.removesuffix(".txt") + "_class_"
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(root + "input.ini", "w") as f:
        f.write("output = mPk\nP_k_max_h/Mpc = 100\nk_per_decade_for_pk = 32\n")
        f.writelines(f"{name} = {v}\n" for name, v in class_parameters.items())
        f.write(f"root = {root}\n")
    subprocess.run([CLASS, root + "input.ini"], check=True, stdout=subprocess.DEVNULL)

    massive_neutrinos = float(value("cosmology_OmegaMNu")) > 0
    k, P = np.loadtxt(root + ("00_pk_cb.dat" if massive_neutrinos else "00_pk.dat"), unpack=True)  # CDM + baryons
    np.savetxt(path, np.column_stack([k, P]), header="k (h/Mpc)   P(k) (Mpc/h)^3, linear CDM+baryon from CLASS")
    for output in set(glob.glob(root + "*")) - {root + "input.ini"}:
        os.remove(output)
    print(f"Computed the input P(k) with CLASS: {path}")
