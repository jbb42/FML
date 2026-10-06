"""The linear input P(k) for scans that change cosmology_* parameters, computed by running the CLASS executable."""
import os
import subprocess
import tempfile

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
# The same for cosmology_model = "w0waCDM", where a fluid with w(a) = w0 + wa (1 - a) replaces the cosmological constant
LUA_TO_CLASS_W0WA = {**LUA_TO_CLASS, "cosmology_w0": "w0_fld", "cosmology_wa": "wa_fld"}


def neutrinos(Neffective, OmegaMNu):
    """CLASS neutrino parameters for COLASolver's neutrinos: 3 species with temperature
    T_CMB (N_eff / 3)^(1/4) (4/11)^(1/3), sharing the mass that gives OmegaMNu."""
    if OmegaMNu == 0:
        return {"N_ur": Neffective}
    return {"N_ur": 0, "N_ncdm": 1, "deg_ncdm": 3, "Omega_ncdm": OmegaMNu,
            "T_ncdm": (Neffective / 3) ** 0.25 * (4 / 11) ** (1 / 3)}


def write_input_power(parameters, path, kmax):
    """Write the linear CDM+baryon P(k) at ic_input_redshift for the cosmology in parameters (the rest from
    parameterfile.lua) to path, as columns k (h/Mpc) and P(k) (Mpc/h)^3, up to at least kmax (h/Mpc). Nothing is done if
    path already goes that far. The CLASS input is kept next to it, in <path>_class.ini."""
    if os.path.exists(path) and np.loadtxt(path)[-1, 0] >= kmax:
        return

    def value(name):
        return parameters.get(name, base_parameter(name))

    model = value("cosmology_model")
    mapping = {"LCDM": LUA_TO_CLASS, "w0waCDM": LUA_TO_CLASS_W0WA}.get(model)
    unknown = [p for p in parameters if p.startswith("cosmology_") and p not in (mapping or {})
               and p not in ("cosmology_model", "cosmology_Neffective", "cosmology_OmegaMNu")]
    if mapping is None or unknown:
        raise SystemExit(f"The input P(k) is set up for cosmology_model = LCDM or w0waCDM, not {model}, "
                         f"and for the parameters in input_power.py, not {unknown}")

    class_parameters = {class_name: value(lua_name) for lua_name, class_name in mapping.items()}
    if model == "w0waCDM":
        class_parameters["Omega_Lambda"] = 0  # The fluid then fills the rest of the energy budget
    class_parameters.update(neutrinos(float(value("cosmology_Neffective")), float(value("cosmology_OmegaMNu"))))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:  # CLASS only accepts short output paths
        ini = "output = mPk\n" + f"P_k_max_h/Mpc = {max(100.0, 1.1 * kmax)}\nk_per_decade_for_pk = 32\n"
        ini += "".join(f"{name} = {v}\n" for name, v in class_parameters.items())
        with open(os.path.join(tmp, "input.ini"), "w") as f:
            f.write(ini + f"root = {tmp}/\n")
        result = subprocess.run([CLASS, os.path.join(tmp, "input.ini")], capture_output=True, text=True)
        if result.returncode != 0:
            raise SystemExit(f"CLASS failed for {path}:\n{ini}{result.stdout[-2000:]}{result.stderr[-2000:]}")
        massive_neutrinos = float(value("cosmology_OmegaMNu")) > 0
        k, P = np.loadtxt(os.path.join(tmp, "00_pk_cb.dat" if massive_neutrinos else "00_pk.dat"), unpack=True)  # CDM + baryons
    np.savetxt(path, np.column_stack([k, P]), header="k (h/Mpc)   P(k) (Mpc/h)^3, linear CDM+baryon from CLASS")
    with open(path.removesuffix(".txt") + "_class.ini", "w") as f:  # The CLASS input, as a record
        f.write(ini)
    print(f"Computed the input P(k) with CLASS: {path}")
