"""The linear input P(k) for the simulations, when a scan changes the cosmology (any cosmology_* parameter).

write_input_power(cosmology, z, path) writes the CDM+baryon P(k) at redshift z as columns k (h/Mpc) and P(k) (Mpc/h)^3,
like ../COLASolver/input/example_power_spectrum_cb_z0.000.txt. It uses CAMB if installed (pip install --user camb),
otherwise CLASS (classy). cosmology holds the cosmology_* values of parameterfile.lua (as floats, without the prefix).
"""
import numpy as np

KMIN, KMAX, NK = 1e-4, 100.0, 700  # k range and number of points of the output, in h/Mpc


def camb_power(c, z):
    import camb
    h = c["h"]
    parameters = camb.set_params(
        H0=100 * h, ombh2=c["Omegab"] * h**2, omch2=c["OmegaCDM"] * h**2, omk=c["OmegaK"],
        mnu=93.14 * c["OmegaMNu"] * h**2, nnu=c["Neffective"], TCMB=c["TCMB_kelvin"],
        As=c["As"], ns=c["ns"], pivot_scalar=c["kpivot_mpc"], WantTransfer=True, redshifts=[z], kmax=1.5 * KMAX * h)
    k, _, P = camb.get_results(parameters).get_matter_power_spectrum(
        minkh=KMIN, maxkh=KMAX, npoints=NK, var1="delta_nonu", var2="delta_nonu")  # delta_nonu: CDM + baryons
    return k, P[0]


def class_power(c, z):
    from classy import Class
    h = c["h"]
    parameters = {"output": "mPk", "P_k_max_1/Mpc": 1.5 * KMAX * h, "z_pk": z, "h": h,
                  "omega_b": c["Omegab"] * h**2, "omega_cdm": c["OmegaCDM"] * h**2, "Omega_k": c["OmegaK"],
                  "T_cmb": c["TCMB_kelvin"], "A_s": c["As"], "n_s": c["ns"], "k_pivot": c["kpivot_mpc"]}
    if c["OmegaMNu"] > 0:  # One massive neutrino, the rest of N_eff massless
        parameters.update({"N_ncdm": 1, "m_ncdm": 93.14 * c["OmegaMNu"] * h**2, "N_ur": c["Neffective"] - 1.0132})
    else:
        parameters["N_ur"] = c["Neffective"]
    cosmo = Class()
    cosmo.set(parameters)
    cosmo.compute()
    k = np.logspace(np.log10(KMIN), np.log10(KMAX), NK)  # h/Mpc
    pk = cosmo.pk_cb_lin if c["OmegaMNu"] > 0 else cosmo.pk_lin  # Without massive neutrinos, CDM+baryon is all matter
    P = np.array([pk(ki * h, z) for ki in k]) * h**3  # CLASS works in 1/Mpc and Mpc^3
    cosmo.struct_cleanup()
    return k, P


def write_input_power(cosmology, z, path):
    try:
        import camb  # noqa: F401
        k, P, code = *camb_power(cosmology, z), "CAMB"
    except ImportError:
        try:
            import classy  # noqa: F401
        except Exception as error:  # Not installed, or e.g. built against another numpy version
            raise SystemExit(f"Changing cosmology_* parameters needs CAMB (pip install --user camb) or CLASS (classy) to "
                             f"compute the input P(k). CAMB is not installed, and classy fails with: {error!r}")
        k, P, code = *class_power(cosmology, z), "CLASS"
    np.savetxt(path, np.column_stack([k, P]), header=f"k (h/Mpc)   P(k) (Mpc/h)^3, linear CDM+baryon at z = {z}, {code}")
    print(f"Computed the input P(k) with {code}: {path}")
