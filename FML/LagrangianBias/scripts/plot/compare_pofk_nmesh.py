#!/usr/bin/env python3
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------------------
# Directory Configurations
# -------------------------------------------------------------
# LagrangianBias root (this script lives in LagrangianBias/scripts/plot)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DIR_V4 = os.path.join(ROOT, "results", "convergence_runs_v1", "f1024_n1024_s40_v4")
DIR_NEW = os.path.join(ROOT, "output", "snapshot_TestSim_z0.000")

FIELD_NAMES = {
    0: r"\mathrm{Matter}",
    1: r"\delta_L",
    2: r"\delta_L^2",
    3: r"s^2",
    4: r"\nabla^2 \delta_L",
}

SPECTRA_PAIRS = [
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4),
    (1, 1), (1, 2), (1, 3), (1, 4),
    (2, 2), (2, 3), (2, 4),
    (3, 3), (3, 4),
    (4, 4),
]

def load_pofk(folder, i, j):
    filename = f"pofk_{i}{j}.txt"
    path = os.path.join(folder, filename)
    if not os.path.exists(path):
        filename = f"pofk_{j}{i}.txt"
        path = os.path.join(folder, filename)
    if not os.path.exists(path):
        return None, None
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1]

print(f"Comparing:")
print(f"  v4  (pofk_nmesh = 512):  {DIR_V4}")
print(f"  New (pofk_nmesh = 1024): {DIR_NEW}")
print("Close each figure window to step to the next spectrum.\n")

for i, j in SPECTRA_PAIRS:
    k_v4, p_v4 = load_pofk(DIR_V4, i, j)
    k_new, p_new = load_pofk(DIR_NEW, i, j)

    if k_v4 is None or k_new is None:
        print(f"Skipping P_{i}{j}: missing file.")
        continue

    # Interpolate to compare over overlapping k-bins
    overlap = (k_new >= k_v4.min()) & (k_new <= k_v4.max())
    k_eval = k_new[overlap]
    p_new_eval = p_new[overlap]
    p_v4_interp = np.interp(k_eval, k_v4, p_v4)

    with np.errstate(divide="ignore", invalid="ignore"):
        frac_diff = (p_new_eval - p_v4_interp) / p_v4_interp

    name_i = FIELD_NAMES.get(i, f"b_{i}")
    name_j = FIELD_NAMES.get(j, f"b_{j}")
    label_title = rf"$P_{{{i}{j}}}(k): {name_i} \times {name_j}$"

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(8, 7), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]}
    )

    # Upper panel: Absolute spectra
    ax1.loglog(k_v4, np.abs(p_v4), label=r"$v4$ ($\mathrm{pofk}=512$)", color="crimson", ls="--", lw=1.8)
    ax1.loglog(k_new, np.abs(p_new), label=r"$\mathrm{New}$ ($\mathrm{pofk}=1024$)", color="navy", lw=1.8)
    ax1.set_ylabel(rf"$|P_{{{i}{j}}}(k)| \ [(\mathrm{{Mpc}}/h)^3]$", fontsize=11)
    ax1.set_title(f"Comparison: {label_title}", fontsize=12)
    ax1.grid(True, which="both", ls=":", alpha=0.5)
    ax1.legend(loc="best", fontsize=10)

    # Lower panel: Fractional error
    ax2.plot(k_eval, frac_diff, color="forestgreen", lw=1.6)
    ax2.axhline(0.0, color="black", lw=0.9, ls="--")
    ax2.axhline(0.01, color="gray", lw=0.7, ls=":")
    ax2.axhline(-0.01, color="gray", lw=0.7, ls=":")
    ax2.axvline(1.0, color="firebrick", lw=0.9, ls=":", label=r"$k = 1.0$")
    ax2.set_xscale("log")
    ax2.set_ylim(-0.08, 0.08)
    ax2.set_ylabel(r"$\frac{P_{1024} - P_{512}}{P_{512}}$", fontsize=11)
    ax2.set_xlabel(r"Wavenumber $k \ [h/\mathrm{Mpc}]$", fontsize=11)
    ax2.grid(True, which="both", ls=":", alpha=0.5)
    ax2.legend(loc="upper left", fontsize=9)

    plt.tight_layout()
    plt.show()

print("Finished all spectra.")
