import numpy as np
import matplotlib.pyplot as plt
import glob
import matplotlib.cm as cm
import os

# Define the directory where your simulation folders are located
# COLASolver root (this script lives in COLASolver/scripts/plot)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(ROOT, "results", "convergence_runs_v2_fixedseeds")
OUTPUT_DIR = os.path.join(ROOT, "figures", "convergence_spectra")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==========================================
# TOGGLE FOR 256 RESOLUTION PARAMETERS
# ==========================================
INCLUDE_256 = False  # Set to True to include 256, False to hide them
# ==========================================

# Enable LaTeX rendering and set publication-quality fonts
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 12,
    "axes.labelsize": 14,
    "legend.fontsize": 9
})

def load_and_average_params(f, n, s, field, is_cross_spectrum=False):
    """Loads all versions for a specific parameter combination and computes mean & std."""
    file_list = glob.glob(os.path.join(DATA_DIR, f"{f}_{n}_{s}_v*", f"pofk_{field}.txt"))
    if not file_list:
        return None, None, None, 0
        
    data_matrix = np.array([np.loadtxt(fn)[:, 1] for fn in file_list])
    k = np.loadtxt(file_list[0])[:, 0]
    
    mean_val = np.mean(data_matrix, axis=0)
    std_val = np.std(data_matrix, axis=0)
    
    if is_cross_spectrum:
        mean_val = np.abs(mean_val)
        
    valid = mean_val > 1e-10
    return k[valid], mean_val[valid], std_val[valid], len(file_list)

# Define the 15 spectra mapping
plot_configs = [
    ('00', r'$P_{00}$ (Matter Auto)', False),
    ('01', r'$P_{01}$ (Matter - $\delta_L$)', False),
    ('11', r'$P_{11}$ ($\delta_L$ Auto)', False),
    ('02', r'$P_{02}$ (Matter - $\delta_L^2$)', False),
    ('12', r'$|P_{12}|$ ($\delta_L - \delta_L^2$)', True),
    ('22', r'$P_{22}$ ($\delta_L^2$ Auto)', False),
    ('03', r'$|P_{03}|$ (Matter - $s^2$)', True),
    ('13', r'$|P_{13}|$ ($\delta_L - s^2$)', True),
    ('23', r'$P_{23}$ ($\delta_L^2 - s^2$)', False),
    ('33', r'$P_{33}$ ($s^2$ Auto)', False),
    ('04', r'$|P_{04}|$ (Matter - $\nabla^2\delta_L$)', True),
    ('14', r'$|P_{14}|$ ($\delta_L - \nabla^2\delta_L$)', True),
    ('24', r'$|P_{24}|$ ($\delta_L^2 - \nabla^2\delta_L$)', True),
    ('34', r'$|P_{34}|$ ($s^2 - \nabla^2\delta_L$)', True),
    ('44', r'$P_{44}$ ($\nabla^2\delta_L$ Auto)', False),
]

# Discover unique parameter combinations
param_combos = set()
for d in glob.glob(os.path.join(DATA_DIR, "f*_n*_s*_v*")):
    parts = os.path.basename(d).split('_')
    # Folder format now f_n_s_v (length 4)
    if len(parts) >= 4:
        param_combos.add((parts[0], parts[1], parts[2]))

# Filter based on toggle
if not INCLUDE_256:
    param_combos = {c for c in param_combos if not (c[0] == 'f256' or c[1] == 'n256')}

# Sort numerically
param_combos = sorted(list(param_combos), key=lambda x: (int(x[0][1:]), int(x[1][1:]), int(x[2][1:])))
print(f"Found {len(param_combos)} unique parameter combinations (INCLUDE_256={INCLUDE_256}).")

# Standard Matplotlib cycle colors (C0 = blue, C1 = orange, C2 = green, etc.)
unique_remaining = sorted(list(set((n, s) for f, n, s in param_combos)))
color_mapping = {rem: f"C{i % 10}" for i, rem in enumerate(unique_remaining)}

# Loop through each spectrum
for field, spectrum_label, apply_abs in plot_configs:
    print(f"Processing spectrum P_{field}: {spectrum_label}")
    
    fig, ax = plt.subplots(figsize=(10, 7))
    k_min_global, k_max_global, plotted_any = np.inf, -np.inf, False
    
    for f, n, s in param_combos:
        k, mean_val, std_val, _ = load_and_average_params(f, n, s, field, is_cross_spectrum=apply_abs)
        
        if k is None or len(k) == 0:
            continue
            
        plotted_any = True
        k_min_global = min(k_min_global, k.min())
        k_max_global = max(k_max_global, k.max())
        
        # Line style based on f resolution
        linestyle = '--' if f == 'f512' else (':' if f == 'f256' else '-')
        
        # Fixed mapping call (removed 'b')
        color = color_mapping[(n, s)]
        param_label = f"$f={f[1:]}, n={n[1:]}, s={s[1:]}$"
        
        ax.plot(k, mean_val, color=color, linestyle=linestyle, lw=0.5, label=param_label)
        ax.fill_between(k, np.maximum(mean_val - std_val, 1e-10), mean_val + std_val, color=color, alpha=0.01)
        
    if not plotted_any:
        print(f"  No data found for P_{field}, skipping.")
        plt.close(fig)
        continue
        
    ax.set_xscale('log')
    ax.set_yscale('log')
    if k_min_global < np.inf:
        ax.set_xlim([k_min_global, k_max_global])
        
    ax.set_ylim([1e-2, 3e4])
    
    ax.set_xlabel(r'Wavenumber $k \ [h/\mathrm{Mpc}]$')
    ax.set_ylabel(r'Power $P_{ij}(k) \ [(\mathrm{Mpc}/h)^3]$')
    ax.set_title(f"Spectrum {spectrum_label} ($N=10$ versions)", pad=15)
    
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left", framealpha=0.9, edgecolor='black', fontsize=8)
    ax.grid(True, which="both", ls=":", alpha=0.4)
    
    plt.tight_layout()
    output_filename = os.path.join(OUTPUT_DIR, f"spectrum_P_{field}.pdf")
    plt.savefig(output_filename, dpi=600, bbox_inches='tight')
    print(f"  Saved as {output_filename}")
    plt.close(fig)

print("All spectra successfully plotted and saved!")

