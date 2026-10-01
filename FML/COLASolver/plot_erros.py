import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# Define directories
DATA_DIR = "saved_spectra2"
OUTPUT_DIR = "error_plots2"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==========================================
# TOGGLE FOR 256 RESOLUTION PARAMETERS
# ==========================================
INCLUDE_256 = False  # Set to True to include 256, False to hide them completely
# ==========================================

# Enable LaTeX rendering and set publication-quality fonts (larger sizes for readability)
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 13,
    "axes.labelsize": 16,
    "axes.titlesize": 16,
    "legend.fontsize": 9,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12
})

def load_and_average_params(f, n, s, field, is_cross_spectrum=False):
    """Loads all versions for a specific parameter combination and computes mean & std."""
    file_list = glob.glob(os.path.join(DATA_DIR, f"{f}_{n}_{s}_v*", f"pofk_{field}.txt"))
    num_files = len(file_list)
    
    if num_files == 0:
        return None, None, None, 0
        
    data_matrix = np.array([np.loadtxt(fn)[:, 1] for fn in file_list])
    k = np.loadtxt(file_list[0])[:, 0]
    
    mean_val = np.mean(data_matrix, axis=0)
    std_val = np.std(data_matrix, axis=0)
    
    if is_cross_spectrum:
        mean_val = np.abs(mean_val)
        
    valid = mean_val > 1e-10
    return k[valid], mean_val[valid], std_val[valid], num_files

# Define the 15 spectra mapping
plot_configs = [
    ('00', r'$P_{00}$ (Matter Auto)', False),
    ('01', r'$P_{01}$ (Matter - $\delta_L$)', True),
    ('11', r'$P_{11}$ ($\delta_L$ Auto)', False),
    ('02', r'$P_{02}$ (Matter - $\delta_L^2$)', True),
    ('12', r'$|P_{12}|$ ($\delta_L - \delta_L^2$)', True),
    ('22', r'$P_{22}$ ($\delta_L^2$ Auto)', False),
    ('03', r'$|P_{03}|$ (Matter - $s^2$)', True),
    ('13', r'$|P_{13}|$ ($\delta_L - s^2$)', True),
    ('23', r'$|P_{23}|$ ($\delta_L^2 - s^2$)', True),
    ('33', r'$P_{33}$ ($s^2$ Auto)', False),
    ('04', r'$|P_{04}|$ (Matter - $\nabla^2\delta_L$)', True),
    ('14', r'$|P_{14}|$ ($\delta_L - \nabla^2\delta_L$)', True),
    ('24', r'$|P_{24}|$ ($\delta_L^2 - \nabla^2\delta_L$)', True),
    ('34', r'$|P_{34}|$ ($s^2 - \nabla^2\delta_L$)', True),
    ('44', r'$P_{44}$ ($\nabla^2\delta_L$ Auto)', False),
]

# Discover parameter combinations (updated for 4-part name: f_n_s_v)
param_combos = set()
for d in glob.glob(os.path.join(DATA_DIR, "f*_n*_s*_v*")):
    dirname = os.path.basename(d)
    parts = dirname.split('_')
    if len(parts) >= 4:
        param_combos.add((parts[0], parts[1], parts[2]))

# Filter based on toggle
if not INCLUDE_256:
    param_combos = {c for c in param_combos if not (c[0] == 'f256' or c[1] == 'n256')}

def numerical_sort_key(combo):
    f, n, s = combo
    return (int(f[1:]), int(n[1:]), int(s[1:]))

param_combos = sorted(list(param_combos), key=numerical_sort_key)
print(f"Found {len(param_combos)} unique parameter combinations (INCLUDE_256={INCLUDE_256}) in '{DATA_DIR}'.")

# Best reference parameter configuration (highest resolution at the end of sorted list)
best_combo = param_combos[-1]
print(f"Reference (best) parameter configuration set to: {best_combo}")

# Setup standard matplotlib color order grouped by remaining parameters (n, s)
unique_remaining = sorted(list(set((n, s) for f, n, s in param_combos)))
color_mapping = {rem: f"C{i % 10}" for i, rem in enumerate(unique_remaining)}

# Loop through each spectrum and generate large, clear error plots
for field, spectrum_label, apply_abs in plot_configs:
    print(f"Generating error plot for P_{field}: {spectrum_label}")
    
    # Load reference data
    k_ref, mean_ref, _, _ = load_and_average_params(*best_combo, field, is_cross_spectrum=apply_abs)
    if k_ref is None:
        print(f"  Warning: Reference data not found for P_{field}, skipping.")
        continue
        
    # Create a large figure for maximum legibility
    fig, ax = plt.subplots(figsize=(14, 9))
    plotted_any = False
    
    for idx, (f, n, s) in enumerate(param_combos):
        k, mean_val, std_val, n_files = load_and_average_params(
            f, n, s, field, is_cross_spectrum=apply_abs
        )
        
        if k is None or len(k) == 0:
            continue
            
        plotted_any = True
        
        # Line style mapping based on f resolution
        if f == 'f512':
            linestyle = '--'
        elif f == 'f256':
            linestyle = ':'
        else:
            linestyle = '-'
            
        color = color_mapping[(n, s)]
        param_label = f"$f={f[1:]}, n={n[1:]}, s={s[1:]}$"
        
        # Calculate fractional error: (P - P_best) / P_best
        if len(k) == len(k_ref) and np.allclose(k, k_ref):
            frac_error = (mean_val - mean_ref) / mean_ref
            ax.plot(k, frac_error, color=color, linestyle=linestyle, lw=1.0, label=param_label)
            
    if not plotted_any:
        plt.close(fig)
        continue
        
    # Configure axes for detailed inspection
    ax.set_xscale('log')
    ax.set_xlim([k_ref.min(), k_ref.max()])
    ax.set_ylim([-0.4, 0.4])  # Clean window for fractional error (+/- 40%)
    
    ax.set_xlabel(r'Wavenumber $k \ [h/\mathrm{Mpc}]$')
    ax.set_ylabel(r'Fractional Error $(P - P_{\mathrm{best}}) / P_{\mathrm{best}}$')
    ax.set_title(f"Fractional Error relative to Best ({best_combo[0]}, {best_combo[1]}, {best_combo[2]}): {spectrum_label}", pad=15)
    
    # Spacious legend outside the plotting area
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", framealpha=0.95, edgecolor='black', fontsize=8)
    ax.grid(True, which="both", ls=":", alpha=0.5)
    
    plt.tight_layout()
    output_path = os.path.join(OUTPUT_DIR, f"error_P_{field}.pdf")
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    print(f"  Saved large error plot to {output_path}")
    plt.close(fig)

print(f"All large error plots successfully saved in '{OUTPUT_DIR}/'!")
