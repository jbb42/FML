import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import argparse
import matplotlib.cm as cm

# Set up command-line arguments
parser = argparse.ArgumentParser(description="Plot all 15 averaged FML COLA power spectra.")
parser.add_argument('--show-all', action='store_true', help="Plot all individual realizations in the background.")
args = parser.parse_args()

# LagrangianBias root (this script lives in LagrangianBias/scripts/plot)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(ROOT, "results", "ensemble_100seeds")
OUTPUT_DIR = os.path.join(ROOT, "figures", "ensemble_100seeds")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Enable LaTeX rendering and set publication-quality fonts
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 12,
    "axes.labelsize": 14,
    "legend.fontsize": 10
})

def load_and_average(field, is_cross_spectrum=False):
    file_list = glob.glob(os.path.join(DATA_DIR, "output_seed_*", f"pofk_{field}.txt"))
    num_files = len(file_list)
    
    if num_files == 0:
        return None, None, None, None, 0
        
    # Load all realizations into a 2D matrix
    data_matrix = np.array([np.loadtxt(f)[:, 1] for f in file_list])
    k = np.loadtxt(file_list[0])[:, 0]
    
    # Calculate both mean and standard deviation from the raw data
    mean_val = np.mean(data_matrix, axis=0)
    std_val = np.std(data_matrix, axis=0)
    
    if is_cross_spectrum:
        mean_val = np.abs(mean_val)
        data_matrix = np.abs(data_matrix)
        
    # Mask out non-positive values to prevent log scale crashes
    valid = mean_val > 1e-10
    
    return k[valid], mean_val[valid], std_val[valid], data_matrix[:, valid], num_files

# Define the exact 1:1 mapping from the Gnuplot script
# Format: (field_id, label, apply_abs, linestyle)
plot_configs = [
    ('00', r'$P_{00}$ (Matter Auto)', False, '-'),
    ('01', r'$P_{01}$ (Matter - $\delta_L$)', False, '-'),
    ('11', r'$P_{11}$ ($\delta_L$ Auto)', False, '-'),
    ('02', r'$P_{02}$ (Matter - $\delta_L^2$)', False, '-'),
    ('12', r'$|P_{12}|$ ($\delta_L - \delta_L^2$)', True, '-'),
    ('22', r'$P_{22}$ ($\delta_L^2$ Auto)', False, '-'),
    ('03', r'$|P_{03}|$ (Matter - $s^2$)', True, '--'),
    ('13', r'$|P_{13}|$ ($\delta_L - s^2$)', True, '--'),
    ('23', r'$P_{23}$ ($\delta_L^2 - s^2$)', False, '--'),
    ('33', r'$P_{33}$ ($s^2$ Auto)', False, '--'),
    ('04', r'$|P_{04}|$ (Matter - $\nabla^2\delta_L$)', True, ':'),
    ('14', r'$|P_{14}|$ ($\delta_L - \nabla^2\delta_L$)', True, ':'),
    ('24', r'$|P_{24}|$ ($\delta_L^2 - \nabla^2\delta_L$)', True, ':'),
    ('34', r'$|P_{34}|$ ($s^2 - \nabla^2\delta_L$)', True, ':'),
    ('44', r'$P_{44}$ ($\nabla^2\delta_L$ Auto)', False, ':'),
]

fig, ax = plt.subplots(figsize=(12, 8))

# Generate 15 visually distinct colors using a colormap
colors = cm.tab20(np.linspace(0, 1, len(plot_configs)))

max_n = 0
k_min_global = np.inf
k_max_global = -np.inf

# Loop through and plot all 15 fields
for (field, label, apply_abs, linestyle), color in zip(plot_configs, colors):
    k, mean_val, std_val, all_data, n_files = load_and_average(field, is_cross_spectrum=apply_abs)
    
    if k is None or len(k) == 0:
        print(f"Skipping P_{field}: No files found.")
        continue
        
    max_n = max(max_n, n_files)
    k_min_global = min(k_min_global, k.min())
    k_max_global = max(k_max_global, k.max())
    
    # Optional: Plot the raw stochastic background lines
    if args.show_all:
        for i in range(all_data.shape[0]):
            ax.plot(k, all_data[i], color=color, alpha=0.05, lw=0.5, linestyle=linestyle)
            
    # Plot the clean, averaged trendline
    lw = 2.5 if linestyle == ':' else 2.0
    ax.plot(k, mean_val, color=color, lw=lw, linestyle=linestyle, label=label)
    
    # Plot the shaded variance band
    lower_bound = np.maximum(mean_val - std_val, 1e-10)
    upper_bound = mean_val + std_val
    ax.fill_between(k, lower_bound, upper_bound, color=color, alpha=0.15)

# Set up axes identically to the Gnuplot layout
ax.set_xscale('log')
ax.set_yscale('log')

# Apply tight x-limits to crop to the data and hardcoded y-limits
if k_min_global < np.inf:
    ax.set_xlim([k_min_global, k_max_global])
ax.set_ylim([1e-2, 3e4])

ax.set_xlabel(r'Wavenumber $k \ [h/\mathrm{Mpc}]$')
ax.set_ylabel(r'Power $P_{ij}(k) \ [(\mathrm{Mpc}/h)^3]$')
ax.set_title(f"Ensemble Averaged Power Spectra with $1\sigma$ Variance ($N={max_n}$)", pad=15)

# Place the legend exactly where it was in Gnuplot (top right, outside the box)
ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left", framealpha=0.9, edgecolor='black')

ax.grid(True, which="both", ls=":", alpha=0.4)

plt.tight_layout()
output_path = os.path.join(OUTPUT_DIR, "ensemble_spectra_variance.pdf")
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"Plot successfully saved as {output_path}")
plt.show()
