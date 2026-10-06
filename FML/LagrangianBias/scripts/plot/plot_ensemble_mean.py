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
        return None, None, None, 0
        
    # Load all realizations into a 2D matrix
    data_matrix = np.array([np.loadtxt(f)[:, 1] for f in file_list])
    k = np.loadtxt(file_list[0])[:, 0]
    
    # Average first to cancel zero-crossing noise, then take absolute value if needed
    mean_val = np.mean(data_matrix, axis=0)
    
    if is_cross_spectrum:
        mean_val = np.abs(mean_val)
        data_matrix = np.abs(data_matrix)
        
    # Mask out non-positive values to prevent log scale crashes
    valid = mean_val > 1e-10
    
    return k[valid], mean_val[valid], data_matrix[:, valid], num_files

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

# Loop through and plot all 15 fields
for (field, label, apply_abs, linestyle), color in zip(plot_configs, colors):
    k, mean_val, all_data, n_files = load_and_average(field, is_cross_spectrum=apply_abs)
    
    if k is None:
        print(f"Skipping P_{field}: No files found.")
        continue
        
    max_n = max(max_n, n_files)
    
    # Optional: Plot the raw stochastic background lines
    if args.show_all:
        for i in range(all_data.shape[0]):
            ax.plot(k, all_data[i], color=color, alpha=0.05, lw=0.5, linestyle=linestyle)
            
    # Plot the clean, averaged trendline
    # Increased lw slightly for dotted lines (':') so they stand out clearly
    lw = 2.5 if linestyle == ':' else 2.0
    ax.plot(k, mean_val, color=color, lw=lw, linestyle=linestyle, label=label)

# Set up axes identically to the Gnuplot layout
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel(r'Wavenumber $k \ [h/\mathrm{Mpc}]$')
ax.set_ylabel(r'Power $P_{ij}(k) \ [(\mathrm{Mpc}/h)^3]$')
ax.set_title(f"Ensemble Averaged Power Spectra ($N={max_n}$)", pad=15)

# Place the legend exactly where it was in Gnuplot (top right, outside the box)
ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left", framealpha=0.9, edgecolor='black')

ax.grid(True, which="both", ls=":", alpha=0.4)

plt.tight_layout()
output_path = os.path.join(OUTPUT_DIR, "ensemble_spectra_1to1.pdf")
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"Plot successfully saved as {output_path}")
plt.show()
