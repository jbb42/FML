# Lagrangian bias spectra

Basis spectra of the hybrid Lagrangian bias expansion from COLA simulations: the matter field and the
weights {$\delta_L$, delta_L^2 - <delta_L^2>, s^2 - <s^2>, nabla^2 delta_L}, evaluated at the particles'
Lagrangian positions q, advected with the particles, and all 15 auto and cross power spectra `pofk_ij.txt`.

## Code

`LagrangianBias.h` (namespace `FML::LAGRANGIANBIAS`) does all the work and is independent of the solver:

| Function | What it does |
|---|---|
| `compute_bias_weights` | Builds the weights at q from the initial field, scaled per mode with D(k, a_out) / D(k, a_ini) |
| `compute_bias_power_spectra` | Deposits the 5 fields at the Eulerian positions and writes `pofk_ij.txt` and `pofk_bias_info.txt` |
| `zero_bias_weights` | Zeros the weights (particles read from file, no linear field) |

`COLASolver` calls it from three places:
- `src/Main.cpp`: the `Particle` carries `bias_weights[]`, `active_bias_index` and `get_mass()`
- `src/Simulation.h`: keeps the initial field (`bias_delta_ini_fourier`) and calls `compute_bias_weights` at each output
- `src/AnalyzeOutput.h`: calls `compute_bias_power_spectra` when `pofk = true`

Build the solver as usual in `../COLASolver` (`make`). The scripts here run `../COLASolver/nbody`.

## Running

All scripts locate the LagrangianBias root from their own path, so they can be run from anywhere.
Simulations are run with the root as working directory, which the run scripts handle by `cd`-ing there,
so the input file path `../COLASolver/input/...` in the parameter files is relative to the root.

## Folder layout (LagrangianBias/)

| Folder | Contents | In git |
|---|---|---|
| `LagrangianBias.h` | The bias weights and spectra (see above) | yes |
| `parameterfile.lua` | Base parameter file used by all run scripts | yes |
| `params/` | Extra parameter files (`parameterfile_2048.lua`: 2048^3 reference run) | yes |
| `scripts/run/` | Batch scripts that run simulations | yes |
| `scripts/plot/` | Plotting and comparison scripts | yes |
| `results/` | Saved `pofk_ij.txt` spectra, one subfolder per campaign | READMEs only |
| `figures/` | Generated plots | no |
| `logs/` | Logs from batch runs | no |
| `output/` | Default output folder of `parameterfile.lua` | no |
| `hostlist.txt` | Nodes for MPI runs | yes |

Parameter files saved next to older results still have the old `COLASolver/` paths. They are only a record
of the settings used (`plot_fofr_boost.py` reads them for that) and are not meant to be rerun as they are.

## Run scripts (`scripts/run/`)

| Script | What it does | Writes to |
|---|---|---|
| `run_ensemble_100.sh` | 100 seeds of `parameterfile.lua` (seed = 12300 + i) | `results/ensemble_100seeds/output_seed_i/` |
| `run_convergence_grid.sh` | force mesh x particles x time steps grid, 10 seeds each | `results/convergence_runs_v2_fixedseeds/` |
| `run_convergence_best.sh` | Highest-resolution point of the grid (1280^3, 40 steps) | `results/convergence_runs_v2_fixedseeds/` |
| `run_gr_vs_fofr.sh` | GR vs f(R) at (f,n) = (512,512), (1024,512), (1024,1024), (1280,1280), 30 steps, fixed amplitude with normal and reversed phases, one seed, identical ICs for both models (`DRY_RUN=1` to only write parameter files) | `results/gr_vs_fofr/<run>/` |
| `run_weekend.sh` | Weekend campaigns, one simulation at a time, f = n = 1024, 64 tasks x 2 threads, outputs z = 2, 1, 0.5, 0, all phase-reversed pairs of GR/F4/F5/F6: `ensemble` (10 seeds, 1024 Mpc/h plus alternately 2048 and 512 Mpc/h, 30 steps), then `steps` (20/30/40 steps, first seed) | `results/weekend/<run>/` |
| `run_resolution_scan.sh` | force mesh x box x particles scan, gnuplot figures only (deletes raw output) | `figures/resolution_scan/` |

The 2048^3 reference is run by hand from the LagrangianBias root:
`mpirun -np <N> ../COLASolver/nbody params/parameterfile_2048.lua` (writes to `results/reference_2048/`).

## Plot scripts (`scripts/plot/`)

| Script | Reads | Writes to |
|---|---|---|
| `plot_ensemble_mean.py` | `results/ensemble_100seeds/` | `figures/ensemble_100seeds/ensemble_spectra_1to1.pdf` |
| `plot_ensemble_variance.py` | `results/ensemble_100seeds/` | `figures/ensemble_100seeds/ensemble_spectra_variance.pdf` |
| `plot_convergence_spectra.py` | `results/convergence_runs_v2_fixedseeds/` | `figures/convergence_spectra/` |
| `plot_convergence_spectra_ref2048.py` | same, plus `results/reference_2048/` | `figures/convergence_spectra_ref2048/` |
| `plot_convergence_errors.py` | `results/convergence_runs_v2_fixedseeds/` | `figures/convergence_errors/` |
| `plot_fofr_boost.py` | a campaign folder, default `results/gr_vs_fofr/` (finished outputs only, so it works mid-campaign) | `figures/fofr_boost/boost_<campaign>_z<z>.pdf` |
| `compare_pofk_nmesh.py` | `results/convergence_runs_v1/f1024_n1024_s40_v4` vs `output/` | interactive windows |
| `plot_pofk_quick.gnuplot` | `output/snapshot_TestSim_z0.000/` | interactive window |

Older figure sets are kept under `figures/` with the date they were made.
