# Lagrangian bias spectra

Basis spectra of the hybrid Lagrangian bias expansion from COLA simulations: the matter field and the
weights {delta_L, delta_L^2 - <delta_L^2>, s^2 - <s^2>, nabla^2 delta_L}, evaluated at the particles'
Lagrangian positions q, advected with the particles, and all 15 auto and cross power spectra `pofk_ij.txt`.

## Code

`LagrangianBias.h` (namespace `FML::LAGRANGIANBIAS`) does all the work and is independent of the solver:

| Function | What it does |
|---|---|
| `compute_bias_weights` | Builds the weights at q from the initial field, scaled per mode with D(k, a_out) / D(k, a_ini) |
| `compute_bias_power_spectra` | Deposits the 5 fields at the Eulerian positions and writes `pofk_ij.txt` and `pofk_bias_info.txt` |
| `reset_bias_weights` | Sets the weights to matter only (particles read from file, no linear field) |

`COLASolver` calls it from three places:
- `src/Main.cpp`: the `Particle` carries `bias_weights[]`, `active_bias_index` and `get_mass()`
- `src/Simulation.h`: keeps the initial field (`bias_delta_ini_fourier`) and calls `compute_bias_weights` at each output
- `src/AnalyzeOutput.h`: calls `compute_bias_power_spectra` when `pofk = true`

Build the solver as usual in `../COLASolver` (`make`). The scripts here run `../COLASolver/nbody`.

## Running simulations (`scripts/run/`)

Every campaign is a short Python file whose loops are the parameters it scans; edit the lists to change them.
`simulations.py` does the work: `run("campaign/name", **parameters)` writes `parameterfile.lua` with those parameters
changed and runs `../COLASolver/nbody` in `results/campaign/name/`, which then holds the parameter file, `log.txt` and
the `snapshot_*` output. Finished runs are skipped, so a campaign can be restarted, and `--dry-run` only writes the
parameter files. `results/<campaign>/current_log.txt` always points to the running simulation's log.

| Campaign | Runs |
|---|---|
| `convergence.py` | force mesh x particles x time steps in 1024 Mpc/h, 10 seeds each |
| `ensemble_100.py` | 100 seeds of `parameterfile.lua` |
| `gr_vs_fofr.py` | GR vs f(R) at four resolutions, one seed as a phase-reversed pair |
| `weekend.py` | GR vs F4/F5/F6 at z = 2, 1, 0.5, 0: 10 seeds in 1024 Mpc/h plus 2048/512 Mpc/h, and 20/30/40 steps |
| `resolution_scan.py` | force mesh x box x particles, each in {256, 512, 1024} |

For example `python3 scripts/run/weekend.py --dry-run`. The 2048^3 reference is run by hand from this folder:
`mpirun -np <N> ../COLASolver/nbody params/parameterfile_2048.lua` (writes to `results/reference_2048/`).

## Plotting (`scripts/plot/`)

Every figure shows all 15 spectra on one page, a panel per spectrum, with at most 5 lines. When more is compared,
the PDF gets several pages, each varying one parameter (`spectra.py` does the loading, plotting and page splitting).

| Script | Writes |
|---|---|
| `plot_runs.py OUTPUT GROUP...` | any runs, e.g. `plot_runs.py ensemble/all_spectra.pdf 'results/ensemble_100seeds/*'` (mean and 1 sigma band), or `--ratio` to divide by the first group |
| `plot_convergence.py` | `figures/convergence/all_spectra.pdf`: spectra and ratios to the best setup, varying f, n or s per page |
| `plot_fofr_boost.py [campaign]` | `figures/fofr_boost/all_boost_<campaign>.pdf`: f(R)/GR per redshift, varying model, L, f, n or s per page |

## Folder layout (LagrangianBias/)

| Folder | Contents | In git |
|---|---|---|
| `LagrangianBias.h` | The bias weights and spectra (see above) | yes |
| `parameterfile.lua` | Base parameter file for all campaigns | yes |
| `params/` | Extra parameter files (`parameterfile_2048.lua`: 2048^3 reference run) | yes |
| `scripts/` | Campaigns and plots, see above | yes |
| `results/<campaign>/<run>/` | Parameter file, log and `snapshot_*` output of each run | READMEs only |
| `figures/` | Generated plots | no |
| `output/` | Default output folder of `parameterfile.lua` | no |
| `hostlist.txt` | Nodes for MPI runs | yes |

Older results (`convergence_runs_v1`, `convergence_runs_v2_fixedseeds`, `ensemble_100seeds`) have only the spectra
of each run, and parameter files saved with older results have `COLASolver/` paths; they are a record of the
settings, not meant to be rerun as they are. Older figure sets are kept under `figures/` with the date they were made.
