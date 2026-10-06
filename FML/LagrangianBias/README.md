# Lagrangian bias spectra

Basis spectra of the hybrid Lagrangian bias expansion from COLA simulations: the matter field and the
weights {delta_L, delta_L^2 - <delta_L^2>, s^2 - <s^2>, nabla^2 delta_L}, evaluated at the particles'
Lagrangian positions q, advected with the particles, and all 15 auto and cross power spectra `pofk_ij.txt`.

## Quick start: scanning any parameter

A scan is a small Python file. Copy `scripts/scan_omega_m.py`, which is all of this:

```python
from scan import scan

scan("omega_m", vary={"cosmology_OmegaCDM": [0.151, 0.251, 0.351]}, seeds=[1001, 1002])
```

and run it from this folder with `python3 scripts/scan_omega_m.py`. It runs every combination of the values in `vary`
(here 3 values x 2 seeds x 2 phases = 12 runs), skips runs that have finished, and writes `figures/scans/omega_m.pdf`:
for each output redshift, pages varying one parameter (mean over seeds and phases, with a 1 sigma band), each followed by
the same divided by the first value.

- `vary`: any parameters of `parameterfile.lua`, with lists of values; several give every combination, e.g.
  `vary={"force_nmesh": [512, 1024], "particle_Npart_1D": [512, 1024]}`
- `fixed`: parameters for all runs, e.g. `fixed={"simulation_boxsize": 2048, "output_redshifts": [1.0, 0.0]}`
- `seeds` (default: the one in `parameterfile.lua`), `phases=(False,)` to skip the reversed phases, `ntasks` (default 64)
- Add values or seeds later and rerun: only the new runs are done.
- `--dry-run` only writes the parameter files; `--plot-only` only remakes the figure (also while runs are going).
- Long scans: `nohup python3 scripts/scan_omega_m.py > scan_omega_m.log 2>&1 &`, and follow the running simulation
  with `tail -F results/scans/omega_m/current_log.txt`.

Runs go to `results/scans/<name>/<values>_seed<seed>_<phase>/`. Changing `cosmology_*` parameters also changes the
initial power spectrum, so a new linear P(k) is computed for each cosmology by running CLASS (`scripts/input_power.py`;
set the executable there or with the environment variable `CLASS`). Each Lua parameter sets one CLASS parameter
(`LUA_TO_CLASS`: `cosmology_h` -> `h`, `cosmology_Omegab` -> `Omega_b`, ...), and `cosmology_Neffective` and
`cosmology_OmegaMNu` give COLASolver's neutrinos (3 species sharing the mass; `N_ur = Neffective` if massless). It is saved, with its CLASS input file, in
`results/scans/<name>/input/`; for the fiducial cosmology it matches `../COLASolver/input/example_power_spectrum_cb_z0.000.txt`
to 0.2%. `cosmology_model = "LCDM"` and `"w0waCDM"` are set up (for w0waCDM, `cosmology_w0` -> `w0_fld`,
`cosmology_wa` -> `wa_fld`, with `Omega_Lambda = 0`), e.g.
`scan("w0", vary={"cosmology_w0": [-1.0, -0.9]}, fixed={"cosmology_model": "w0waCDM", "cosmology_wa": 0.0})`.
Parameters you don't vary keep their `parameterfile.lua` values.
- Redshift: CLASS computes P(k) at `ic_input_redshift` (0), and FML scales it back to `ic_initial_redshift` with its own
  growth factor, so the runs reproduce the CLASS P(k) at z = 0. To use CLASS at the initial redshift instead, add
  `fixed={"ic_input_redshift": 20.0}` (equal to `ic_initial_redshift`).
- k range: the table goes from about 1e-5 h/Mpc to at least 100 h/Mpc, and further when a run's IC grid needs it
  (sqrt(3) times the Nyquist frequency pi particle_Npart_1D / simulation_boxsize); FML holds P(k) constant beyond it.
- Amplitude: set by `cosmology_As`. To fix sigma_8 instead, add `fixed={"ic_sigma8_normalization": True, "ic_sigma8": 0.81}`.

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

## Campaigns (`scripts/run_*.py`)

Every campaign is a short Python file whose loops are the parameters it scans; edit the lists to change them.
`simulations.py` does the work: `run("campaign/name", **parameters)` writes `parameterfile.lua` with those parameters
changed and runs `../COLASolver/nbody` in `results/campaign/name/`, which then holds the parameter file, `log.txt` and
the `snapshot_*` output. Finished runs are skipped, so a campaign can be restarted, and `--dry-run` only writes the
parameter files. `results/<campaign>/current_log.txt` always points to the running simulation's log.

| Campaign | Runs |
|---|---|
| `run_convergence.py` | force mesh x particles x time steps in 1024 Mpc/h, 10 seeds each |
| `run_ensemble_100.py` | 100 seeds of `parameterfile.lua` |
| `run_gr_vs_fofr.py` | GR vs f(R) at four resolutions, one seed as a phase-reversed pair |
| `run_weekend.py` | GR vs F4/F5/F6 at z = 2, 1, 0.5, 0: 10 seeds in 1024 Mpc/h plus 2048/512 Mpc/h, and 20/30/40 steps |
| `run_resolution_scan.py` | force mesh x box x particles, each in {256, 512, 1024} |

For example `python3 scripts/run_weekend.py --dry-run`. The 2048^3 reference is run by hand from this folder:
`mpirun -np <N> ../COLASolver/nbody params/parameterfile_2048.lua` (writes to `results/reference_2048/`).

## Plotting (`scripts/plot_*.py`)

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
