# Lagrangian bias spectra

Basis spectra of the hybrid Lagrangian bias expansion from COLA simulations: the matter field and the
weights {delta_L, delta_L^2 - <delta_L^2>, s^2 - <s^2>, nabla^2 delta_L}, evaluated at the particles'
Lagrangian positions q, advected with the particles, and all 15 auto and cross power spectra `pofk_ij.txt`.

## Scripts

| Script | What it does |
|---|---|
| `scripts/runs.py` | Runs COLASolver with any parameters of `parameterfile.lua` changed: `run()` for one run, `scan()` for a set |
| `scripts/plot_page.py` | All 15 spectra on one page, at most 5 lines per page (more are split into pages varying one parameter) |
| `scripts/plot_single.py` | One file per spectrum (`P_00.pdf`, ...), with every curve in one figure |
| `scripts/spectra.py` | Loads the spectra of runs and combines them (means, ratios, f(R)/GR boosts) into figures for the two above |
| `scripts/run_*.py`, `scan_*.py` | The campaigns: the runs, then the figures |

Every run script takes `--dry-run` (only write the parameter files and print the commands) and `--plot-only` (run
nothing, only make the figures from the finished runs). Finished runs are skipped, so a script can be restarted or
extended, and `results/<campaign>/current_log.txt` points to the log of the running simulation (`tail -F`).

One run from the command line (lists as Python lists, `true`/`false` or `True`/`False`):

```
python3 scripts/runs.py tests/mesh1024 pofk_nmesh=1024 cosmology_h=0.7
```

Plots of any runs from the command line (a run folder or a quoted glob pattern per line, averaged):

```
python3 scripts/plot_page.py ensemble/all_spectra.pdf 'results/ensemble_100seeds/*'
python3 scripts/plot_single.py ensemble 'results/ensemble_100seeds/*'
python3 scripts/plot_page.py scans/mesh.pdf results/scans/planck/fiducial_GR_seed1001_normal 'results/scans/planck_mesh/*' --ratio
```

## Quick start: scanning parameters

A scan is a small Python file: `scan()` does the runs, `scan_figures()` makes the figures, and `save_pages()` and
`save_single()` write them. `scripts/scan_cosmology.py` changes Omega_m and A_s one at a time by -50%, -10%, +10% and
+50% around the fiducial cosmology of `parameterfile.lua`, with two seeds and both phases:

```python
FRACTIONS = [-0.5, -0.1, 0.1, 0.5]
Omega_b = fiducial("cosmology_Omegab")
Omega_m = fiducial("cosmology_OmegaCDM") + Omega_b

cosmology = scan("cosmology", seeds=[1001, 1002], one_at_a_time={
    "cosmology_OmegaCDM": [round(Omega_m * (1 + fraction) - Omega_b, 6) for fraction in FRACTIONS],
    "cosmology_As": around_fiducial("cosmology_As", FRACTIONS),
})
figures = scan_figures(cosmology)
save_pages(figures, "scans/cosmology.pdf")
save_single(figures, "scans/cosmology")
```

Run it from this folder with `python3 scripts/scan_cosmology.py` (1 fiducial + 4 + 4 cosmologies x 2 seeds x 2 phases
= 36 runs). For each output redshift, `scans/cosmology.pdf` has a page per parameter (mean over seeds and phases, with a
1 sigma band), then the same divided by the fiducial run; `scans/cosmology/P_ij.pdf` has all of them in one figure.

- `one_at_a_time={parameter: [values]}` changes one parameter at a time, the others staying fiducial, plus one
  fiducial run. The fiducial values come from `fixed`, or else from `parameterfile.lua`: `scripts/scan_planck.py`
  centres on Planck 2018 and steps Omega_m, H_0, A_s, w0 and wa by -5, -1, +1 and +5 sigma. Add parameters as lines,
  e.g. `"cosmology_ns": around_fiducial("cosmology_ns", [-0.05, 0.05])`.
- `vary={parameter: [values]}` instead runs every combination, e.g. `vary={"force_nmesh": [512, 1024], "particle_Npart_1D": [512, 1024]}`.
- `fixed={parameter: value}` for all runs, e.g. `fixed={"simulation_boxsize": 2048, "output_redshifts": [1.0, 0.0]}`.
- `seeds` (default: the one in `parameterfile.lua`), `phases=["normal"]` to skip the reversed phases, `ntasks` (default 64).
- Everything else comes from `parameterfile.lua`: 1024 Mpc/h, 1024^3 particles and force mesh, 30 time steps, ICs at
  z = 20, output at z = 0, GR. Change the defaults there, or per scan with `fixed`.
- `gravity=["F5"]` (or `["F4", "F5", "F6"]`) also runs every point in these f(R) models and in GR, all from the same
  initial conditions (`ic_use_gravity_model_GR = true`), and adds the boosts f(R)/GR to the figures. Without it, scans
  run the gravity model of `parameterfile.lua` (GR).
- Long scans: `nohup python3 scripts/scan_cosmology.py > scan_cosmology.log 2>&1 &`, or in tmux.

Runs go to `results/scans/<name>/<values>[_<gravity>]_seed<seed>_<phase>/`, where `<values>` are the parameters that
differ from the fiducial, without the `cosmology_` prefix (e.g. `w0-0.6_F5_seed1001_normal`), or `fiducial`. Run folder
names hold only the varied parameters, so use a new scan name when changing `fixed`. Changing `cosmology_*` parameters
also changes the initial power spectrum, so a new linear P(k) is computed for each cosmology by running CLASS
(`write_input_power` in `scripts/runs.py`; set the executable there or with the environment variable `CLASS`). Each Lua
parameter sets one CLASS parameter
(`LUA_TO_CLASS`: `cosmology_h` -> `h`, `cosmology_Omegab` -> `Omega_b`, ...), and `cosmology_Neffective` and
`cosmology_OmegaMNu` give COLASolver's neutrinos (3 species sharing the mass; `N_ur = Neffective` if massless). It is saved, with its CLASS input file, in
`results/scans/<name>/input/` (for a single `run()`, in `results/<folder>/input_power.txt`); for the fiducial cosmology it matches `../COLASolver/input/example_power_spectrum_cb_z0.000.txt`
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

## Campaigns

| Script | Runs | Figures |
|---|---|---|
| `run_convergence.py` | force mesh x particles x time steps in 1024 Mpc/h, 10 seeds each | `convergence/`: spectra and ratios to the best setup |
| `run_ensemble_100.py` | 100 seeds at 512^3, 10 steps | `ensemble/`: mean and scatter |
| `run_gr_vs_fofr.py` | GR vs f(R) at four resolutions, one seed as a phase-reversed pair | `fofr_boost/all_boost_gr_vs_fofr.pdf`, `fofr_boost/gr_vs_fofr/` |
| `run_weekend.py` | GR vs F4/F5/F6 at z = 2, 1, 0.5, 0: 10 seeds in 1024 Mpc/h plus 2048/512 Mpc/h, and 20/30/40 steps | `fofr_boost/all_boost_weekend.pdf`, `fofr_boost/weekend/` |
| `run_resolution_scan.py` | force mesh x box x particles, each in {256, 512, 1024} | `resolution_scan/` |
| `scan_cosmology.py`, `scan_planck.py` | the scans above | `scans/<name>.pdf`, `scans/<name>/` |
| `scan_planck_mesh.py` | the Planck fiducial with `pofk_nmesh = 1024` | `scans/planck_mesh.pdf` |
| `plot_fixes.py` | (no runs) the -k^2 rule for nabla^2 delta_L, smoothed boosts and normalised differences | `scans/planck_fixes.pdf`, `scans/planck_fixes/` |

The 2048^3 reference is run by hand from this folder:
`mpirun -np <N> ../COLASolver/nbody params/parameterfile_2048.lua` (writes to `results/reference_2048/`).

## Folder layout (LagrangianBias/)

| Folder | Contents | In git |
|---|---|---|
| `LagrangianBias.h` | The bias weights and spectra (see above) | yes |
| `parameterfile.lua` | Base parameter file for all scans and campaigns | yes |
| `params/` | Extra parameter files (`parameterfile_2048.lua`: 2048^3 reference run) | yes |
| `scripts/` | Run and plot scripts, see above | yes |
| `results/<campaign>/<run>/` | Parameter file, log and `snapshot_*` output of each run | READMEs only |
| `figures/` | Generated plots | no |
| `output/` | Default output folder of `parameterfile.lua` | no |
| `hostlist.txt` | Nodes for MPI runs | yes |

Older results (`convergence_runs_v1`, `convergence_runs_v2_fixedseeds`, `ensemble_100seeds`) have only the spectra
of each run, and parameter files saved with older results have `COLASolver/` paths; they are a record of the
settings, not meant to be rerun as they are. Older figure sets are kept under `figures/` with the date they were made.
