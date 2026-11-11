# Hybrid bias spectra

Basis spectra of the hybrid Lagrangian bias expansion from COLA simulations: the matter field and the
weights {delta_L, delta_L^2 - <delta_L^2>, s^2 - <s^2>, nabla^2 delta_L}, evaluated at the particles'
Lagrangian positions q, advected with the particles, and all 15 auto and cross power spectra, written by
`FML/SimulationCode` to `snapshot_*/hybridbias_z<z>.txt`. (Until November 2026 this project was `LagrangianBias/` and ran
its own `LagrangianBias.h` in COLASolver; see "Runs made with COLASolver" below.)

Quick test (3 s, 128^3 particles in 256 Mpc/h, from this folder, after building `../SimulationCode`):

```
python3 scripts/runs.py tests/quick simulation_boxsize=256.0 particle_Npart_1D=128 force_nmesh=128 pofk_nmesh=128 "timestep_nsteps=[10]" ntasks=8
```

## Scripts

| Script | What it does |
|---|---|
| `scripts/runs.py` | Runs SimulationCode with any parameters of `parameterfile.lua` changed: `run()` for one run, `scan()` for a set |
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
- Everything else comes from `PROJECT` in `scripts/runs.py` (1024 Mpc/h, 1024^3 particles and force mesh, 30 time
  steps, the bias spectra) and `parameterfile.lua` (ICs at z = 20, output at z = 0, GR, ...). Change the project
  defaults in `PROJECT`, or per scan with `fixed`.
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
`cosmology_OmegaMNu` give SimulationCode's neutrinos (3 species sharing the mass; `N_ur = Neffective` if massless;
massive neutrinos need `ic_type_of_input = "class"` or `"internal"`, not this P(k)). It is saved, with its CLASS input file, in
`results/scans/<name>/input/` (for a single `run()`, in `results/<folder>/input_power.txt`); for the fiducial cosmology it matches `../SimulationCode/input/example_power_spectrum_cb_z0.000.txt`
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

The spectra are computed by `FML/SimulationCode` (upstream `FMLDevelop`), with the `hybridbias_*` parameters
(`../ComputePowerSpectra/HybridBiasSpectra.h`). The weights are built at each output from the initial field at the
particles' lattice positions, scaled per mode with D(k, a_out) / D(k, a_ini), and deposited at the Eulerian positions.
`snapshot_*/hybridbias_z<z>.txt` has the columns k, number of modes, the 15 P_ij (in the order of `PAIRS` in
`spectra.py`) and the linear P_cb(k); the Poisson shot noise of each is in the header. Optional, in the same folder:

| Parameter | File | What it is |
|---|---|---|
| `hybridbias_zeldovich` | `hybridbias_zeldovich_z<z>.txt` | Zel'dovich prediction of the 15 spectra |
| `hybridbias_oneloop` | `hybridbias_oneloop_z<z>.txt` | One-loop LPT (CLEFT, no counterterms); `hybridbias_oneloop_kernels = "model"` for MG kernels |
| `hybridbias_zeldovich_field` | `hybridbias_zeldovichfield_z<z>.txt` | The Zel'dovich field of the same initial conditions |
| both of the first and third | `hybridbias_cv_z<z>.txt` | Control variate P - beta(k) (P_zeldovichfield - P_zeldovich) (Kokron et al. 2022), with beta(k) as an extra column |

The predictions and the control variate need `hybridbias_kcut_over_nyquist > 0` (a spherical cut of the linear field;
1.0 is the SimulationCode default, `parameterfile.lua` has 0.0 as COLASolver). The predictions take about a minute per
output and are the same for every seed of a cosmology, so compute them in one run per cosmology.
`hybridbias_zeldovich_field_only = true` adds the Zel'dovich field to an existing run without time stepping.

Build the solver in `../SimulationCode` with `make -j16` (the library paths are in `../Makefile.local`). The scripts here
run `../SimulationCode/nbody`; `./nbody --parameters` lists every parameter.

## Runs made with COLASolver (before November 2026)

`parameterfile.lua` is SimulationCode's own parameter file, with the old COLASolver values of this project and the
reasons for the changes as comments (the old file is `params/parameterfile_colasolver.lua`). The project's settings for
every run (box, particles, force grid, steps, the bias spectra) are in `PROJECT` at the top of `scripts/runs.py`, with
the values of the old runs, so new runs reproduce them; a run's own parameters override `PROJECT`. The lattice offset
and gradient kernel are the SimulationCode defaults, which make no difference when `force_nmesh = particle_Npart_1D`.
The scripts read the old runs (15 `pofk_ij.txt` files) as before and count them as finished, so finished scans are not
rerun. Same seed, 128^3 particles, 10 steps, the two codes side by side:

| Runs | Matter P(k) | The 15 spectra |
|---|---|---|
| `force_nmesh = particle_Npart_1D` (most runs) | 2e-4 | 1e-3 (nabla^2 delta_L; the rest ~1e-4) |
| `force_nmesh = 2 particle_Npart_1D` | 2e-3 | 1.5e-2 (nabla^2 delta_L) |

(differences relative to sqrt(P_ii P_jj), k < 0.5 h/Mpc). With a force grid finer than the particle grid the old setup
(particles on the grid nodes, continuous gradient) is itself inaccurate (about +10% at k = 0.2 h/Mpc for 2x, see
`../SimulationCode/README.md`), and SimulationCode's fix of the Nyquist mode of the gradient changes it; such runs can't
be reproduced exactly. Rerun them under a new scan name rather than mixing old and new runs (add
`ic_lattice_offset = "none"`, `force_gradient_kernel = "continuous"` to a run to get as close as possible to the old ones).

### Force grid and particle number at the production spacing

Scaled-down runs with the spacings of the production runs (256 Mpc/h with 128^3 or 256^3 particles is 512^3 or 1024^3
in 1024 Mpc/h; force grid 256^3 = 1024^3 there), one seed, 30 steps, SimulationCode defaults (lattice offset and
discrete gradient on finer grids), against 1024^3 particles on a 3072^3 force grid (P_ij ratio - 1, %):

| Setup (as 1024 Mpc/h) | k (h/Mpc) | P_11 | P_delta delta | P_delta2 delta2 | P_s2 s2 | P_1 s2 | P_1 nabla2 | P_nabla2 nabla2 |
|---|---|---|---|---|---|---|---|---|
| 512^3 on 1024^3 | 0.03-0.1 | 0.0 | 0.0 | -1.8 | -2.7 | -1.9 | +9.5 | -13 |
| | 0.2-0.4 | -0.7 | -0.7 | -2.6 | -5.6 | +3.8 | -14 | -30 |
| | 0.4-0.8 | -3.2 | -2.6 | -2.6 | -5.5 | +3.0 | -19 | -34 |
| 1024^3 on 1024^3 (old setup) | 0.03-0.1 | -0.1 | +0.4 | +0.2 | +0.1 | +7.1 | +85 | +166 |
| | 0.2-0.4 | -1.1 | -0.3 | +0.2 | +1.2 | +2.6 | -15 | -29 |
| | 0.4-0.8 | -3.8 | -2.5 | -1.1 | +1.9 | +1.3 | -21 | -30 |
| 1024^3 on 2048^3 | 0.03-0.1 | 0.0 | 0.0 | 0.0 | 0.0 | +0.1 | +3.7 | +4.8 |
| | 0.2-0.4 | 0.0 | -0.1 | -0.3 | -0.3 | +0.1 | -3.2 | -5.9 |
| | 0.4-0.8 | -0.2 | -0.3 | -0.6 | -0.4 | +0.6 | -4.9 | -6.0 |

1024^3 on 2048^3 is within 1% in all but the nabla^2 delta_L spectra (5%); the time stepping took 2.5 times as long as
on 1024^3 (the 3072^3 reference 7 times). 512^3 on 1024^3 (a third of the time) and 1024^3 on 1024^3 are about equally
good for the matter and delta_L spectra (to 1% below k = 0.4 h/Mpc); the first is worse in delta^2 and s^2 (2-6%) and
has 8 times the Poisson shot noise (V/N = 8 (Mpc/h)^3) and half the lattice Nyquist frequency, the second has large
errors in the nabla^2 delta_L spectra at low k. (A coarser test, 128^3 in 512 Mpc/h on 128^3, gave +7% in P_delta delta
at low k as the SimulationCode README finds; at the production spacing it is 0.4%.) Memory: see "What fits in 500 GB"
below.

### Convergence against a finer reference (particle number and force grid)

Convergence: the same scaled-down setup (256 Mpc/h, 30 steps, seed 1001, R = 1 Mpc/h Gaussian smoothing so that all
particle numbers have the same operators) against 384^3 particles on a 768^3 grid (1536^3 on 3072^3 in 1024 Mpc/h);
P_ij ratio - 1 in %, for k = 0.3-0.6 / 0.6-1.0 h/Mpc:

| Setup (as 1024 Mpc/h) | P_11 | P_delta delta | P_s2 s2 | P_1 nabla2 |
|---|---|---|---|---|
| 512^3 on 1024^3 | -1.8 / -6.0 | -0.9 / -4.1 | +1.9 / +4.0 | -4.4 / -9.8 |
| 1024^3 on 1024^3 (old setup) | -2.3 / -6.7 | -0.5 / -3.8 | +2.1 / -0.3 | -5.6 / -10 |
| 768^3 on 1536^3 | -0.5 / -2.3 | -0.2 / -1.5 | +0.3 / -0.2 | -2.1 / -4.3 |
| 1024^3 on 1536^3 | -0.7 / -2.6 | 0.0 / -1.3 | +0.5 / -0.3 | -2.5 / -4.8 |
| 1024^3 on 2048^3 | -0.1 / -0.7 | +0.2 / -0.2 | 0.0 / -0.1 | -0.7 / -1.7 |

Below k = 0.3 h/Mpc all are within 1% except the nabla^2 delta_L spectra at k < 0.1 (up to +55% for 1024^3 on 1024^3).
The force grid matters more than the particle number (768^3 and 1024^3 on 1536^3 are about the same).

Scale-dependent mu(a,k) (`gravity_model_mu_parametrization = "BZ"` or a table): for BZ {4/3, 5 Mpc/h, 4} (mu = 1.07 at
k = 0.1 and 1.32 at k = 1 h/Mpc today, 1.0 at z > 1) with 30 steps, standard and scale-dependent COLA
(`simulation_use_scaledependent_cola = true`) agree to 1e-4 in all 15 spectra at z = 0, and scale-dependent COLA takes
twice as long and more memory (table below); models whose scale dependence starts earlier, or fewer steps, may need
scale-dependent COLA (check once).

### What fits in 500 GB (measured on euclid22)

1024^3 particles in 1024 Mpc/h, 64 MPI tasks on euclid22 (502 GB, AMD EPYC 7543, 128 threads), 2 time steps (the
memory peak is in the first step and at the output), under a watchdog; node memory in use at the peak (490 GB were free
before), none of them used swap. Time per step from the 2 steps (rough), the 30-step estimate without the output:

| Model, force grid | Precision | Memory | Per step | 30 steps |
|---|---|---|---|---|
| GR / mu(a), 1536^3 | double | 308 GB | 28 s | ~15 min |
| mu(a) (mu0), 2048^3 | double | **456 GB** | 55 s | ~27 min |
| mu(a), 2048^3, with the Zel'dovich field and prediction | double | 453 GB | | +15 min per output (see below) |
| f(R) F5, 1536^3 | double | 329 GB | 58 s | ~29 min |
| nDGP H0 rc = 1, 1536^3 | double | 315 GB | 60 s | ~30 min |
| mu(a,k) BZ, 1536^3, standard COLA | double | 333 GB | 29 s | ~15 min |
| mu(a,k) BZ, 1536^3, scale-dependent COLA | double | 388 GB | 62 s | ~31 min |
| mu(a), 2048^3 | float | 320 GB | 37 s | ~18 min |
| nDGP, 2048^3 | float | 353 GB | 103 s | ~51 min |
| mu(a,k), 2048^3, standard COLA | float | 354 GB | 39 s | ~20 min |
| f(R) F5, 2048^3 | float | 353 GB | 93 s | ~46 min |

(For comparison: the COLASolver runs with 1280^3 particles on 1280^3 peaked at ~540 GB and swapped.)

- Every model fits with 1024^3 particles on a 2048^3 force grid in single precision, and GR / mu(a) also in double
  precision. Single precision changes the 15 spectra by < 6e-5 (k < 0.1 h/Mpc) and ~1e-6 above (same seed, 2048^3) and
  is faster. The single-precision build is `../SimulationCode/nbody_float` (CMake, FFTW in float from
  /mn/stornext/u3/hansw/local, without FFTW threads); run it with the scripts as `NBODY=../SimulationCode/nbody_float
  python3 scripts/...`. Rebuild it with
  `cmake -S ../SimulationCode -B <dir> -DCMAKE_BUILD_TYPE=Release -DFML_FFTW_PRECISION=float -DFML_USE_FFTW_THREADS=OFF
  -DFML_BUILD_TESTS=OFF -DFFTW_ROOT=/mn/stornext/u3/hansw/local -DCMAKE_PREFIX_PATH="/mn/stornext/u3/hansw/local;$HOME/local"
  -DCMAKE_CXX_COMPILER=mpicxx && cmake --build <dir> -j16`.
- Both binaries were built on rastaban (Intel Broadwell, -march=native) and run on euclid22; building them on euclid22
  (-march=native for its AMD CPU) may be faster.
- The control variate: the Zel'dovich field takes ~80 s per output and no extra memory, but the analytic Zel'dovich
  prediction ~14 min per output. The prediction is the same for every seed of a cosmology, so turn on
  `hybridbias_zeldovich_field` in every run and `hybridbias_zeldovich` in one run per cosmology, and form
  P_cv = P - beta(k) (P_zeldovichfield - P_zeldovich) afterwards.
- Many MPI tasks on a small grid fail (e.g. 64 tasks on 128^3: fewer slices per task than the assignment kernel needs).

### Production validation (euclid22)

1024^3 particles in 1024 Mpc/h, 30 steps over the outputs z = 2, 1, 0.5, 0, 64 tasks, single precision, bias spectra on
512^3 with the cut at 1.0, Zel'dovich field and prediction at every output, seed 1230 as the 2048^3 reference run
(`results/reference_2048`, COLASolver, 2048^3 particles on 2048^3); P_11 / reference - 1 at z = 0, for k = 0.3-0.6 /
0.6-1.0 / 1.0-1.5 h/Mpc:

| Force grid | P_11 | Run time (of which the 4 Zel'dovich predictions) | Memory |
|---|---|---|---|
| 1024^3 (the old setup) | -2.4 / -6.3 / -13.5% | 40 min (36 min) | ~300 GB |
| 2048^3 | -0.06 / -0.24 / -0.63% | 51 min (36 min) | ~300 GB |

(The other spectra differ from the reference also by definition: without smoothing the operators depend on the
particle lattice, 2048^3 there, and the reference had no cut.) The Zel'dovich prediction takes ~8 min per output with
hybridbias_nmesh 512, the field ~1 min, the spectra ~1 min. Native builds on euclid22 (`../build_euclid22.sh`, AMD
-march=native) run exactly as fast as the builds from rastaban (memory and FFT bound). Other tests on euclid22: the
corners mu0 = 1 with (w0, wa) = (-0.5, -2), mu0 = -0.5 with (-1.3, 0.5), and nDGP H0 rc = 0.5 run without problems
(512^3 on 1024^3: 100-280 s for 30 steps); time stepping: COLA with Quinn steps linear in a (the default) is the best
per step (growth-factor spacing is worse with COLA; BullFrog without COLA needs as many steps); nDGP: the boost
converges with resolution as for mu0 (512^3: <= 0.9%), and the approximate screening overestimates it by ~1% in P_11
against the exact multigrid solution (~7% of the MG effect for H0 rc = 1), at about half the cost.

### Smoothing of the operators (SimulationCode, November 2026)

Convergence with the particle number at fixed smoothing R (Gaussian, `hybridbias_smoothing_scale`; cut 1.0): 1024^3
against 1536^3 particles (force grids 2x; 1/64-scale box, same seed, 30 steps), in % for k = 0.3-0.6 / 0.6-1.0 h/Mpc:

| R (Mpc/h) | P_11 | P_delta2 delta2 | P_s2 s2 | P_1 nabla2 | P_nabla2 nabla2 |
|---|---|---|---|---|---|
| 0 | -0.1 / -0.7 | -1.0 / -1.4 | -1.7 / -2.4 | -5.1 / -10 | -15 / -21 (-11 at k < 0.3) |
| 0.5 | -0.1 / -0.7 | +0.1 / -0.3 | +0.2 / +0.6 | -2.3 / -4.6 | -1.7 / -1.2 |
| 1 | -0.1 / -0.7 | +0.2 / 0.0 | 0.0 / -0.1 | -0.7 / -1.7 | -0.2 / -0.3 |
| 2 | -0.1 / -0.7 | +0.4 / +0.4 | -0.2 / -0.7 | -0.1 / -0.8 | 0.0 / -0.4 |

Without smoothing the operators depend on the modes the particle lattice holds, so the nabla^2 delta_L spectra don't
converge even on large scales; with R = 1 Mpc/h all 15 converge to <= 0.3% below k = 0.6 and <= 1.7% to k = 1, and the
rest is the matter spectrum itself (resolution of the dynamics). With 512^3 particles the gain is larger still
(P_nabla2 nabla2 -57% -> +0.4% at k = 0.6-1). The smoothing is part of the operator definition (the bias parameters
refer to the smoothed field), so the theory must use the same R.

### Control variates: is the mean of the Zel'dovich field the prediction?

64 Zel'dovich fields (32 seeds, both phases, fixed amplitudes; 1/64-scale box, hybridbias_nmesh 256) against the
prediction, as the bias the control variate would put into the result, beta(k) (mean - prediction) / prediction:
with the cut at 1.0 and no smoothing P_nabla2 nabla2 is biased (+1.9% at k = 0.3-0.6, 8.6 sigma; +0.9% at 0.6-1) and
P_delta2 delta2 slightly (+0.3% at 0.6-1); with the cut at 0.5, or with R = 1 Mpc/h (cut 1.0 or 0.5), all 15 spectra
agree within ~2 sigma (<= 0.07% above k = 0.3). Gaussian initial conditions behave as fixed amplitudes. So use R = 1
Mpc/h, or the cut at 0.5, with control variates.

### Theory predictions per model

The Zel'dovich and one-loop (CLEFT, no counterterms) predictions are available for every model: with
`hybridbias_oneloop_kernels = "model"` the one-loop uses the model's own LPT kernels (the log names them: GR, mu(k,a),
nDGP with screening, f(R) with screening). Time for one output (512 Mpc/h, 256^3, bias grid 256, 32 tasks on rastaban):
GR ~7 min, nDGP ~10 min, mu ~23 min (the general mu(k,a) kernels are solved numerically); the EdS kernels with the
model's linear P(k) are much faster. The predicted MG/GR boost of P_11 agrees with the simulation to ~1 percentage point at
k = 0.05 h/Mpc (mu0 = 0.3: 7.4 vs 7.9%; nDGP H0 rc = 1: 13.2 vs 14.0%) and falls well below it at k > 0.2.

### Cheap runs for many cosmologies (Latin hypercube)

Scaled tests (256 Mpc/h with 256^3 particles on 512^3 = 1024^3 on 2048^3 in 1024 Mpc/h; same seed; mu0 = 0.3 against
GR; R = 1 Mpc/h), against 60 time steps at full resolution; error in % for k = 0.3-0.6 / 0.6-1.0 h/Mpc in P_11 (the
delta and delta^2 spectra are similar, s^2 and nabla^2 delta_L smaller):

| Run | Absolute P_ij | Boost P_mu / P_GR (the boost itself is +15 / +20%) |
|---|---|---|
| 30 steps | -0.4 / -0.6 | -0.4 / -1.0 |
| 20 steps | -1.5 / -2.6 | -0.9 / -2.1 |
| 15 steps | -2.6 / -5.1 | -1.3 / -3.1 |
| 10 steps | -5.2 / -11 | -2.2 / -5.3 |
| 512^3 on 1024^3, 30 steps | -2.1 / -5.9 | -0.3 / -0.9 (delta^2: -1.0 / -2.2) |
| 512^3 on 1024^3, 15 steps | -3.0 / -7.3 | -0.8 / -1.9 |

- Fewer time steps do not cancel in the boost: keep 30 (or test growth-factor spacing / BullFrog, see
  ../SimulationCode/README.md).
- The resolution largely cancels in the boost (5-10 times smaller error than in the spectra), so a boost can be
  computed from 512^3 on 1024^3 runs (an eighth of the work and ~60 GB) for k < 1 h/Mpc. Tested for mu0 only;
  screened models (nDGP, f(R)) depend more on resolution, check once.
- MG and GR runs from the same initial conditions share most of their sample variance, so the boost scatters much
  less than the spectra, but the Zel'dovich control variate still helps it on large scales. 8 seeds of 128^3 on 256^3
  in 256 Mpc/h (mu0 = 0.3 vs GR; scatter of one run, std in %, without -> with the control variate, default beta):

  | k (h/Mpc) | Boost P_11 | Boost P_delta delta | Boost P_s2 s2 | Absolute GR P_delta delta |
  |---|---|---|---|---|
  | 0.02-0.1 | 0.17 -> 0.08 | 0.19 -> 0.06 | 0.27 -> 0.11 | 5.3 -> 0.7 |
  | 0.1-0.3 | 0.31 -> 0.24 | 0.40 -> 0.30 | 0.35 -> 0.22 | 5.5 -> 2.4 |
  | > 0.3 | unchanged | unchanged | unchanged | 4.0 -> 3.3, then unchanged |

  (a 4-9 times smaller variance of the boost at k < 0.1, 1.3-2.6 times at 0.1-0.3; the means agree to < 0.1%. In a
  1024 Mpc/h box every scatter is ~8 times smaller.) The Zel'dovich field costs little; the analytic prediction is
  needed once per cosmology and model.
- With many MPI tasks on few particles, clustering can overfill a task's particle buffer ("Receive count exceeds
  allocated particle capacity"): use particle_allocation_factor = 1.5 for such runs.
- The analytic Zel'dovich prediction (needed per cosmology for the control variate) took 52 s / 76 s with
  hybridbias_nmesh 128 / 256 in the scaled box and ~14 min with 1024 at full size, so its cost grows with the number
  of k bins; use the smallest hybridbias_nmesh that reaches the k you need.

mu-Sigma (Planck, DESI: mu = 1 + mu0 Omega_DE(a)/Omega_DE0, Sigma likewise): only mu changes the dynamics, so the
simulation is `gravity_model = "mu"` with `gravity_model_mu_parametrization = "DE"` and `gravity_model_mu_mu0`, with GR's
memory (1024^3 on 2048^3 fits). Sigma (the lensing potential Phi + Psi) changes no particle or P_ij; it enters only
lensing observables built from the spectra, e.g. P_g kappa ~ Sigma(a) P_1a in the Limber integral.

Other changes that move numbers:
- f(R) on a w0waCDM background: the background curvature R(a) now comes from that background (COLASolver used LCDM's),
  so the F5 boosts at the w0, wa points of `scan_planck.py` change. GR runs and LCDM backgrounds are unaffected.
- `hybridbias_kcut_over_nyquist > 0` (needed for the predictions and the control variate) changes the spectra, mainly
  the nabla^2 delta_L ones, compared with the old runs, which had no cut.
- `gravity_model = "Geff"` is now `"mu"` (`gravity_model_mu_parametrization = "DE"`, `"BZ"` or `"table"`).

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
`mpirun -np <N> ../SimulationCode/nbody params/parameterfile_2048.lua` (writes to `results/reference_2048/`).

## Folder layout (HybridBias/)

| Folder | Contents | In git |
|---|---|---|
| `parameterfile.lua` | SimulationCode's parameter file, with the old values and the changes as comments; every run starts from it, with `PROJECT` from `scripts/runs.py` on top | yes |
| `params/` | `parameterfile_2048.lua`: 2048^3 reference run; `*_colasolver.lua`: the old COLASolver files | yes |
| `scripts/` | Run and plot scripts, see above | yes |
| `results/<campaign>/<run>/` | Parameter file, log and `snapshot_*` output of each run | READMEs only |
| `figures/` | Generated plots | no |
| `output/` | Default output folder of `parameterfile.lua` | no |
| `hostlist.txt` | Nodes for MPI runs | yes |

Older results (`convergence_runs_v1`, `convergence_runs_v2_fixedseeds`, `ensemble_100seeds`) have only the spectra
of each run, and parameter files saved with older results have `COLASolver/` or `LagrangianBias/` paths; they are a
record of the settings, not meant to be rerun as they are. Older figure sets are kept under `figures/` with the date they were made.
