# The Python scripts, line by line

A walkthrough of `scripts/`, in dependency order: the run engine, the data layer, the two plotters, then the run
scripts that sit on top. Line numbers refer to the versions of the files when this was written; if the files change,
the function names still lead to the right place.

---

## 1. [runs.py](scripts/runs.py): runs the simulations

### Header and module constants (lines 1–32)

- **1–14, docstring.** What the file does, plus the two ways to use it: the command line, or `run()`/`scan()` from a
  script.
- **15–25, imports.**
  - `ast` parses command-line values.
  - `itertools` builds parameter combinations.
  - `re` edits the Lua file.
  - `subprocess` starts `mpirun` and CLASS.
  - `tempfile` gives CLASS a short folder.
  - `dataclass` makes the small `Run`/`Scan` records.
  - `numpy` reads CLASS output.
- **28, `ROOT`.** `__file__` is `.../LagrangianBias/scripts/runs.py`. One `dirname` gives `scripts/`, the second gives
  `LagrangianBias/`. All paths are built from this, so the scripts work from any directory.
- **29, `NBODY`.** The solver executable, `../COLASolver/nbody` next to `LagrangianBias`.
- **30, `BASE_PARAMETER_FILE`.** The Lua file every run is copied from.
- **31, `DRY_RUN`, `PLOT_ONLY`.** True if those flags are anywhere on the command line. They're read once at import,
  so every function sees them without passing them around.
- **32, `counts`.** A tally of what happened, printed by `summary()`.

### Reading the base file (35–51)

- **35–38, `base_parameter(name)`.** Finds the parameter's line with a regular expression:
  - `^\s*{name}\s*=\s*` matches the start of the line, optional indentation, the name and `=`.
  - `(.*?)` captures the value lazily, i.e. as short as possible.
  - `\s*(?:--.*)?$` allows trailing spaces and a Lua comment `-- ...`, which are not captured.
  - `re.M` makes `^` and `$` match at every line, not just the start and end of the file.
  - `.group(1)` is the captured value, and `.strip('"')` removes quotes, so `"GR"` becomes `GR`.
  - It always returns a string, e.g. `"1024"` or `"2.215e-9"`.
- **41–46, `fiducial(name)`.** The same value, converted to a float if possible, otherwise left as text (`"GR"`).
  Used for the fiducial values in scans.
- **49–51, `around_fiducial(name, fractions)`.** For each fraction f it returns fiducial × (1 + f), formatted with
  `:.6g` (6 significant digits) and turned back into a float. That keeps folder names short and avoids things like
  `0.30000000000000004`.

### Writing Lua (54–62)

**`lua(value)`** converts a Python value to Lua syntax for the parameter file.
- `bool` must be checked first, because `True` is also an `int` in Python. It becomes `true`/`false`.
- `str` gets quotes.
- `list`/`tuple` become `{a, b}`; this calls itself on each element, so nested values work.
- Everything else (numbers) uses `str()`.

### CLASS input power spectrum (65–132)

- **66, `CLASS`.** The executable path; the environment variable `CLASS` wins if set.
- **69–79, `LUA_TO_CLASS`.** The one-to-one dictionary from Lua names to CLASS names. `ic_input_redshift → z_pk` makes
  CLASS output P(k) at the redshift FML expects as input.
- **81, `LUA_TO_CLASS_W0WA`.** The same mapping plus `w0`/`wa`. `{**A, ...}` copies A and adds keys.
- **84–90, `neutrinos()`.** Converts FML's neutrino description to CLASS's:
  - No mass: all of `N_eff` goes into massless species (`N_ur`).
  - Massive: one CLASS species with degeneracy 3 (`deg_ncdm = 3`, three equal masses) and FML's temperature ratio
    `T_ncdm`.
- **93–132, `write_input_power(parameters, path)`.**
  - **97–98, `value(name)`.** A small inner function that takes the value from `parameters` if given, otherwise from
    the base file. This is the "override, else base file" rule in one line.
  - **101, `kmax`.** √3 × π N / L, the corner of the IC grid. FML's spline needs P(k) up to there.
  - **103–109, model check.** Picks the dictionary for the cosmology model. Any `cosmology_*` parameter that CLASS
    wouldn't receive counts as `unknown`, and the script stops (`SystemExit`) rather than silently ignoring it.
  - **111, CLASS parameters.** Builds `{CLASS name: value}` for every mapped parameter.
  - **112–113, w0waCDM.** Sets `Omega_Lambda = 0` so the fluid replaces Λ.
  - **114, neutrinos.** Adds the neutrino parameters.
  - **115–116, ini text.** The CLASS input as text:
    - `output = mPk` asks for the matter power spectrum.
    - `P_k_max_h/Mpc` is at least 100, and 10% beyond `kmax` if larger.
    - `k_per_decade_for_pk = 32` sets the sampling density.
    - Then one line per parameter.
  - **117–119, cache.** If the output and an identical saved `_class.ini` already exist, CLASS has already computed
    exactly this, so it returns. Comparing the text is simpler and safer than comparing parameter values.
  - **121–128, running CLASS** inside a temporary folder, deleted automatically after the `with` block:
    - **121–123:** writes `input.ini` with `root = <tmp>/`. CLASS refuses long output paths, which is why it runs in a
      short temporary folder.
    - **124:** runs CLASS and captures its output.
    - **125–126:** on failure, stops and shows the last 2000 characters of CLASS's output.
    - **127–128:** reads `00_pk_cb.dat` (CDM + baryons) if neutrinos are massive, else `00_pk.dat`. With massless
      neutrinos those are the same thing, and CLASS doesn't write the `_cb` file.
  - **129–132, saving.** Writes two columns plus a header, saves the ini text next to it (for the record and the
    cache), and prints.

### Running (135–188)

- **135–136, `is_finished(folder)`.** A run is finished if `snapshot_TestSim_z0.000/pofk_44.txt` exists. That is the
  last spectrum written at the last output, z = 0.
- **139–180, `run(name, ntasks=64, **parameters)`.** `**parameters` collects every `key=value` you pass into a dict.
  - **142:** the run folder is `results/<name>`.
  - **143–144:** with `--plot-only`, returns the folder at once, so the run scripts still know every folder.
  - **145–147:** finished runs are counted as skipped and returned.
  - **148:** creates the folder.
  - **149–152, input P(k).** If you changed a cosmology parameter but gave no input file (a direct `run()`), it
    computes one into the run folder and points FML to it. `scan()` always passes `ic_input_filename`, so it never
    gets here.
  - **154–156, the parameter file.** Reads the base file and, for every parameter plus `output_folder`, replaces its
    line:
    - `^(\s*){key}\s*=.*$` matches the whole line and captures the indentation.
    - `\g<1>` puts the indentation back.
    - `{key} = {lua(value)}` writes the new value.
    - A parameter that isn't in the file is never matched, so it is silently ignored.
  - **157–158:** saves it as `results/<name>/parameterfile.lua`.
  - **160:** the command, `mpirun -np 64 nbody <file>`.
  - **161–164:** counts the run; with `--dry-run`, prints the command and stops here.
  - **165–169, log link.** Prints the start and replaces the symlink `results/<campaign>/current_log.txt` with one
    pointing at this run's log. `lexists` also catches a broken old link.
  - **171, environment.**
    - Threads per task = hardware threads ÷ tasks: 128/64 = 2 on euclid22.
    - `I_MPI_PIN_DOMAIN=omp` pins each task to its own cores.
    - `**os.environ` comes last, so your own settings win.
  - **172–174, running.** Runs the solver from `ROOT`, because the base file has relative input paths, and writes all
    its output to `log.txt`.
  - **175–180, checking.** Checks the result with `is_finished` rather than the exit code, since MPI exit codes are
    unreliable. Prints the time, or an error, and returns the folder.
- **183–188, `summary()`.** Prints the tally, except with `--plot-only`, when nothing was attempted. With `--dry-run`
  it ends the script (`sys.exit()`), so no figures are made from runs that don't exist.

### Helpers for campaigns (191–207)

- **191–195, `gravity_parameters(model)`.** `"GR"` gives `{gravity_model: "GR"}`. `"F5"` gives `f(R)` with
  `fofr0 = 10^-5`, where `int(model[1:])` is the number after the F.
- **198–207, `steps_per_interval(nsteps, redshifts)`.** COLASolver needs one step count per interval between outputs.
  - **202:** the scale factors: initial, then each output.
  - **203:** the exact share of each interval, nsteps × Δa / (1 − a_ini), in proportion to the change in a.
  - **204:** rounds each down.
  - **205–206:** the missing steps go one each to the intervals with the largest fractional remainders, so the total
    is exactly nsteps. For 30 steps at z = 2, 1, 0.5, 0 this gives {9, 5, 5, 11}.

### Scans (210–270)

- **210–225, `Run` and `Scan`.** `@dataclass` makes simple record classes with named fields.
  - `Run`: one simulation's folder, varied values, gravity model, seed and phase.
  - `Scan`: everything the figures need.
- **228–229, `short()`.** Removes the `cosmology_` prefix for folder names and labels.
- **232–270, `scan(...)`.**
  - **243:** `list(dict)` gives its keys, the names of the varied parameters.
  - **244–245, `vary`.** `itertools.product` gives every combination, e.g. 2 × 3 values → 6 tuples. There's no
    fiducial.
  - **246–249, `one_at_a_time`.**
    - The reference tuple holds each parameter's fiducial value, from `fixed` if given, else the base file.
    - Points: the reference first, then for each parameter i and each of its values, the reference with position i
      replaced. `reference[:i] + (value,) + reference[i+1:]` builds that tuple.
    - Values equal to the fiducial are skipped, to avoid a duplicate fiducial run.
  - **250:** with `gravity=["F5"]` the models are `["GR", "F5"]`; without it, `[None]`, meaning use the base file's
    model.
  - **252–268, the main loop over points.**
    - **253, `changed`.** The parameters that differ from the reference. With `vary` there is no reference, so
      everything counts as changed.
    - **254, `settings`.** `fixed` overridden by this point's values.
    - **255–260, input P(k).** If any cosmology parameter is involved, it computes one per cosmology:
      - The file name uses only the changed cosmology parameters (`pofk_h0.7006.txt`, or `pofk_fiducial.txt`).
      - Runs differing only in gravity, seed or phase share it, and the cache in `write_input_power` stops
        recomputation.
      - With `--plot-only` CLASS isn't run, but the path is still set, so the parameters (and thus folder names) stay
        consistent.
    - **261, folder label.** E.g. `OmegaCDM0.2295_h0.7`, or `fiducial`.
    - **262–263, gravity settings.** For each model, its parameters plus `ic_use_gravity_model_GR = True`, so all
      models start from identical initial conditions.
    - **264–268, runs.** For each seed (default: the base file's) and phase, calls `run()` with:
      - name `scans/<scan>/<label>_<model>_seed<seed>_<phase>`;
      - all the settings;
      - the seed;
      - `ic_reverse_phases = (phase == "reversed")`.

      Then records a `Run`.
  - **269–270:** prints the tally and returns the `Scan`.

### Command line (273–286)

`if __name__ == "__main__"` runs only when you execute `python3 runs.py`, not when another script imports it.
- **274–281, `parse()`.**
  - Turns `true`/`false` into booleans.
  - `ast.literal_eval` safely turns `1024`, `0.7` and `[1.0, 0.0]` into Python values.
  - Anything that isn't a Python literal, such as `f(R)`, stays a string.
- **283:** drops the `--` flags, then takes the first argument as the run name and the rest as assignments.
- **284:** `"key=value".split("=", 1)` splits at the first `=` only, giving the dict of parameters.
- **285–286:** runs and summarises. `ntasks=128` works too, because it matches `run()`'s own argument.

---

## 2. [spectra.py](scripts/spectra.py): loads and combines spectra (no plotting)

- **1–14, docstring.** Defines the two data formats used everywhere:
  - **curve:** a dict holding `k`, `P` with shape [15, nk] (all 15 spectra stacked), an optional `error`, and a
    `label`.
  - **figure:** a dict holding a title, curves (either `{values: curve}` with parameter `names`, or a plain list) and
    plot options.
- **24, `PAIRS`.** The 15 (i, j) with i ≤ j in file order: (0,0), (0,1) … (4,4). `P[p]` is spectrum `PAIRS[p]`.
- **25, `NAMES`.** The LaTeX names of the 5 fields, for panel titles.
- **28–35, `load_run(folder, z)`.**
  - `lru_cache` remembers every result, so a run used in several figures is read only once.
  - Finds `snapshot_*_z0.000`; if there's no such output, or it isn't complete (no `pofk_44.txt`), returns
    `(None, None)`.
  - Otherwise reads the 15 files. k is the first column of the first file; `P` stacks the second columns into a
    [15, nk] array.
- **38–41, `redshifts(folders)`.**
  - Collects all snapshot folder names of these runs.
  - `rsplit("_z", 1)[1]` takes the text after the last `_z`, e.g. `0.500`, as a float.
  - A set removes duplicates, sorted highest first.
- **44–52, `mean_curve(folders, z)`.**
  - Loads every folder and keeps the finished ones.
  - Mean over runs (`axis=0` averages across runs, leaving [15, nk]).
  - Standard deviation as `error` if more than one run.
  - `info` says how many went in.
- **55–80, `boosts(runs, model, z)`.** The f(R)/GR boost.
  - **60, GR lookup.** A table from (values, seed, phase) to the GR folder.
  - **61, `pairs`.** A two-level dictionary, values → seed → list of (GR, f(R)) spectra. `defaultdict` creates empty
    entries on first use.
  - **62–66, matching.** For every run of this model, finds its GR twin (same values, seed and phase). Only if both
    have finished is the pair stored.
  - **68–73, per seed.** Averages GR and f(R) over the phases first, then divides. This is the step that cancels the
    sign-flipping odd spectra. `np.where(signs equal, ratio, NaN)` leaves out points where the two have opposite
    signs.
  - **74–79, mean over seeds.**
    - `nanmean`/`nanstd` ignore the NaNs.
    - The error is the standard error, std/√n with `ddof=1` (sample standard deviation).
    - The warning filter silences numpy's "mean of empty slice" where every seed is NaN.
- **83–107, `scan_figures(scan)`.** For each redshift:
  - **90–94, spectra.** For each point in run order, the mean over its seeds and phases in the first model (GR).
    `dict.fromkeys` removes duplicates while keeping the order.
  - **95:** a figure of these spectra.
  - **96–103, ratios.**
    - With at least 2 points, divides everything by the fiducial (or, without one, by the smallest values tuple).
    - `np.interp` evaluates the reference at each curve's k, so it also works if the k grids differ.
    - The title says what it is relative to.
  - **104–106:** one boost figure per f(R) model.
- **110–133, `pages_varying_one_parameter(curves, names, max_lines=5)`.** Used only by the 15-panel pages.
  - For each parameter i, groups the curves by the values of all other parameters, so within a group only parameter
    i varies.
  - Groups with fewer than 2 curves have nothing to compare and are skipped.
  - Big groups are cut into chunks of 5. Each chunk is a page titled "Varying h, OmegaCDM = …", with each curve
    labelled by its varying value.
  - `shown` remembers what appeared, so curves that never shared a page (differing in 2+ parameters from all others)
    go on "Other" pages at the end, with full labels.
- **136–139, `label(values, names, fiducial)`.** "f = 512, n = 1024". With fiducial values it lists only the
  parameters that differ, or "fiducial" if none do.
- **142–163, `group_figures(groups, zs, ratio)`.** The command-line case.
  - For each redshift and each glob pattern: the mean curve, labelled with the pattern and run count.
  - With `ratio`, removes the first group as reference and divides the others by it, interpolated. The error is
    scaled by |P_ref|.

---

## 3. [plot_page.py](scripts/plot_page.py): 15 spectra on one page

- **OPTIONS.** The figure keys that are plot options, copied into function arguments.
- **`plot_page(curves, title, ratio, kmin, kmax, line, ylabel)`.**
  - **Grid.** Builds a 5×5 grid (`sharex` links the x axes) and turns every panel off; only the upper triangle is
    turned back on later.
  - **Each curve.**
    - The colour is the n-th of matplotlib's 10 default colours, unless the curve brings its own (`color`/`ls`, e.g.
      the black reference).
    - `shown` masks the k range.
    - Each of the 15 spectra goes into panel (i, j), as |P| unless it's a ratio, with the label; the legend later
      collects these.
    - The error band is P ± error.
  - **Panels.** For each of the 15: log x axis and title δ × δ. Tick labels only on the bottom panel of each column
    (the diagonal, i == j). Ratios get a horizontal line and `set_ratio_range`; spectra get a log y axis.
  - **Legend and labels.** The legend is taken from the first panel, since every panel has the same lines, and placed
    in the empty lower-left triangle. Then the overall title and axis labels.
- **`set_ratio_range(ax, line)`.** Sets the y range to the 2nd–98th percentile of all line values, always including
  the reference line, plus 5% margins. One spike at a zero crossing therefore can't squash the panel.
  `get_lines()[:-1]` leaves out the horizontal line itself, which is the last one added.
- **`save(figures, filename, verbose)`.** Writes the figures as pages of `figures/<filename>`, closing each to free
  memory. Does nothing if empty.
- **`save_pages(figures, filename)`.** For each figure, splits its curves:
  - `{values: curve}` goes through `pages_varying_one_parameter`;
  - a plain list is cut into chunks of 5.

  Every page gets the figure's `extra` curves (e.g. the reference) and the joined title. Everything goes into one
  PDF.
- **`__main__`.** `argparse` reads `OUTPUT GROUP... [--ratio] [--z ...]` and calls `group_figures` → `save_pages`.

---

## 4. [plot_single.py](scripts/plot_single.py): one figure per spectrum

- **18, `STYLE`.** Serif fonts with Computer Modern maths, like the old LaTeX plots, applied only here.
- **19, `LINESTYLES`.** After 10 colours the line style changes, so up to 40 curves stay distinguishable.
- **22–49, `plot_single(curves, p, ...)`.** One spectrum, `PAIRS[p]`, for all curves:
  - **28–29, style.** Colour `n % 10` and line style `n // 10`, unless the curve sets its own.
  - **30–36:** the line and a faint band (`alpha=0.05`, so many bands don't merge into a blob).
  - **37–42:** log or ratio axes, as in `plot_page`.
  - **43:** the title, e.g. "P₀₂: δ_m × δ_L²", plus the figure title on a second line.
  - **44–48:** axis labels, the legend outside on the right (`bbox_to_anchor=(1.02, 1)` puts it just right of the
    axes), a light grid, and `tight_layout` so the legend fits.
- **52–68, `save_single(figures, folder)`.**
  - **55–62.** For each figure, turns `{values: curve}` into a labelled list (only what differs from the fiducial)
    and keeps the options. Unlike `save_pages`, nothing is split: all curves stay together.
  - **63–66.** Inside the style context, for each of the 15 spectra, makes one page per figure and saves
    `figures/<folder>/P_ij.pdf`.
  - **67–68:** prints one summary line.
- **71–78, `__main__`.** Same interface as `plot_page.py`, but OUTPUT is a folder.

---

## 5. The run scripts

They all follow the same pattern:
1. imports;
2. settings in CAPITALS;
3. loops calling `run()`, or one `scan()`;
4. `summary()`;
5. building figure dicts;
6. `save_pages` + `save_single`.

### [run_convergence.py](scripts/run_convergence.py)

- **14:** the 12 combinations of f, n ∈ {512, 1024} and s ∈ {10, 20, 40}, plus (1280, 1280, 40).
- **17–21:** for each setup, 10 runs. 128 tasks for n > 1024 (memory); seed 1230 + version. The folders are stored per
  setup.
- **25:** adds the hand-made 2048 run by its path. It is not passed to `run()`, so it can never be started by
  accident.
- **26:** the mean curve per setup. `:=` (walrus) assigns and tests in one go, which drops setups with no finished
  runs.
- **27, best setup.** `max` over tuples compares f first, then n, then s, so the highest resolution wins.
- **28:** every setup divided by the best one.
- **29–30:** the 2048³ reference as a solid black curve.
- **33–37:** two figures. The spectra have the reference as `extra` on every page; then the ratios.

### [run_weekend.py](scripts/run_weekend.py)

- **18–21:** models, redshifts, seeds, and the list that collects `Run`s.
- **24–33, `run_models`.**
  - For the 4 models × 2 phases: the old folder-name scheme, all the parameters, and the step split from
    `steps_per_interval`.
  - It records a `Run` with values (L, f, n, s), so `boosts` can pair f(R) with GR.
- **36:** the non-flag arguments choose `ensemble` and/or `steps`; default both.
- **37–40:** for each seed, the 1024 box, plus 2048 or 512 alternately (`i % 2`).
- **41–43:** 20/30/40 steps for the first seed.
- **46–49:** one figure per redshift, low z first. Its curves come from `boosts()` per f(R) model, keyed
  `(model, L, f, n, s)` so pages can vary any of those. Only k = 0.1–1 is plotted.

### [run_gr_vs_fofr.py](scripts/run_gr_vs_fofr.py)

The same pattern for 4 resolutions × 2 phases × {GR, f(R)}, seed 1234. The f(R) runs keep the base file's fofr0
(1e-5), so they're recorded as `"F5"`.

### [run_ensemble_100.py](scripts/run_ensemble_100.py)

100 runs at 512³ and 10 steps, seeds 12301–12400. One mean curve with its scatter.

### [run_resolution_scan.py](scripts/run_resolution_scan.py)

27 runs, f × L × n. One curve per setup, with pages varying f, L or n.

### [scan_planck.py](scripts/scan_planck.py)

- **31–33:** the Planck 2018 numbers. Ω_b = ω_b/h².
- **35–37, `PLANCK`.** The fiducial w0waCDM cosmology, at module level so other scripts can import it.
  A_s = exp(3.044) × 10⁻¹⁰, rounded to 6 digits.
- **39:** `__main__` guard; importing `PLANCK` doesn't start the scan.
- **40–46:** one-at-a-time steps of −5, −1, +1, +5 σ:
  - Ω_m through Ω_cdm at fixed Ω_b;
  - A_s in ln A_s;
  - wa capped at +3σ.

  `gravity=["F5"]` gives GR + F5 for each.
- **47–49:** figures, then both plotters.

### [scan_cosmology.py](scripts/scan_cosmology.py)

The same pattern around the base file's cosmology: ±10% and ±50% in Ω_m and A_s, two seeds.

### [scan_planck_mesh.py](scripts/scan_planck_mesh.py)

Imports `PLANCK` and scans `pofk_nmesh = [1024]`, GR only, normal phase. That was the single mesh test.

### [plot_fixes.py](scripts/plot_fixes.py)

No runs.
- **`phase_mean`:** the fiducial GR or F5 spectra, averaged over the phases.
- **`with_laplacian_rule`:** replaces the 5 spectra with index 4 (∇²δ_L) by −k² P_{0,i}, and ∇²δ_L × ∇²δ_L by
  k⁴ P₀₀.
- **`savgol`:** for each point, fits a cubic through the 11 points around it (shifted inward at the ends) and takes
  its value at that point.
- **`normalised_difference`:** ΔP_ij / √(P_ii P_jj), using the GR auto spectra.
- **`measured_and_fixed`:** builds the blue "measured" and red "fixed" curve pair.
- **Figures:** three pages (spectra, boosts, normalised differences), saved with both plotters.

---

## Changing parameters: what applies where

- Every run starts from a copy of [parameterfile.lua](parameterfile.lua), with only the parameters passed to `run()`
  (or a scan's `fixed`/`vary`/`one_at_a_time`) replaced. An edit to the base file therefore applies to every new run,
  unless the run or scan sets that parameter itself. Several campaigns pin their settings (e.g. the ensemble at 512³
  and 10 steps).
- Finished runs are skipped, and folder names only hold the varied parameters. After editing the base file, use a
  new scan or campaign name, or a campaign ends up with mixed settings.
- Changing cosmology only in the base file does not recompute the input P(k). CLASS runs only when a run or scan sets
  `cosmology_*` parameters itself; otherwise the runs keep the base file's `ic_input_filename`. Set cosmology
  parameters in the run or scan, e.g. `fixed={"cosmology_h": 0.7}`.
- Replacement is by name, anywhere in the file, including indented lines inside `if` blocks. A misspelled name
  matches nothing and is silently ignored.
