#!/bin/bash
# GR vs f(R) bias spectra at several resolutions, as phase-reversed pairs, all from the same seed.
#
#   resolutions (force_nmesh, particle_Npart_1D) x ic_reverse_phases {false,true} x gravity {GR, f(R)}
#
# Both gravity models start from identical initial conditions (ic_use_gravity_model_GR = true): the
# LCDM input P(k) is scaled back to the initial redshift with GR growth in both runs.
# f(R) parameters (fR0, n, screening) are taken from parameterfile.lua.
#
# Usage:  scripts/run/run_gr_vs_fofr.sh            run everything (finished runs are skipped)
#         DRY_RUN=1 scripts/run/run_gr_vs_fofr.sh  only write the parameter files and print the commands

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)" # LagrangianBias root
NBODY="$(cd "${BASE_DIR}/../COLASolver" && pwd)/nbody" # The COLASolver executable
RESULTS_DIR="${BASE_DIR}/results/gr_vs_fofr"
BASE_PARAM_FILE="${BASE_DIR}/parameterfile.lua"

# Settings shared by all runs
SEED=1234
NSTEPS=30
FIX_AMPLITUDE=true
BOXSIZE=1024
P_NMESH=512

# Resolutions as "force_nmesh:particle_Npart_1D", and MPI tasks per particle number
RESOLUTIONS=("512:512" "1024:512" "1024:1024" "1280:1280")
declare -A NTASKS=([512]=64 [1024]=64 [1280]=128)

REVERSE_PHASES=(false true)
GRAVITY_MODELS=("GR" "f(R)")

mkdir -p "$RESULTS_DIR"
cd "$BASE_DIR" || exit 1 # input paths in the parameter file are relative to the root

nfailed=0
for res in "${RESOLUTIONS[@]}"; do
  fmesh=${res%%:*}
  npart=${res##*:}
  for rev in "${REVERSE_PHASES[@]}"; do
    for gravity in "${GRAVITY_MODELS[@]}"; do

      gravity_label=$([ "$gravity" == "GR" ] && echo "gr" || echo "fofr")
      fix_label=$([ "$FIX_AMPLITUDE" == "true" ] && echo "fixed" || echo "gauss")
      rev_label=$([ "$rev" == "true" ] && echo "reversed" || echo "normal")
      RUN_ID="${gravity_label}_f${fmesh}_n${npart}_s${NSTEPS}_${fix_label}_${rev_label}_seed${SEED}"
      OUT_DIR="${RESULTS_DIR}/${RUN_ID}"
      PARAM_FILE="${OUT_DIR}/parameterfile.lua"

      echo "=================================================="
      if [ -f "${OUT_DIR}/snapshot_TestSim_z0.000/pofk_bias_info.txt" ]; then
        echo "Skipping $RUN_ID (already finished)"
        continue
      fi
      echo "Preparing $RUN_ID"
      mkdir -p "$OUT_DIR"

      # The parameter file is kept next to the results for reference
      sed -e "s/^[[:space:]]*gravity_model[[:space:]]*=.*/gravity_model = \"${gravity}\"/" \
          -e "s/^[[:space:]]*ic_use_gravity_model_GR[[:space:]]*=.*/ic_use_gravity_model_GR = true/" \
          -e "s/^[[:space:]]*ic_fix_amplitude[[:space:]]*=.*/ic_fix_amplitude = ${FIX_AMPLITUDE}/" \
          -e "s/^[[:space:]]*ic_reverse_phases[[:space:]]*=.*/ic_reverse_phases = ${rev}/" \
          -e "s/^[[:space:]]*ic_random_seed[[:space:]]*=.*/ic_random_seed = ${SEED}/" \
          -e "s/^[[:space:]]*force_nmesh[[:space:]]*=.*/force_nmesh = ${fmesh}/" \
          -e "s/^[[:space:]]*particle_Npart_1D[[:space:]]*=.*/particle_Npart_1D = ${npart}/" \
          -e "s/^[[:space:]]*pofk_nmesh[[:space:]]*=.*/pofk_nmesh = ${P_NMESH}/" \
          -e "s/^[[:space:]]*simulation_boxsize[[:space:]]*=.*/simulation_boxsize = ${BOXSIZE}/" \
          -e "s/^[[:space:]]*timestep_nsteps[[:space:]]*=.*/timestep_nsteps = {${NSTEPS}}/" \
          -e "s|^[[:space:]]*output_folder[[:space:]]*=.*|output_folder = \"${OUT_DIR}\"|" \
          "$BASE_PARAM_FILE" > "$PARAM_FILE"

      CMD=(mpirun -np "${NTASKS[$npart]}" "$NBODY" "$PARAM_FILE")
      if [ -n "$DRY_RUN" ]; then
        echo "DRY_RUN: ${CMD[*]}"
        continue
      fi

      echo "Running $RUN_ID with ${NTASKS[$npart]} tasks"
      "${CMD[@]}" > "${OUT_DIR}/log.txt" 2>&1
      MPI_STATUS=$?

      if [ $MPI_STATUS -eq 0 ] && [ -f "${OUT_DIR}/snapshot_TestSim_z0.000/pofk_bias_info.txt" ]; then
        echo "Finished $RUN_ID"
      else
        echo "ERROR: $RUN_ID failed (MPI exit code $MPI_STATUS), see ${OUT_DIR}/log.txt"
        nfailed=$((nfailed + 1))
      fi

    done
  done
done

echo "=================================================="
echo "All runs done ($nfailed failed). Results in $RESULTS_DIR"
