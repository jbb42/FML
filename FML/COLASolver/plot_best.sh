#!/bin/bash

# Define absolute base path
BASE_DIR="/mn/stornext/u3/jonasbbe/pc/Dokumenter/FML/FML/COLASolver"
SPECTRA_DIR="${BASE_DIR}/saved_spectra2" # Central directory to save spectra

# Create output directory
mkdir -p "$SPECTRA_DIR"

# Fixed parameters
P_NMESH=512
BOXSIZE=1024

# Parameter grids
VALUES=(1280)
STEPS_VALUES=(40) # Adjust if you only want 10
VERSIONS=10             # Number of runs per combination for averaging

# Base Lua configuration file
BASE_PARAM_FILE="${BASE_DIR}/parameterfile.lua"

for f_nmesh in "${VALUES[@]}"; do
  for npart in "${VALUES[@]}"; do
    for nsteps in "${STEPS_VALUES[@]}"; do
      for run_idx in $(seq 1 $VERSIONS); do
        
        # Generate a unique seed for this realization
        SEED=$(( run_idx + 1230 ))
        
        # Define specific output folder
        RUN_ID="f${f_nmesh}_p${P_NMESH}_b${BOXSIZE}_n${npart}_s${nsteps}_v${run_idx}"
        OUT_DIR="${BASE_DIR}/output_${RUN_ID}"
        SNAP_DIR="${OUT_DIR}/snapshot_TestSim_z0.000/"
        TEMP_LUA="${BASE_DIR}/temp_param_${RUN_ID}.lua"
        
        echo "=================================================="
        echo "Preparing parameter file for: $RUN_ID"
        
        # Update Lua parameters on the fly, including Lua array syntax for timesteps
        sed -e "s/^[[:space:]]*force_nmesh[[:space:]]*=.*/force_nmesh = ${f_nmesh}/" \
            -e "s/^[[:space:]]*pofk_nmesh[[:space:]]*=.*/pofk_nmesh = ${P_NMESH}/" \
            -e "s/^[[:space:]]*simulation_boxsize[[:space:]]*=.*/simulation_boxsize = ${BOXSIZE}/" \
            -e "s/^[[:space:]]*particle_Npart_1D[[:space:]]*=.*/particle_Npart_1D = ${npart}/" \
            -e "s/^[[:space:]]*timestep_nsteps[[:space:]]*=.*/timestep_nsteps = {${nsteps}}/" \
            -e "s/^[[:space:]]*ic_random_seed[[:space:]]*=.*/ic_random_seed = ${SEED}/" \
            -e "s|^[[:space:]]*output_folder[[:space:]]*=.*|output_folder = \"${OUT_DIR}\"|" \
            "$BASE_PARAM_FILE" > "$TEMP_LUA"

        echo "Running FML nbody for: $RUN_ID"
        
        # 1. Run FML nbody solver
        cd "$BASE_DIR" || exit
        mpirun -np 128 ./nbody "$TEMP_LUA"
        
        MPI_STATUS=$?
            
        if [ $MPI_STATUS -eq 0 ]; then
            # 2. Save the spectra text files directly
            echo "Saving raw spectra for $RUN_ID..."
            RUN_SPECTRA_DIR="${SPECTRA_DIR}/${RUN_ID}"
            mkdir -p "$RUN_SPECTRA_DIR"
            cp "${SNAP_DIR}"pofk_*.txt "$RUN_SPECTRA_DIR/"

            # 3. Clean up the heavy output data folder to save disk space
            echo "Cleaning up heavy snapshot data in ${OUT_DIR}..."
            rm -rf "$OUT_DIR"
        
        else
            echo "=================================================="
            echo "ERROR: Run $RUN_ID crashed with MPI exit code $MPI_STATUS!"
            echo "Skipping file save. Raw data folder kept for debugging."
            echo "=================================================="
        fi

        # 4. Clean up the temporary Lua file
        rm -f "$TEMP_LUA"
        
      done
    done
  done
done

echo "Batch execution complete."
