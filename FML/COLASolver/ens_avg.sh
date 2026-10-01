#!/bin/bash

PARAM_FILE="parameterfile.lua"
BASE_SEED=12300
ENSEMBLE_SIZE=100

# ==========================================
# PRE-FLIGHT CHECKS
# ==========================================
echo "Running pre-flight checks..."

# 1. Check if the parameter file exists
if [ ! -f "$PARAM_FILE" ]; then
    echo "ERROR: Parameter file '$PARAM_FILE' not found."
    exit 1
fi

# 2. Check if the executable is present and has run permissions
if [ ! -x "./nbody" ]; then
    echo "ERROR: Executable './nbody' not found or is missing 'chmod +x' permissions."
    exit 1
fi

# 3. Check write permissions for the base output directory
mkdir -p output
if [ ! -w "output" ]; then
    echo "ERROR: Cannot write to the 'output' directory. Check your storage quota."
    exit 1
fi

echo "All checks passed. Starting ensemble..."
echo "=========================================="

# ==========================================
# MAIN ENSEMBLE LOOP
# ==========================================
for i in $(seq 1 $ENSEMBLE_SIZE); do
    CURRENT_SEED=$((BASE_SEED + i))
    echo "Starting Realization $i/$ENSEMBLE_SIZE (Seed: $CURRENT_SEED)"
    
    # Update seed in the Lua parameter file
    sed -i "s/ic_random_seed = .*/ic_random_seed = $CURRENT_SEED/" $PARAM_FILE
    
    # Run solver
    mpirun -np 64 ./nbody $PARAM_FILE
    
    # TRIPWIRE: Verify FML actually created the snapshot directory and files
    # If the simulation crashed (e.g., MPI buffer overflow), this catches it and aborts.
    if ! ls output/snapshot_TestSim_z*/*.txt 1> /dev/null 2>&1; then
        echo "CRITICAL ERROR: No text output detected for Seed $CURRENT_SEED."
        echo "The FML solver likely crashed. Aborting the ensemble loop."
        exit 1
    fi
    
    # Safely create the seed directory and move the data
    mkdir -p output_seed_${i}
    mv output/snapshot_TestSim_z*/*.txt output_seed_${i}/
    
    # Clean up the empty snapshot directory for the next run
    rm -rf output/snapshot_TestSim_z*
done

echo "100-seed ensemble completed successfully."
