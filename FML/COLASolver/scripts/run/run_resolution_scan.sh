#!/bin/bash

# Define absolute base path to avoid relative path confusion
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)" # COLASolver root
PLOTS_DIR="${BASE_DIR}/figures/resolution_scan"

# Create a central directory for all plots
mkdir -p $PLOTS_DIR

# Fixed parameter
P_NMESH=512

# Parameter values: powers of two from 256 to 1024
VALUES=(256 512 1024)

# Your base Lua configuration file
BASE_PARAM_FILE="${BASE_DIR}/parameterfile.lua"

for f_nmesh in "${VALUES[@]}"; do
  for boxsize in "${VALUES[@]}"; do
    for npart in "${VALUES[@]}"; do
      
      # Define specific output folder for this parameter combination
      RUN_ID="f${f_nmesh}_p${P_NMESH}_b${boxsize}_n${npart}"
      OUT_DIR="${BASE_DIR}/output_${RUN_ID}"
      
      SNAP_DIR="${OUT_DIR}/snapshot_TestSim_z0.000/"
      TEMP_LUA="${BASE_DIR}/temp_param_${RUN_ID}.lua"
      
      echo "=================================================="
      echo "Preparing parameter file for: $RUN_ID"
      
      # Use sed to copy the base file and replace the variables on the fly
      # NOTE: Ensure 'output_folder' matches the exact variable name in your lua file!
      sed -e "s/^[[:space:]]*force_nmesh[[:space:]]*=.*/force_nmesh = ${f_nmesh}/" \
          -e "s/^[[:space:]]*pofk_nmesh[[:space:]]*=.*/pofk_nmesh = ${P_NMESH}/" \
          -e "s/^[[:space:]]*simulation_boxsize[[:space:]]*=.*/simulation_boxsize = ${boxsize}/" \
          -e "s/^[[:space:]]*particle_Npart_1D[[:space:]]*=.*/particle_Npart_1D = ${npart}/" \
          -e "s|^[[:space:]]*output_folder[[:space:]]*=.*|output_folder = \"${OUT_DIR}\"|" \
          $BASE_PARAM_FILE > $TEMP_LUA

      echo "Running FML nbody for: $RUN_ID"
      
      # 1. Run FML nbody solver
      cd $BASE_DIR
      mpirun -np 64 ./nbody $TEMP_LUA
      
      # Capture the exit status of mpirun to detect crashes (like the node 42 buffer error)
      MPI_STATUS=$?
          
      if [ $MPI_STATUS -eq 0 ]; then
          echo "Generating plots for $RUN_ID..."

          # 2. Plot: Calculated vs Theoretical (k^2 scaled) cross-spectra
          gnuplot <<- EOF
              set terminal pngcairo size 1800,1200 font "Sans,12"
              set output "${PLOTS_DIR}/plot_k2_comparison_${RUN_ID}.png"
              set logscale xy
              set xlabel 'Wavenumber k [h/Mpc]'
              set ylabel 'Power P_{ij}(k) [(Mpc/h)^3]'
              set key top right box outside
              
              R_star = 0.75
              R2 = R_star**2
              dir = '${SNAP_DIR}'
              
              plot \
                  dir.'pofk_04.txt' using 1:(abs(\$2)) with lines lw 2 lc 'red' title 'Calculated |P_{04}|', \
                  dir.'pofk_01.txt' using 1:(abs((\$1**2) * R2 * \$2)) with lines dt 2 lw 2 lc 'red' title 'Theory k^2 R_*^2 |P_{01}|', \
                  \
                  dir.'pofk_14.txt' using 1:(abs(\$2)) with lines lw 2 lc 'blue' title 'Calculated |P_{14}|', \
                  dir.'pofk_11.txt' using 1:(abs((\$1**2) * R2 * \$2)) with lines dt 2 lw 2 lc 'blue' title 'Theory k^2 R_*^2 P_{11}', \
                  \
                  dir.'pofk_24.txt' using 1:(abs(\$2)) with lines lw 2 lc 'forest-green' title 'Calculated |P_{24}|', \
                  dir.'pofk_12.txt' using 1:(abs((\$1**2) * R2 * \$2)) with lines dt 2 lw 2 lc 'forest-green' title 'Theory k^2 R_*^2 |P_{12}|', \
                  \
                  dir.'pofk_34.txt' using 1:(abs(\$2)) with lines lw 2 lc 'orange' title 'Calculated |P_{34}|', \
                  dir.'pofk_13.txt' using 1:(abs((\$1**2) * R2 * \$2)) with lines dt 2 lw 2 lc 'orange' title 'Theory k^2 R_*^2 |P_{13}|', \
                  \
                  dir.'pofk_44.txt' using 1:(abs(\$2)) with lines lw 2 lc 'purple' title 'Calculated P_{44}', \
                  dir.'pofk_11.txt' using 1:(abs((\$1**4) * (R2**2) * \$2)) with lines dt 2 lw 2 lc 'purple' title 'Theory k^4 R_*^4 P_{11}'
EOF

          # 3. Plot: All exact fields
          gnuplot <<- EOF
              set terminal pngcairo size 1800,1200 font "Sans,12"
              set output "${PLOTS_DIR}/plot_exact_fields_${RUN_ID}.png"
              set logscale xy
              set xlabel 'Wavenumber k [h/Mpc]'
              set ylabel 'Power P_{ij}(k) [(Mpc/h)^3]'
              set key top right box outside
              
              dir = '${SNAP_DIR}'
              
              plot \
                  dir.'pofk_00.txt' using 1:2 with lines lw 2 title 'P_{00} (Matter Auto)', \
                  dir.'pofk_01.txt' using 1:2 with lines lw 2 title 'P_{01} (Matter - \\delta_L)', \
                  dir.'pofk_11.txt' using 1:2 with lines lw 2 title 'P_{11} (\\delta_L Auto)', \
                  dir.'pofk_02.txt' using 1:2 with lines lw 2 title 'P_{02} (Matter - \\delta_L^2)', \
                  dir.'pofk_12.txt' using 1:(abs(\$2)) with lines lw 2 title '|P_{12}| (\\delta_L - \\delta_L^2)', \
                  dir.'pofk_22.txt' using 1:2 with lines lw 2 title 'P_{22} (\\delta_L^2 Auto)', \
                  dir.'pofk_03.txt' using 1:(abs(\$2)) with lines lw 2 dt 2 title '|P_{03}| (Matter - s^2)', \
                  dir.'pofk_13.txt' using 1:(abs(\$2)) with lines lw 2 dt 2 title '|P_{13}| (\\delta_L - s^2)', \
                  dir.'pofk_23.txt' using 1:2 with lines lw 2 dt 2 title 'P_{23} (\\delta_L^2 - s^2)', \
                  dir.'pofk_33.txt' using 1:2 with lines lw 2 dt 2 title 'P_{33} (s^2 Auto)', \
                  dir.'pofk_04.txt' using 1:(abs(\$2)) with lines lw 2 dt 3 title '|P_{04}| (Matter - \\nabla^2\\delta_L)', \
                  dir.'pofk_14.txt' using 1:(abs(\$2)) with lines lw 2 dt 3 title '|P_{14}| (\\delta_L - \\nabla^2\\delta_L)', \
                  dir.'pofk_24.txt' using 1:(abs(\$2)) with lines lw 2 dt 3 title '|P_{24}| (\\delta_L^2 - \\nabla^2\\delta_L)', \
                  dir.'pofk_34.txt' using 1:(abs(\$2)) with lines lw 2 dt 3 title '|P_{34}| (s^2 - \\nabla^2\\delta_L)', \
                  dir.'pofk_44.txt' using 1:2 with lines lw 2 dt 3 title 'P_{44} (\\nabla^2\\delta_L Auto)'
EOF

          # 4. Clean up the heavy output data folder to save disk space
          echo "Cleaning up raw data in ${OUT_DIR}..."
          rm -rf $OUT_DIR
      
      else
          echo "=================================================="
          echo "ERROR: Run $RUN_ID crashed with MPI exit code $MPI_STATUS!"
          echo "Skipping plotting. Raw data folder kept for debugging."
          echo "=================================================="
      fi

      # 5. Clean up the temporary Lua file
      rm -f $TEMP_LUA
      
    done
  done
done

echo "Batch execution complete."
