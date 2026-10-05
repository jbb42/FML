#!/bin/bash
# Weekend campaign: GR vs f(R) bias spectra, all with fixed amplitudes and phase-reversed pairs, and
# identical initial conditions for every gravity model (ic_use_gravity_model_GR = true).
# Every run has f = n = 1024 and 30 steps unless stated otherwise, with outputs at z = 2, 1, 0.5, 0.
#
#   ensemble: 10 rounds, one seed each, GR, F4, F5 and F6, 30 steps:
#               1024 Mpc/h every round (production), plus 2048 Mpc/h (low k) on odd and
#               512 Mpc/h (high k) on even rounds, for stitching the boxes together
#   steps   : GR, F4, F5 and F6 at 20, 30 and 40 steps in 1024 Mpc/h, first seed only
#
# Every seed is run as a pair (normal and reversed phases). Simulations run one at a time, in this order.
# The total is sized to take roughly a weekend (~30-42 h). Runs shared between campaigns are only done once
# (the 30-step runs of the first seed are part of both), and finished runs are skipped on restart.
#
# Usage:  scripts/run/run_weekend.sh [campaign ...]                 default: ensemble steps
#         DRY_RUN=1 scripts/run/run_weekend.sh [campaign ...]       only write the parameter files and list the runs
#         tail -F results/weekend/current_log.txt                   follow the running simulation, switching automatically

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)" # COLASolver root
RESULTS_DIR="${BASE_DIR}/results/weekend"
BASE_PARAM_FILE="${BASE_DIR}/parameterfile.lua"

# Setup shared by all runs
N=1024
NTASKS=64
NTHREADS=2 # OpenMP threads per MPI task (all runs fit on one node, so the ranks inherit these)
export OMP_NUM_THREADS=$NTHREADS
export I_MPI_PIN_DOMAIN=omp # Intel MPI: give each rank its own NTHREADS cores (ignored by other MPIs)
PROD_BOX=1024
PROD_NSTEPS=30
POFK_NMESH=512 # k_Nyquist = 0.79, 1.57, 3.14 h/Mpc in the 2048, 1024, 512 Mpc/h boxes
OUTPUT_REDSHIFTS="2.0, 1.0, 0.5, 0.0"
SEEDS=($(seq 1001 1010)) # one ensemble round per seed
GRAVITY_MODELS=(GR F4 F5 F6)

mkdir -p "$RESULTS_DIR"
cd "$BASE_DIR" || exit 1 # input/ paths in the parameter file are relative to the root

#==================================================================================
# Split nsteps over the output intervals in proportion to the change in a, so the step size
# matches a single-output run with the same total (timestep_scalefactor_spacing = "linear")
#==================================================================================
steps_per_interval() {
    local nsteps=$1 zini=$2
    echo "$OUTPUT_REDSHIFTS" | awk -v n="$nsteps" -v zini="$zini" -F', *' '{
        aprev = 1.0 / (1.0 + zini); tot = 1.0 - aprev
        for (i = 1; i <= NF; i++) { a = 1.0 / (1.0 + $i); x[i] = n * (a - aprev) / tot; s[i] = int(x[i]); used += s[i]; aprev = a }
        # Largest remainder rounding so the steps add up to n
        while (used < n) { best = 1; for (i = 1; i <= NF; i++) if (x[i] - s[i] > x[best] - s[best]) best = i; s[best]++; used++; x[best] = s[best] }
        out = s[1]; for (i = 2; i <= NF; i++) out = out ", " s[i]; print out
    }'
}
ZINI=$(grep -E "^[[:space:]]*ic_initial_redshift[[:space:]]*=" "$BASE_PARAM_FILE" | sed 's/.*=[[:space:]]*\([0-9.eE+-]*\).*/\1/')

#==================================================================================
# Run one simulation: run_sim <gravity GR|F4|F5|F6> <box> <nsteps> <reverse true|false> <seed>
#==================================================================================
START_TIME=$(date +%s)
nrun=0; nskip=0; nfailed=0
run_sim() {
    local grav=$1 box=$2 nsteps=$3 rev=$4 seed=$5
    local gravity_model="GR" fr0_line="" glabel="gr"
    if [ "$grav" != "GR" ]; then
        gravity_model="f(R)"; glabel="fofr${grav}"
        fr0_line="s/^[[:space:]]*gravity_model_fofr_fofr0[[:space:]]*=.*/  gravity_model_fofr_fofr0 = 1e-${grav#F}/"
    fi
    local rlabel=$([ "$rev" == "true" ] && echo "reversed" || echo "normal")
    local run_id="${glabel}_L${box}_f${N}_n${N}_s${nsteps}_fixed_${rlabel}_seed${seed}"
    local out_dir="${RESULTS_DIR}/${run_id}"
    local param_file="${out_dir}/parameterfile.lua"

    if [ -f "${out_dir}/snapshot_TestSim_z0.000/pofk_bias_info.txt" ]; then
        nskip=$((nskip + 1)); return
    fi

    mkdir -p "$out_dir"

    sed -e "s/^[[:space:]]*gravity_model[[:space:]]*=.*/gravity_model = \"${gravity_model}\"/" \
        -e "${fr0_line:-s/^\$//}" \
        -e "s/^[[:space:]]*ic_use_gravity_model_GR[[:space:]]*=.*/ic_use_gravity_model_GR = true/" \
        -e "s/^[[:space:]]*ic_fix_amplitude[[:space:]]*=.*/ic_fix_amplitude = true/" \
        -e "s/^[[:space:]]*ic_reverse_phases[[:space:]]*=.*/ic_reverse_phases = ${rev}/" \
        -e "s/^[[:space:]]*ic_random_seed[[:space:]]*=.*/ic_random_seed = ${seed}/" \
        -e "s/^[[:space:]]*simulation_boxsize[[:space:]]*=.*/simulation_boxsize = ${box}/" \
        -e "s/^[[:space:]]*force_nmesh[[:space:]]*=.*/force_nmesh = ${N}/" \
        -e "s/^[[:space:]]*particle_Npart_1D[[:space:]]*=.*/particle_Npart_1D = ${N}/" \
        -e "s/^[[:space:]]*pofk_nmesh[[:space:]]*=.*/pofk_nmesh = ${POFK_NMESH}/" \
        -e "s/^[[:space:]]*output_redshifts[[:space:]]*=.*/output_redshifts = {${OUTPUT_REDSHIFTS}}/" \
        -e "s/^[[:space:]]*timestep_nsteps[[:space:]]*=.*/timestep_nsteps = {$(steps_per_interval "$nsteps" "$ZINI")}/" \
        -e "s|^[[:space:]]*output_folder[[:space:]]*=.*|output_folder = \"${out_dir}\"|" \
        "$BASE_PARAM_FILE" > "$param_file"

    local cmd=(mpirun -np "$NTASKS" ./nbody "$param_file")
    nrun=$((nrun + 1))
    if [ -n "$DRY_RUN" ]; then
        echo "DRY_RUN [$CAMPAIGN] ${cmd[0]} ${cmd[1]} ${cmd[2]} ${cmd[3]} $run_id"
        return
    fi

    echo "$(date '+%F %T') [$CAMPAIGN] Running $run_id"
    # Always points to the running simulation's log: follow it with  tail -F results/weekend/current_log.txt
    ln -sfn "${out_dir}/log.txt" "${RESULTS_DIR}/current_log.txt"
    local t0=$(date +%s)
    "${cmd[@]}" > "${out_dir}/log.txt" 2>&1
    local status=$?
    local duration=$(( $(date +%s) - t0 ))
    if [ $status -eq 0 ] && [ -f "${out_dir}/snapshot_TestSim_z0.000/pofk_bias_info.txt" ]; then
        echo "$(date '+%F %T') [$CAMPAIGN] Finished $run_id in $((duration / 60)) min"
    else
        echo "$(date '+%F %T') [$CAMPAIGN] ERROR: $run_id failed (MPI exit code $status), see ${out_dir}/log.txt"
        nfailed=$((nfailed + 1))
    fi
}

# Both phases of a seed for a list of gravity models: run_pairs <box> <nsteps> <seed> <gravity...>
run_pairs() {
    local box=$1 nsteps=$2 seed=$3; shift 3
    for grav in "$@"; do
        for rev in false true; do
            run_sim "$grav" "$box" "$nsteps" "$rev" "$seed"
        done
    done
}

#==================================================================================
# Campaigns
#==================================================================================
campaign_ensemble() {
    local i=0
    for seed in "${SEEDS[@]}"; do
        run_pairs $PROD_BOX $PROD_NSTEPS $seed "${GRAVITY_MODELS[@]}"
        if [ $((i % 2)) -eq 0 ]; then
            run_pairs 2048 $PROD_NSTEPS $seed "${GRAVITY_MODELS[@]}"
        else
            run_pairs 512 $PROD_NSTEPS $seed "${GRAVITY_MODELS[@]}"
        fi
        i=$((i + 1))
    done
}
campaign_steps() {
    for nsteps in 20 30 40; do run_pairs $PROD_BOX $nsteps ${SEEDS[0]} "${GRAVITY_MODELS[@]}"; done
}

CAMPAIGNS=("$@")
[ ${#CAMPAIGNS[@]} -eq 0 ] && CAMPAIGNS=(ensemble steps)
for CAMPAIGN in "${CAMPAIGNS[@]}"; do
    if ! declare -f "campaign_${CAMPAIGN}" > /dev/null; then
        echo "Unknown campaign '$CAMPAIGN' (choose from: ensemble steps)"; exit 1
    fi
    echo "=================================================="
    echo "$(date '+%F %T') Campaign $CAMPAIGN"
    "campaign_${CAMPAIGN}"
done

echo "=================================================="
echo "Done after $(( ($(date +%s) - START_TIME) / 3600 )) h: $nrun run, $nskip already finished, $nfailed failed. Results in $RESULTS_DIR"
