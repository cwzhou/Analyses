#!/bin/bash

# ============================
# Local parallel execution with per-run logging (no date argument)
# ============================

# Script to run (CR01 or CR02)
SCRIPT=${1:-CR01.Simulation_Run.R}  # default to CR01

MAX_JOBS=4
job_count=0

# Function to run one simulation
run_sim() {
    l=$1
    k=$2
    j=$3
    i=$4
    h=$5
    g=$6
    
    # Create a unique logfile name based on parameters
    logfile="output_endpoint${1}_censor${g}_beta${h}_prop${i}_size${j}_crit${k}_cause${l}.txt"
    
    # Run Rscript and redirect stdout & stderr to logfile
    Rscript "$SCRIPT" 1 "$g" 1 "$h" "$i" "$j" "$k" "$k" "$l" > "$logfile" 2>&1
    
    # Optional: print progress
    echo "Finished simulation: l=$l k=$k j=$j i=$i h=$h g=$g"
}

# Nested loops (same as your cluster version)
for l in {1..2}
do
  for k in 1
  do
    for j in {1..2}
    do
      for i in {1..2}
      do
        for h in {1..2}
        do
          for g in {1..2}
          do
            # Start simulation in background
            run_sim "$l" "$k" "$j" "$i" "$h" "$g" &
            ((job_count++))

            # Limit number of parallel jobs
            if (( job_count >= MAX_JOBS )); then
              wait
              job_count=0
            fi
          done
        done
      done
    done
  done
done

# Wait for any remaining background jobs
wait

echo "All simulations completed!"
