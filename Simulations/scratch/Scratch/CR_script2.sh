#!/bin/bash

# Get today's date in YYYY-MM-DD format
today=$(date +"%Y-%m-%d")

# Define the directory path where the flag file should be located
flags_dir="./CR/output/flags"

# Maximum time to wait for completion signal (in seconds)
max_wait_time=$((10 * 24 * 60 * 60))  # 5 days in seconds

# Start time
start_time=$(date +%s)

# Check if the completion flag file exists for today's date
while [ ! -f "$flags_dir/completion_flag_${today}" ]; do
    current_time=$(date +%s)
    elapsed_time=$((current_time - start_time))

    # Check if timeout threshold exceeded
    if [ $elapsed_time -ge $max_wait_time ]; then
        echo "Timeout threshold exceeded. Exiting."
        exit 1
    fi

    sleep 1
done

# Once signal is received, proceed with plotting
sbatch -p general -N 1 --mem 7120 -n 1 -t 6-9:00:00 --mail-type=fail --mail-user=cwzhou@email.unc.edu --wrap="Rscript CR02.Simulation_Summary.R > $flags_dir/plots_output_${today}.txt"

