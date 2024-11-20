#!/bin/bash

# Get the number of CPU cores and calculate half of it
ncpu=$(nproc)
max_jobs=$((ncpu / 2))

# Create an array to hold the commands
commands=()

# Define the pt-hat values
pthat_values="20 50 100 150 200 300"

# Add the commands to the array
for pthat in $pthat_values
do
    commands+=("./test_embedding.py run3_tstruct.yaml local_run3_files.txt -e 10000 --pt-hat-min ${pthat} -o output_pthat${pthat}_lhc_run3.root --lhc-run 3")
    commands+=("./test_embedding.py run2_tstruct.yaml local_run2_files.txt -e 10000 --pt-hat-min ${pthat} -o output_pthat${pthat}_lhc_run2.root --lhc-run 2")
done

# Export the commands array
export commands

# Use gnuparallel to run the commands in parallel with a maximum number of jobs
parallel -j $max_jobs ::: "${commands[@]}"