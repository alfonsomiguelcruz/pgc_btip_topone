#!/bin/bash

# Parameter setting
nsam=100 # no. of samples
nreps=500 # no. of repetitions or simulations
nsites=10000 # loci or number of sites


# Define theta and rho values
theta=(50 500 5000)
rho=(4 12 36 72)



# Ensure Sims directory exists
mkdir -p Sims

# Loop through each combination of theta and rho
for theta in "${theta[@]}"; do
    for rho in "${rho[@]}"; do
        # Define output file name
        output_file="Sims/sim_${theta}_${rho}.txt"

        # Run the command with theta and rho values
        ./ms $nsam $nreps -t $theta -r $rho $nsites > $output_file
    done
done

# make sciprt executable: chmod +x "test sim.sh""
# run script: ./"test sim.sh"
# make sure script file is in the msdir path