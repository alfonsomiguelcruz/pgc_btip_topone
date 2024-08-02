#!/bin/bash

# Parameter setting
nsam=10 # no. of samples
nreps=5 # no. of repetitions or simulations
nsites=10000 # loci or number of (recomb?) sites
segsites=1000 #length per sample

# Define theta and rho values
theta=(50)
rho=(2 4 12 36 72)



# Ensure Sims directory exists
mkdir -p TestSims1

# Loop through each combination of theta and rho
for theta in "${theta[@]}"; do
    for rho in "${rho[@]}"; do
        # Define output file name
        output_file="TestSims1/sims_${theta}_${rho}.txt"

        # Run the command with theta and rho values
        ./ms $nsam $nreps -t $theta -s $segsites -r $rho $nsites > $output_file
    done
done

# make sciprt executable: chmod +x test_sims
# run script: ./test_sims.sh
# make sure script file is in the msdir path