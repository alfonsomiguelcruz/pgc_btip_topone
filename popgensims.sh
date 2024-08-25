#!/bin/bash

# Parameter setting
nsam=(100 1000 10000)   # no. of samples and the hypothesize pop size (as assumed)
nreps=500               # no. of repetitions or simulations
nsites=2501             # loci or number of (recomb?) sites
segsites=300            # length per sample

# Define theta and rho values
theta=(50 500 5000) # 4 N_0 \mu
rho=(4 12 36 72 144) # 4 N_0 r

# Ensure Sims directory exists
mkdir -p final_sims

# Loop through each combination of theta and rho and nsam
for nsam in "${nsam[@]}"; do
    for theta in "${theta[@]}"; do
        for rho in "${rho[@]}"; do
            # Define output file name
            output_file="final_sims/sims_n${nsam}_t${theta}_r${rho}.txt"

            # Run the command with theta and rho values
            ./ms $nsam $nreps -t $theta -s $segsites -r $rho $nsites > $output_file
        done
    done
done

# make script executable: chmod +x popgensims.sh
# run script: ./popgensims.sh
# make sure script file is in the msdir path