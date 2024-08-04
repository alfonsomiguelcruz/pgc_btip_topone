A folder for all the simulation outputs we have.

Current outputs are generated using the `ms` package.

Folders are arranged according to date of simulation in the format `Sims_mmddyy`


## About the `ms` package

Command line: `./ms nsam nreps -t $\theta$ -s segsites -r $\rho$ nsites`

Definition of terms:
- `nsam` is the number of chromosomes in each sample (i.e. no. of lines)
- `nreps` is the number of replicate samples to generate (i.e. no. of groups/blocks)
- $\theta$ is the mutation parameter where $\theta = 4N_0 \mu$
  - $N_0$ is the effective diploid population size
  - $\mu$ is the neutral mutation rate for the entire segment
- `segsites` is the (fixed) number of segregating sites (i.e. length of each line)
- $\rho$ is the recombination parameter where $\rho = 4N_0 r$
  - $r$ is the probability of cross-over (recombination) per generation between the ends of the locus being simulated
- `nsites` is the number if base pairs in the locus

Sample usage:
To generate 100 samples of 5 chromosomes with $\theta = 3.0$:
Command: `./ms 5 100 -t 3.0 > outfile`

Output format:
- 2nd line is the number of polymorphic sites
- Segment of the chromosome being simulated is represented by the interval (0,1)
  - Where polymorphic sites (where mutation has occurred somewhere in the genealogy of the sample) is identified
- in biallelic format
  - 0 if the allele is in ancestral state (no mutation)
  - 1 if the allele is in mutated state
