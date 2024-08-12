"""
A Python script that provides the hamming distance matrices of a group of simulations.
"""

from scipy.spatial import distance
import re 
import numpy as np
import random as rand
from ripser import Rips, ripser

class TARG:
    def __init__(self, samples, simulations,
                 groups, segsites, maxdim) -> None:
        self.SAMPLES = samples
        self.SIMULATIONS = simulations
        self.GROUPS = groups
        self.SEGSITES = segsites
        self.MAXDIM = maxdim


    """
    Collect the sequences from one group of simulations.

    One text file consists of a group of simulations, with
    a defined mutation rate (theta) and recombinationr rate (rho).

    The function returns:
        simulations - a 2-dimensional list of dimension (SIMULATION, SAMPLES)
                      containing all simulations each with a list of
                      samples in biallelic (0 and 1) format.
    """
    def get_sample_sequences(self, fname):       
        sim = []

        fp = open(fname, "r")
        simulations = []
        for ln in fp:
            if re.match("^[01]+$", ln):
                sim.append(ln.replace("\n", ""))
                if len(sim) == self.SAMPLES:
                    simulations.append(sim)
                    sim = []
        fp.close()

        return simulations


    """
    Collect the Hamming distance matrices across all simulations.

    Parameters:
        simulations - a 2-dimensional list of dimension (SIMULATION, SAMPLES)
                      containing all simulations each with a list of
                      samples in biallelic (0 and 1) format.

    The function returns:
        hdm_sims - a 3-dimensional list of dimension (SIMULATION, SAMPLES, SAMPLES)
                   containing all simulations each with a 2D Hamming Distance matrix
                   of dimension (SAMPLES, SAMPLES)
    """
    def get_hdmatrices(self, simulations):
        hdm_sims = []
        for sim in simulations:
            size = len(sim)

            # Construct the Hamming Distance Matrix
            mat = np.zeros((size, size))
            for i in range(0, size):
                for j in range(i, size):
                    mat[i, j] = round(distance.hamming(
                            np.array(list(sim[i]), dtype=int),
                            np.array(list(sim[j]), dtype=int)
                        ) * self.SEGSITES)
                    mat[j, i] = mat[i, j]
            hdm_sims.append(mat)

        return hdm_sims
    

    """ 
    To model 'sparsity', we simply take a sample (select number of lines)
    from the samples of the simulation (we assume it comprises the population).
    Replicate sparsity in the population (i.e. simulation block).

    FORMAT:
    sequences[group or file][simulation][chromosome or line seq or sample]
    
    NOTE: Because we assume sparsity to be from an accurate population,
    we only need to consider one simulation.

    PSUEDOCODE:
    1. for each simulation/sample
        - for each entry in sparsity_list
            - take a random sample of the percentage on the list
            - append this random sample to the mother list of sparse rand samples

    Parameters:
        group - a 2-dimensional list containing __ simulations,
                each simulation with a list of __ samples
                in biallelic (0 and 1) format.

    The function returns:
        sparse_samples - a 2-dimensional list of dimension (10, SAMPLESIZE) 
                        containing the 10 sampled lists according to increments 
                        of 10%.

    """
    def get_sparse_samples(self, group):
        sparsity_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
        # (SIMULATIONS, SPARSITY, SAMPLESIZE)
        # sparse_list = []

        # for i in range(self.SIMULATIONS):

        # # (SPARSITY, SAMPLESIZE)
        # sparse_samples = []

        # for p in sparsity_list:
        #     sample_size = int(self.SAMPLES * p)
        #     sampled_seq = rand.sample(group[0], sample_size)
        #     sparse_samples.append(sampled_seq)

        # # sparse_list.append(sparse_samples)


        sparsampled_list = [] #over all simulations
        for sim_idx in range(self.SIMULATIONS):

            sparse_samples = [] # samples over varying sparsities

            for percentage in sparsity_list:
                sample_size = int(self.SAMPLES * percentage)
                sampled_seq = rand.sample(group[sim_idx], sample_size)
                sparse_samples.append(sampled_seq)

            sparsampled_list.append(sparse_samples)

        return sparsampled_list
    
    """
    We will now add a noise matrix to the Hamming Distance Matrix to replicate stochasticity.
    Returns a 4-d list of dimension (VARLIST, SIMULATIONS, SAMPLES, SAMPLES)
    """
    def add_stochasticity(self, hdmatrices):
        # to iterate over increasing values of noise variance at increments of 0.1
        var_list = [round(x / 10, 2) for x in range(0,21)]

        rand.seed(123)

        # # 21 * GROUPS * SIMULATIONS * SAMPLES * SAMPLES
        # noise_list = []
        # for var in var_list:
        #     noise_list.append(np.random.normal(0, var, (self.SIMULATIONS, self.SAMPLES, self.SAMPLES)))

        # # 21 * GROUPS * SIMULATIONS * SAMPLES * SAMPLES
        # noisymat_list = []
        # for noise in noise_list:
        #     noisymat_list.append(hdmatrices + noise)

        noisymat_list = []
        for sim_hd in hdmatrices:
            noise_list = []
            for var in var_list:
                error_matrix = np.random.normal(0, var, (self.SAMPLES, self.SAMPLES)) # generates error matrix
                np.fill_diagonal(error_matrix, 0) # Sets diagonal entries to 0
                noise_list.append(error_matrix) # appends revised error mat

            noisy_mat_sim = []
            for noise in noise_list:
                noisy_mat_sim.append(sim_hd + noise)

            noisymat_list.append(noisy_mat_sim)

        return noisymat_list
    

    def get_betti_numbers(self, homologies):
        bettis = np.array([])
        for b in homologies:
            bettis = np.append(bettis, len(b))
        
        return bettis
    

    def get_barcode_lengths(self, homologies):
        barcode_lens = []

        for i in range(0, self.MAXDIM + 1):
            homologies[i] = homologies[i][np.argsort(homologies[i][:, 1])]
            finite_lengths = homologies[i][:, 1] - homologies[i][:, 0]
            finite_lengths = finite_lengths[np.isfinite(finite_lengths)]
            barcode_lens.append(finite_lengths)
            
        # Replace empty arrays with an array containing zero
        if len(finite_lengths) == 0:
            finite_lengths = np.array([0])

        return barcode_lens


    def get_barcode_length_statistics(self, barcode_lens):
        barcode_length_vars  = []
        barcode_length_means = []
        
        for l in barcode_lens:
            barcode_length_means.append(np.mean(l))
            barcode_length_vars.append(np.var(l))

        return barcode_length_means, barcode_length_vars
    


    def get_homologies(self, hdlist):
        rips = Rips(maxdim=self.MAXDIM, do_cocycles=True)

        homologies = []
        for sim in hdlist:
            hom_grps = []

            for j in range(len(hdlist[sim])):
                hom_grps.append(rips.fit_transform(sim[j], distance_matrix=True))

            homologies.append(hom_grps)

        return homologies