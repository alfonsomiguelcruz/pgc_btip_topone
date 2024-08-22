import re
import numpy as np
import random as rand
from scipy.spatial import distance

class SimSeq:
    def __init__(self) -> None:
        pass

    """
    Collect the sequences from one group of simulations.

    One text file consists of a group of simulations, with
    a defined mutation rate (theta) and recombination rate (rho).

    Parameters:
        fname   - a string containing the filename of the text file
        params  - a dictionary of additional parameters from the TopONE object

    Returns:
        simulations - a 2-dimensional list of dimension (SIMULATION, SAMPLES)
                      containing all simulations each with a list of
                      samples in biallelic (0 and 1) format.
    """
    def get_sample_sequences(self, fname, params):       
        sim = []

        fp = open(fname, "r")
        simulations = []
        for ln in fp:
            if re.match("^[01]+$", ln):
                sim.append(ln.replace("\n", ""))
                if len(sim) == params["SAMPLES"]:
                    simulations.append(sim)
                    sim = []
        fp.close()

        return simulations
    

    """
    Produces the hamming distance matrices of each simulation, where
    each simulation has a list of sequences in biallelic format.

    Parameters:
        simulations - a 2-D list of dimension (SIMULATION, SAMPLES)
                        containing all simulations each with a list of
                        in biallelic formatted sequences.

    Returns:
        hdm_sims - a 3-D list of dimension (SIMULATION, SAMPLES, SAMPLES)
                       where each simulation has a 2-D matrix of dimension
                       (SAMPLES, SAMPLES).
    """
    def get_hdmatrices(self, simulations, params):
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
                        ) * params["SEGSITES"])
                    mat[j, i] = mat[i, j]
            hdm_sims.append(mat)

        return hdm_sims
    

    """ 
    Takes a fraction of sequences from the population of sequences.

    Parameters:
        sequences - a 2D list of dimension (SIMULATION, SAMPLES)
                        where each simulation contains the population of
                        biallelic formatted sequences.

    Returns:
        sparse_samples - a 3-D list of dimension (10, SIMULATIONS, SAMPLES),
                            where each element contains a fraction of
                            the populations in each simulation with
                            varying values of sparsity.
    """
    def sparsity_sampling(self, group, params):
        # (10, SIMULATIONS, SAMPLES)
        sparsampled_list = [] #over all simulations

        for p in params["SPALIST"]:
            sparse_samples = []
            for s in range(params["SIMULATIONS"]):
                sample_size = int(params["SAMPLES"] * p)
                sampled_seq = rand.sample(group[s], sample_size)
                sparse_samples.append(sampled_seq)
            sparsampled_list.append(sparse_samples)
            

        return sparsampled_list
    

    """
    Adds varying stochasticity or noise to the Hamming Distance (HD) matrices
    present in every simulation.

    Parameters:
        hdmatrices - a 3-D list of dimension (SIMULATION, SAMPLES, SAMPLES)
                        where each simulation has a 2-D matrix of dimension
                        (SAMPLES, SAMPLES).

    Returns:
        returned_data - a 4-D list of dimension (21, SIMULATION, SAMPLES, SAMPLES)
                            where 21 HD matrices are made for every variance value,
                            done for each simulation.
    """
    def add_stochasticity(self, hdmatrices, params):
        # Prepare the stochasticity matrix
        rand.seed(123)
        noise_list = []
        # for var in params["VARLIST"]:
        #     noise_list.append(np.random.normal(0,
        #                                        var,
        #                                        (params["SIMULATIONS"],
        #                                         params["SAMPLES"],
        #                                         params["SAMPLES"])))

        # Add the stochasticity matrix to the HD matrices
        # (21, SIMULATIONS, SAMPLES, SAMPLES)
        # Var 0.0 = Sim[0], Sim[1], ...
        # Var 1.0 = Sim[0], Sim[1], ...
        noisymat_list = []
        sims = len(hdmatrices)
        # (SIMULATION, SAMPLES, SAMPLES) +
        for noise in params["VARLIST"]:
            noise_list = []
            for sim_hd in hdmatrices:
                min_val = -np.min(sim_hd)

                error_matrix = np.random.normal(0,
                                                np.sqrt(noise),
                                                (params["SAMPLES"], params["SAMPLES"]))

                np.fill_diagonal(error_matrix, 0)
                error_matrix = np.maximum(error_matrix, min_val)
                noise_list.append(sim_hd + error_matrix)
            noisymat_list.append(noise_list)

        # for noise in noise_list:
        #     noisy_sim_list = []
        #     for i in range(0, sims):
        #     # Ensure all diagonals are 0
        #         mat = hdmatrices[i] + noise[i]
        #         idx = np.diag_indices_from(mat)
        #         mat[idx] = 0
        #         noisy_sim_list.append(mat)
        #     noisymat_list.append(noisy_sim_list)

        return noisymat_list