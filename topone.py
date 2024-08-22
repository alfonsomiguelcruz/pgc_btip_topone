import re 
import numpy as np
import pandas as pd
import random as rand
import matplotlib.pyplot as plt
from scipy.spatial import distance
from ripser import Rips
from simseq import SimSeq
from virseq import VirSeq
from plot_diagrams import plot_diagrams

class TopONE:
    def __init__(self,
                 samples,
                 seqtype="SIM",
                 maxdim=2,
                 simulations=None,
                 segsites=None,
                ):
        """
        Initializes an instance of the class
        --------------------
        Parameters:

        samples     : int
            Number of samples or sequences in the input file
        
        seqtype     : String, default = "SIM"
            Defines the type of sequence used in the class, which can be:
            SIM - use of simulated sequences in biallelic (0 and 1) format
            VIR - use of viral sequences in ATCG format

        maxdim      : int, default = 2
            Highest homology dimension to be used in the computation
            of homologies, default set at 2, or the H2 homology group

        simulations : int, default = None
            Number of simulations used in an input file.
            Applicable only for "SIM" seqtype.

        segsites    : int, default = None
            Number of segregating sites in an input file.
            TODO: change name of segsites to seq_length
        """
        # Assign class variables to instantiation parameters
        self.SAMPLES = samples
        self.SEQTYPE = seqtype
        self.MAXDIM = maxdim
        self.SIMULATIONS = simulations
        self.SEGSITES = segsites

        # Create varlist and spalist for adding noise,
        # and taking fractions of the population
        self.varlist = np.around(np.linspace(0, 100, 21),
                                 decimals=2).tolist()
        self.spalist = np.around(np.linspace(0.1, 1, 10),
                                 decimals=2).tolist()
        
        # Instantiate Ripser for persistent homology computations
        # using the maxdim parameter
        self.ripser = Rips(maxdim=maxdim)

        # If the given seqtype by the user is neither SIM nor VIR
        # raise an error
        if self.SEQTYPE not in ["SIM", "VIR"]:
            raise Exception("Seqtype mode does not exist. " +
                            "Please choose between \"SIM\" and \"VIR\".")
        
        # Instantiate the helper functions from SimSeq and VirSeq
        self.simseq = SimSeq()
        self.virseq = VirSeq()

        # Place all parameters in a dictionary `params`
        self.params = {
            "SAMPLES": self.SAMPLES,
            "SEQTYPE": self.SEQTYPE,
            "MAXDIM": self.MAXDIM,
            "SIMULATIONS": self.SIMULATIONS,
            "SEGSITES": self.SEGSITES,
            "VARLIST": self.varlist,
            "SPALIST": self.spalist
        }


    def fit_transform(self, data, is_dist_mat=True):
        """
        Produces the birth and death time pairs in each homology,
        sorting the pairs according to increasing death time.
        --------------------
        Parameters:

        data        : numpy array
            The data used to produce the birth and death time
            pairs, which can be a distance matrix or a pair of numbers

        is_dist_mat : boolean, default = True
            A boolean value that dictates whether the data passed
            is a distance matrix or not

            
        Returns:

        hom : list
            The birth and death time pairs of each homology
            from H_0 to H_maxdim.
        """
        hom = self.ripser.fit_transform(data,
                                        distance_matrix=is_dist_mat)
        hom[0] = np.delete(hom[0], -1, axis=0)
        
        for i in range(self.ripser.maxdim + 1):
            hom[i] = hom[i][np.argsort(hom[i][:, 1])]

        return hom


    def switch_seqtype(self):
        """
        Switches the sequence type from SIM to VIR or vice versa.
        """
        if self.SEQTYPE == "SIM":
            self.SEQTYPE = "VIR"
        else:
            self.SEQTYPE = "SIM"

        print(f"Switched SEQTYPE to {self.SEQTYPE}.")


    def get_sample_sequences(self, fname):
        """
        Gets the sample sequences from an input file.

        For seqtype SIM, the sequences in each simulation are obtained
        from a generated population, stored in a text file.

        For seqtype VIR, the sequences are stored in a FASTA file.
        --------------------
        Parameters:

        fname   : String
            The filename of the input file.

            
        Returns:
        
        sequences   : list
            The list of sequences taken from the input files.
        """
        sequences = None
        if self.SEQTYPE == "SIM":
            sequences = self.simseq.get_sample_sequences(fname, self.params)
        else:
            sequences = self.virseq.get_sample_sequences(fname, self.params)

        return sequences


    def get_hdmatrices(self, sequences):
        """
        TODO: Convert Format like the ones above
        Produces the hamming distance matrix of a collection of sequences.

        Parameters:
            sequences - SEQTYPE VIR:
                            a 1-D list of dimension (SAMPLES) containing the
                            processed biological sequences from viral samples.
                        
                        SEQTYPE SIM:
                            a 2-D list of dimension (SIMULATION, SAMPLES)
                            containing all simulations each with a list of
                            in biallelic (0 and 1) formatted sequences.

        Returns:
            hdmatrices - SEQTYPE VIR:
                                a 2-D matrix of dimension (SAMPLES, SAMPLES).
                        
                        SEQTYPE SIM:
                            a 3-D list of dimension (SIMULATION, SAMPLES, SAMPLES)
                            where each simulation has a 2-D matrix of dimension
                            (SAMPLES, SAMPLES).
        """
        hdmatrices = None

        if self.SEQTYPE == "SIM":
            hdmatrices = self.simseq.get_hdmatrices(sequences,
                                                    self.params)
        else:
            hdmatrices = self.virseq.get_hdmatrices(sequences,
                                                    self.params)

        return hdmatrices
    

    def get_all_homologies(self, hdmatrices):
        """
        TODO: Convert Format like the ones above
        Gets all the homologies per variance/sparsity per simulation
        (21/10, SIMULATIONS, MAXDIM+1, X)
        """
        all_homologies = []
        for sim in hdmatrices:
            hom_grps = []
            for i in range(len(sim)):
                hom_grps.append(self.fit_transform(sim[i], True))
            all_homologies.append(hom_grps)
        
        return all_homologies
    

    
    def get_all_topoquants(self, homologies):
        """
        TODO: Convert Format like the ones above
        Gets all the topological quantities per variance per simulation
        (21, SIMULATIONS, MAXDIM+1)
        """
        betti_nums = []
        betti_mean_lens = []
        betti_vars_lens = []

        # 21
        for i in range(len(homologies)):
            sim_betti_nums = []
            sim_betti_mean_lens = []
            sim_betti_vars_lens = []
            
            # SIMULATIONS
            for j in range(len(homologies[i])):
                hom = homologies[i][j]

                bn = self.get_betti_numbers(hom)
                bl = self.get_barcode_lengths(hom)
                bml, blv = self.get_barcode_length_statistics(bl)

                sim_betti_nums.append(bn)
                sim_betti_mean_lens.append(bml)
                sim_betti_vars_lens.append(blv)
            
            betti_nums.append(sim_betti_nums)
            betti_mean_lens.append(sim_betti_mean_lens)
            betti_vars_lens.append(sim_betti_vars_lens)

        return betti_nums, betti_mean_lens, betti_vars_lens
                

    
    def sparsity_sampling(self, sequences):
        """
        TODO: Convert Format like the ones above
        Takes a fraction of sequences from the population of sequences.

        Parameters:
            sequences - SEQTYPE VIR:
                            a 1-D list containing the list of viral sequences.
                        
                        SEQTYPE SIM:
                            a 2D list of dimension (SIMULATION, SAMPLES)
                            where each simulation contains the population of
                            biallelic sequences.

        Returns:
            sparse_samples - SEQTYPE VIR:
                                a 2-D list of dimension (10, SAMPLES),
                                where each element contains a fraction of
                                the population of sequences in varying
                                values of sparsity.
                        
                            SEQTYPE SIM:
                                a 3-D list of dimension (10, SIMULATIONS, SAMPLES),
                                where each element contains a fraction of
                                the populations in each simulation with
                                varying values of sparsity.
        """
        sparse_samples = None

        if self.SEQTYPE == "SIM":
            sparse_samples = self.simseq.sparsity_sampling(sequences,
                                                           self.params)
        else:
            sparse_samples = self.virseq.sparsity_sampling(sequences,
                                                           self.params)

        return sparse_samples
    

    
    def add_stochasticity(self, hdmatrices):
        """
        TODO: Convert Format like the ones above
        Adds the stochasticity or noise to the Hamming Distance (HD) matrices

        Parameters:
            hdmatrices - SEQTYPE VIR:
                            a 2-D HD matrix of dimension (SAMPLES, SAMPLES)
                        
                        SEQTYPE SIM:
                            a 3-D list of dimension (SIMULATION, SAMPLES, SAMPLES)
                            where each simulation has a 2-D matrix of dimension
                            (SAMPLES, SAMPLES).

        Returns:
            hdmatrices_stoch - SEQTYPE VIR:
                                    a 3-D list of dimension (21, SAMPLES, SAMPLES),
                                    where 21 HD matrices are made for every variance value.
                        
                            SEQTYPE SIM:
                                a 4-D list of dimension (21, SIMULATION, SAMPLES, SAMPLES)
                                where 21 HD matrices are made for every variance value,
                                done for each simulation.
        """
        hdmatrices_stoch = None

        if self.SEQTYPE == "SIM":
            hdmatrices_stoch = self.simseq.add_stochasticity(hdmatrices,
                                                          self.params)
        else:
            hdmatrices_stoch = self.virseq.add_stochasticity(hdmatrices,
                                                          self.params)
            
        return hdmatrices_stoch
    

    
    def get_betti_numbers(self, homologies):
        """
        TODO: Convert Format like the ones above
        Computes for the Betti numbers of every homology

        Parameters:
            homologies - a 3-D list, where each element is a homology
                        with a 2-D list of its birth and death time pairs

        Returns:
            betti_numbers - a 1-D list, where the kth element of the list
                            represents the kth Betti number of homology k
        """
        betti_numbers = np.array([])
        for b in homologies:
            betti_numbers = np.append(betti_numbers, len(b))
        
        return betti_numbers
    

    
    def get_barcode_lengths(self, homologies):
        """
        TODO: Convert Format like the ones above
        Computes for the barcode lengths of every homology.
        The barcode lengths are computed as the difference between
        the death time from the birth time.

        Parameters:
            homologies - a 3-D list, where each element is a homology
                        with a 2-D list of its birth and death time pairs

        Returns:
            barcode_lengths - a 2-D list, where each element is a homology
                            with a 1-D list of the barcode lengths
        """
        barcode_lengths = []

        for i in range(0, self.ripser.maxdim + 1):
            barcode_lengths.append(homologies[i][:, 1] - homologies[i][:, 0])

        return barcode_lengths


    
    def get_barcode_length_statistics(self, barcode_lengths):
        """
        TODO: Convert Format like the ones above
        Computes for the barcode lengths of every homology.
        The barcode lengths are computed as the difference between
        the death time from the birth time.

        Parameters:
            barcode_lengths - a 2-D list, where each element is a homology
                            with a 1-D list of the barcode lengths

        Returns:
            barcode_length_means - a 1-D list, where the kth element
                                of the list represents the mean of
                                the barcode lengths in homology k

            barcode_length_vars - a 1-D list, where the kth element
                                of the list represents the variance of
                                the barcode lengths in homology k
        """
        barcode_length_means = []
        barcode_length_vars  = []
        
        for l in barcode_lengths:
            barcode_length_means.append(np.mean(l))
            barcode_length_vars.append(np.var(l))

        return barcode_length_means, barcode_length_vars
    

    def plot(self, hom, title="Updated Persistence Diagram"):
        """
        TODO: Convert Format like the ones above
        """
        plot_diagrams.plot(hom, title)


    def plot_multiple_diagrams(self, varshom):
        """
        TODO: Convert Format like the ones above
        """
        fig, axs = plt.subplots(nrows=21, ncols=5, figsize=(25,105))

        # Loop over the range of variances and simulations
        for var in range(len(self.varlist)):
            for sim in range(self.params["SIMULATIONS"]):
                self.ripser.plot(varshom[var][sim], ax=axs[var, sim], show=False)

        plt.show()


    def plot_barcode_diagram(self, betti_numbers, hom, title="Barcode Diagram"):
        """
        TODO: Convert Format like the ones above
        """
        plot_diagrams.plot_barcode_diagrams(betti_numbers, hom, title)


    def plot_topoquant_diagram(self, betti_numbers,
                               barcode_mean_length,
                               barcode_vars_length,
                               save_as_pdf=False,
                               pdf_fname=None):
        """
        TODO: Convert Format like the ones above
        """
        plot_diagrams.plot_topoquant_diagram(self.params,
                                             betti_numbers,
                                             barcode_mean_length,
                                             barcode_vars_length,
                                             save_as_pdf,
                                             pdf_fname)
        
        
    def create_topquant_dataframe(self, betti_numbers,
                                  barcode_mean_length,
                                  barcode_vars_length,
                                  list_type,
                                  fname):
        """
        TODO: Convert Format like the ones above
        """
        val = len(self.params[list_type])
        
        if list_type == "SPALIST":
            data_type = "sparsity"
        else:
            data_type = "variance"

        cols = ["simulation",
                data_type,
                "h0",
                "h1",
                "h2",
                "avg_h0",
                "avg_h1",
                "avg_h2",
                "var_h0",
                "var_h1",
                "var_h2"]
        df = pd.DataFrame(columns=cols)

        # Per simulation, per variance value
        # 21 or 10
        for v in range(val):
            # SIMULATIONS
            for i in range(self.SIMULATIONS):
                if list_type == "SPALIST":
                    list_val = self.spalist[v]
                else:
                    list_val = self.varlist[v]

                row = pd.DataFrame([{"simulation": i,
                                     data_type: list_val,
                                     "h0": betti_numbers[v][i][0],
                                     "h1": betti_numbers[v][i][1],
                                     "h2": betti_numbers[v][i][2],
                                     "avg_h0": barcode_mean_length[v][i][0],
                                     "avg_h1": barcode_mean_length[v][i][1],
                                     "avg_h2": barcode_mean_length[v][i][2],
                                     "var_h0": barcode_vars_length[v][i][0],
                                     "var_h1": barcode_vars_length[v][i][1],
                                     "var_h2": barcode_vars_length[v][i][2],
                                    }])
                df = pd.concat([df, row], ignore_index = True)

        df.sort_values(by=["simulation", data_type], ascending=[True, True])
        df.to_csv(fname, index=False)