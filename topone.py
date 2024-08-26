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
        Produces the hamming distance matrix of a collection of sequences.
        --------------------
        Parameters:

        sequences   :   list
            [SEQTYPE - VIR] a 1-D list of dimension (SAMPLES) containing
            the processed biological sequences from viral samples.
                        
            [SEQTYPE - SIM] a 2-D list of dimension (SIMULATION, SAMPLES)
            containing all simulations each with a list of in biallelic (0 and 1)
            formatted sequences.

            
        Returns:

        hdmatrices  :   list
            [SEQTYPE - VIR] a 2-D matrix of dimension (SAMPLES, SAMPLES).
                        
            [SEQTYPE - SIM] a 3-D list of dimension (SIMULATION, SAMPLES, SAMPLES)
            where each simulation has a 2-D matrix of dimension (SAMPLES, SAMPLES).
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
        Used exclusively for SEQTYPE SIM.

        Gets the homologies for each Hamming distance matrix per simulation.
        Computing homologies per simulation is repeated 21 times for each variance value
        in "varlist" to add noise to the matrices, or 10 times for each sparsity value in
        "spalist" to apply to the population of sequences.
        --------------------
        Parameters:

        hdmatrices  :   list
            A 3-D list of dimension (SIMULATIONS, SAMPLES, SAMPLES), containing
            the Hamming distance matrices of per simulation. For noise-added matrices,
            the matrix dimensions are (SAMPLES, SAMPLES). For matrices from sparse populations,
            the dimensions are (SAMPLES*X, SAMPLES*X), where X is a value in "spalist".

        
        Returns:

        all_homologies  :   list
            A 4-D list of dimension (21/10, SIMULATIONS, MAXDIM+1, X), where
            for every value of variance (21) or sparsity (10), each simulation has a
            list of (MAXDIM+1) homologies, with each homology group having X pairs of
            birth and death times.
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
        Used exclusively for SEQTYPE SIM.

        Gets the topological quantities (Betti numbers, barcode length mean, barcode length variance)
        in each list of homologies per simulation. Computing topological quantities is
        repeated 21 times for every variance value in "varlist" or 10 times for every
        sparsity value in "spalist".
        --------------------
        Parameters:

        homologies  :   list
            A 4-D list of dimension (21/10, SIMULATIONS, MAXDIM+1, X), where
            for every value of variance (21) or sparsity (10), each simulation has a
            list of (MAXDIM+1) homologies, with each homology group having X pairs of
            birth and death times.
        
        
        Returns:

        betti_nums      :   list
            A 3-D list of dimension (21/10, SIMULATIONS, MAXDIM+1), where for every
            value of variance (21) or sparsity (10), each simulation has a
            list of length (MAXDIM+1) containing the Betti numbers of each homology group.

        bcode_mean_lens :   list
            A 3-D list of dimension (21/10, SIMULATIONS, MAXDIM+1), where for every
            value of variance (21) or sparsity (10), each simulation has a list of length
            (MAXDIM+1) containing the barcode length mean for each homology group.

        bcode_vars_lens :   list
            A 3-D list of dimension (21/10, SIMULATIONS, MAXDIM+1), where for every
            value of variance (21) or sparsity (10), each simulation has a list of length
            (MAXDIM+1) containing the barcode length variance for each homology group.
        """

        betti_nums = []
        bcode_mean_lens = []
        bcode_vars_lens = []

        # 21/10
        for i in range(len(homologies)):
            sim_betti_nums = []
            sim_bcode_mean_lens = []
            sim_bcode_vars_lens = []
            
            # SIMULATIONS
            for j in range(len(homologies[i])):
                hom = homologies[i][j]

                # MAXDIM + 1
                bn = self.get_betti_numbers(hom)
                bl = self.get_barcode_lengths(hom)
                bml, blv = self.get_barcode_length_statistics(bl)

                sim_betti_nums.append(bn)
                sim_bcode_mean_lens.append(bml)
                sim_bcode_vars_lens.append(blv)
            
            betti_nums.append(sim_betti_nums)
            bcode_mean_lens.append(sim_bcode_mean_lens)
            bcode_vars_lens.append(sim_bcode_vars_lens)

        return betti_nums, bcode_mean_lens, bcode_vars_lens
                

    
    def sparsity_sampling(self, sequences):
        """
        Used exclusively for SEQTYPE SIM.

        Takes a fraction of sequences from the population of sequences, where each
        fraction or percentage to apply to the sequences is taken from "spalist", a
        list containing fractional values from 0.1 to 1.0 in 0.1 increments.
        --------------------
        Parameters:
        
        sequences   :   list
            A 2D list of dimension (SIMULATION, SAMPLES)
            where each simulation contains the population of biallelic sequences.

            
        Returns:
        
        sparse_samples  :   list
            A 3-D list of dimension (10, SIMULATIONS, SAMPLES),
            where each element contains a fraction of the populations in each
            simulation with varying values of sparsity.
        """

        sparse_samples = None

        if self.SEQTYPE == "SIM":
            sparse_samples = self.simseq.sparsity_sampling(sequences,
                                                           self.params)
        else:
            print("Error! This is used for simulated sequence types only.")

        return sparse_samples
    

    
    def add_stochasticity(self, hdmatrices):
        """
        Used exclusively for SEQTYPE SIM.

        Adds stochasticity or noise to the Hamming distance matrices. The noise is taken
        from a Normal distribution of mean = 0 and standard deviation of sqrt(var), where
        var is a value in the list "varlist" and is added to the original Hamming distance
        matrices.
        --------------------
        Parameters:

        hdmatrices  :   list
            A 3-D list of dimension (SIMULATION, SAMPLES, SAMPLES)
            where each simulation has a 2-D matrix of dimension (SAMPLES, SAMPLES).

            
        Returns:
        
        hdmatrices_stoch    :   list
            A 4-D list of dimension (21, SIMULATION, SAMPLES, SAMPLES)
            where 21 HD matrices are made for every variance value, for each simulation.
        """

        hdmatrices_stoch = None

        if self.SEQTYPE == "SIM":
            hdmatrices_stoch = self.simseq.add_stochasticity(hdmatrices,
                                                          self.params)
        else:
            print("Error! This is used for simulated sequence types only.")
            
        return hdmatrices_stoch
    

    
    def get_betti_numbers(self, homologies):
        """
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
        Gets the means and variances of the barcode lengths for every homology.
        --------------------
        Parameters:
        barcode_lengths :   list
            A 2-D list of dimension (MAXDIM+1, X), containing the
            X barcode lengths all (MAXDIM+1) homologies.

            
        Returns:
        barcode_length_means    :   list
            A 1-D list with (MAXDIM+1) elements, where the kth element
            is the barcode length mean of the kth homology.

        barcode_length_vars     :   list
            A 1-D list with (MAXDIM+1) elements, where the kth element
            is the variance of the barcode lengths in homology H_k.
        """

        barcode_length_means = []
        barcode_length_vars  = []
        
        for l in barcode_lengths:
            barcode_length_means.append(np.mean(l))
            barcode_length_vars.append(np.var(l))

        return barcode_length_means, barcode_length_vars
    

    def plot(self, hom, title="Persistence Diagram"):
        """
        Creates a persistence diagram plot, containing the birth and death time
        pairs of each homology group in the topology, and the frequency distributions
        of the birth and death times for each homology, except for the H0 birth times.
        --------------------
        Parameters:

        hom     :   list
            A 2-D list of dimension (MAXDIM+1, X), containing the (MAXDIM+1) homology
            groups, each with its list of birth and death time pairs.

        title   :   String, default = "Persistence Diagram"
            A string containing the title of the persistence diagram.
        """

        plot_diagrams.plot(hom, title)


    def plot_barcode_diagram(self, betti_nums, hom, title="Barcode Diagram"):
        """
        Creates a barcode diagram, that visualizes the lengths of the barcodes for
        every homology group, which is the difference of the death times from the birth times.
        --------------------
        Parameters:

        betti_nums  :   list
            A 1-D list of length (MAXDIM+1), where the kth element is the
            kth Betti number of the kth homology.

        hom         :   list
            A 2-D list of dimension (MAXDIM+1, X), containing the (MAXDIM+1) homology
            groups, each with its list of birth and death time pairs.

        title       :   String, default = "Barcode Diagram"
            A string containing the title of the barcode diagram.
        """
        plot_diagrams.plot_barcode_diagrams(betti_nums, hom, title)


    def plot_topoquant_diagram(self, betti_numbers,
                               barcode_mean_length,
                               barcode_vars_length,
                               list_type,
                               save_as_pdf=False,
                               pdf_fname=None):
        """
        Creates the plots that displays the values of the topological quantities
        when variance increases or sparsity increases.
        --------------------
        Parameters:

        betti_nums          :   list
            A 3-D list of dimension (21/10, SIMULATIONS, MAXDIM+1), where for every
            value of variance (21) or sparsity (10), each simulation has a
            list of length (MAXDIM+1) containing the Betti numbers of each homology group.

        barcode_mean_length :   list
            A 3-D list of dimension (21/10, SIMULATIONS, MAXDIM+1), where for every
            value of variance (21) or sparsity (10), each simulation has a list of length
            (MAXDIM+1) containing the barcode length mean for each homology group.

        barcode_vars_length :   list
            A 3-D list of dimension (21/10, SIMULATIONS, MAXDIM+1), where for every
            value of variance (21) or sparsity (10), each simulation has a list of length
            (MAXDIM+1) containing the barcode length variance for each homology group.

        list_type           :   String
            A string that determines whether the dataframe will use values from
            the varlist ("VARLIST") or spalist ("SPALIST").

        save_as_pdf         :   Boolean, default = False
            A Boolean value dictating whether the plot should be saved in
            a PDF file or be displayed as is.

        pdf_fname           :   String, default = None
            A String containing the filename of the plots in a PDF file format.
        """
        plot_diagrams.plot_topoquant_diagram(self.params,
                                             betti_numbers,
                                             barcode_mean_length,
                                             barcode_vars_length,
                                             list_type,
                                             save_as_pdf,
                                             pdf_fname)
        
        
    def create_topquant_dataframe(self, betti_numbers,
                                  barcode_mean_length,
                                  barcode_vars_length,
                                  list_type,
                                  fname):
        """
        Creates the dataframe for the topological quantities per simulation
        and per value of sparsity or variance.
        --------------------
        Parameters:
        
        betti_nums          :   list
            A 3-D list of dimension (21/10, SIMULATIONS, MAXDIM+1), where for every
            value of variance (21) or sparsity (10), each simulation has a
            list of length (MAXDIM+1) containing the Betti numbers of each homology group.

        barcode_mean_length :   list
            A 3-D list of dimension (21/10, SIMULATIONS, MAXDIM+1), where for every
            value of variance (21) or sparsity (10), each simulation has a list of length
            (MAXDIM+1) containing the barcode length mean for each homology group.

        barcode_vars_length :   list
            A 3-D list of dimension (21/10, SIMULATIONS, MAXDIM+1), where for every
            value of variance (21) or sparsity (10), each simulation has a list of length
            (MAXDIM+1) containing the barcode length variance for each homology group
        
        list_type           :   String
            A string that determines whether the dataframe will use values from
            the varlist ("VARLIST") or spalist ("SPALIST").

        fname               :   String
            The filename of the dataframe.
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