import numpy as np
import random as rand
from scipy.spatial import distance
from Bio import SeqIO

class VirSeq:
    def __init__(self) -> None:
        pass

    """
    Take a FASTA filename, and produce the list of sequences
    """
    def get_fasta_samples(self, fname):
        fasta_sequences = SeqIO.parse(open(fname),'fasta')
        seqs = []
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            seqs.append(sequence)
        
        return seqs


    """
    Processes the raw files from FASTA
    Returns a list of sequences
    """
    def get_sample_sequences(self, fname, params):
        ambuigities = ['R', 'Y', 'S', 'W', 'K', 'M',
                       'B', 'D', 'H', 'V', 'N']
        
        seqs = self.get_fasta_samples(fname)
        max_len = len(seqs[0])
        j = 0
        while j < max_len:
            i = 0
            deleted = False

            while i < len(seqs) and not(deleted):
                # If some letter in a sample is in the ambiguities
                if seqs[i][j] in ambuigities:
                    # Across each sample, delete the column
                    for k in range(0, len(seqs)):
                        seqs[k] = seqs[k][:j] + seqs[k][(j+1):]
                    deleted = True
                    i = 0
                    max_len -= 1
                    j -= 1

                i += 1
            j +=1

        return seqs
    

    """
    Produces the hamming distance matrices of a list of
    viral sample sequences.

    Parameters:
        simulations - a 1-D list of dimension (SAMPLES) containing the
                        processed biological sequences from viral samples.

    The function returns:
        mat         - a 2-D matrix of dimension (SAMPLES, SAMPLES).
    """
    def get_hdmatrices(self, sequences, params):
        # Construct the Hamming Distance Matrix
        strlen = len(sequences[0])
        samples = len(sequences)

        mat = np.zeros((samples, samples))
        for i in range(0, samples):
            for j in range(i, samples):
                mat[i, j] = round(distance.hamming(
                        list(sequences[i]),
                        list(sequences[j])
                    ) * strlen)
                mat[j, i] = mat[i, j]

        return mat
    
    
    """ 
    Takes a fraction of sequences from the population of sequences.

    Parameters:
        sequences - a 1-D list containing the list of viral sequences.

    Returns:
        sparse_samples - a 2-D list of dimension (10, SAMPLES),
                            where each element contains a fraction of
                            the population of sequences in varying
                            values of sparsity.
    """
    def sparsity_sampling(self, group, params):
        # (SPARSITY, SAMPLESIZES)
        sparse_samples = []

        for p in params["SPALIST"]:
            sample_size = int(params["SAMPLES"] * p)
            sampled_seq = rand.sample(group[0], sample_size)
            sparse_samples.append(sampled_seq)

        return sparse_samples
    

    """
    Adds varying stochasticity or noise to the Hamming Distance (HD) matrices
    present in every simulation.

    Parameters:
        hdmatrices - a 2-D HD matrix of dimension (SAMPLES, SAMPLES)

    Returns:
        returned_data - a 3-D list of dimension (21, SAMPLES, SAMPLES),
                        where 21 HD matrices are made for every variance value.
    """
    def add_stochasticity(self, hdmatrices, params):
        # Prepare the stochasticity matrix
        noise_list = []
        for var in params["VARLIST"]:
            noise_list.append(np.random.normal(0,
                                               var,
                                               (params["SAMPLES"],
                                                params["SAMPLES"])))

        # Add the stochasticity matrix to the HD matrices
        noisymat_list = []
        for noise in noise_list:
            # Ensure all diagonals are 0
            noise_mat = hdmatrices + noise
            idx = np.diag_indices_from(noise_mat)
            noise_mat[idx] = 0
            noisymat_list.append(noise_mat)

        return noisymat_list