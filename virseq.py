import numpy as np
import random as rand
from scipy.spatial import distance
from Bio import SeqIO

class VirSeq:
    def __init__(self) -> None:
        pass

    """
    Obtains the sequences from a FASTA file.
    --------------------
    Parameters:

    fname   :   String
        A string of the filename of the FASTA file.

        
    Returns:
    
    seqs    :   list
        A list of strings containing the biological sequences.
    """
    def get_fasta_samples(self, fname):
        fasta_sequences = SeqIO.parse(open(fname),'fasta')
        seqs = []
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            seqs.append(sequence)
        
        return seqs


    """
    Processes the raw sequences from a FASTA file, removing a column across
    all sequences if an unknown nucleotide is found in any of the sequences.
    --------------------
    Parameters:

    fname   :   String
        The filename of the FASTA file.

    params  :   dictionary
        The parameters containing attributes from the TopONE object.


    Returns:

    seqs    :   list
        A list containing the processed sequences.
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
    Creates the Hamming distance matrices from a list of
    equal-length sequences.
    --------------------
    Parameters:

    simulations :   list
        A 1-D list of dimension containing the processed
        equal-length biological sequences.

    params  :   dictionary
        The parameters containing attributes from the TopONE object.

                    
    Returns:
    
    mat :   list
        A 2-D matrix of dimension (SAMPLES, SAMPLES) where each element
        mat[i][j] represents the Hamming distance of sequences i and j.
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