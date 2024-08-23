from topone import TopONE
import argparse
import logging

"""
TODO: Create Choice for user to input
- Our HD matrices
    - xbc1
    - xbe
    - xbz
- Their HD matrices
    - recombinant lineage name
- Their own aligned sequences
    - recombinant lineage name
"""
parser = argparse.ArgumentParser()
parser.add_argument("--xbc1", action="store_true")
parser.add_argument("--xbe", action="store_true")
parser.add_argument("--xbz", action="store_true")
parser.add_argument("--hdmat", action="store_true")

args = parser.parse_args()

"""
TODO: Choice only for using our hd matrices
"""
def get_filepaths(lineage):
    paths = []
    if lineage == "xbc.1":
        countries = ["cn", "ph", "sk", "sg", "us"]
    elif lineage == "xbe":
        countries = ["uk", "us"]
    else:
        countries = ["dk", "at", "ge"]

    for c in countries:
        paths.append(f"inputs/{lineage}/{c}_recom_aligned.csv")
        paths.append(f"inputs/{lineage}/{c}_nonrecom_aligned.csv")
        paths.append(f"inputs/{lineage}/{c}_mixed_aligned.csv")

    return paths


"""
TODO: filter down the choices according to user input
"""
def get_recom_dataframe(lineage):
    """
    Extract sequences from FASTA (r, nr, m)
    HDMat
    Fit Transform
    Append HDMat..?
    Append Homologies

    Create dataframe
    """
    topone = TopONE(samples=-1,
                    seqtype="VIR",
                    maxdim=2,
                    simulations=-1,
                    segsites=-1)

    # filepath = fname
    # seqs = topone.get_sample_sequences(fname)

    # hdmat = topone.get_hdmatrices(seqs)

    # hom = topone.fit_transform(hdmat)

    paths = get_filepaths(lineage)

    """
    TODO: HDMATRICES INPUT
    for each path
        get the hdmatrices
        apply fit_transform
        append homologies

    create the dataframe to outputs



    TODO: SEQUENCES INPUT
    see virseq.py
    """


def main():
    pass


if __name__  == "__main__":
    main()