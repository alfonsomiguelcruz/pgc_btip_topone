from topone import TopONE
import argparse
import logging
import numpy as np
import pandas as pd

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
parser.add_argument("--input", help="Input Type")
parser.add_argument("--mats", help="Custom Matrices")
parser.add_argument("--seqs", help="Custom Sequences")

args = parser.parse_args()


def get_custom_filepaths(lineage, fname="countries.csv"):
    paths = []

    df = pd.read_csv(fname)
    codes = df["country_codes"].tolist()
    countries = df["country"].tolist()

    for c in codes:
        c_path = []
        c_path.append(f"inputs/{lineage}/{c}_recom_aligned.csv")
        c_path.append(f"inputs/{lineage}/{c}_nonrecom_aligned.csv")
        c_path.append(f"inputs/{lineage}/{c}_mixed_aligned.csv")
        paths.append(c_path)

    return paths, codes, countries


"""
TODO: Choice only for using our hd matrices
"""
def get_sample_filepaths(lineage):
    paths = []

    if lineage == "xbc.1":
        codes = ["cn", "ph", "sk", "sg", "us"]
        countries = ["China",
                     "Philippines",
                     "South Korea",
                     "Singapore",
                     "United States of America"]
    elif lineage == "xbe":
        codes = ["uk", "us"]
        countries = ["United Kingdom",
                     "United States of America"]
    else:
        codes = ["dk", "at", "ge"]
        countries = ["Denmark",
                     "Austria",
                     "Germany"]


    for c in codes:
        c_path = []
        c_path.append(f"inputs/{lineage}/{c}_recom_aligned.csv")
        c_path.append(f"inputs/{lineage}/{c}_nonrecom_aligned.csv")
        c_path.append(f"inputs/{lineage}/{c}_mixed_aligned.csv")
        paths.append(c_path)

    return paths, codes, countries


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
    gene_types = ["recom", "nonrecom", "mixed"]

    topone = TopONE(samples=-1,
                    seqtype="VIR",
                    maxdim=2,
                    simulations=-1,
                    segsites=-1)

    # filepath = fname
    # seqs = topone.get_sample_sequences(fname)

    # hdmat = topone.get_hdmatrices(seqs)

    # hom = topone.fit_transform(hdmat)

    paths, codes, countries = get_sample_filepaths(lineage)
    print(paths)
    # (COUNTRIES, GENETYPES, HOMOLOGIES, PAIRSPERHOM)
    countries_homologies = []

    # Create List
    for i in range(0, len(paths)):
        c_homs = []
        for j in range(0, 3):
            hdmat = np.genfromtxt(paths[i][j], delimiter=',')
            hom = topone.fit_transform(hdmat)
            c_homs.append(hom)
        countries_homologies.append(c_homs)


    cols = np.array(["country", "country_code", "gene_type",
                     "homology", "b_time", "d_time"])
    df = pd.DataFrame(columns=cols)

    
    # 5 countries
    for cn in range(0, len(countries)):
        # 3 gene types
        for gt in range(0, 3):
            # 3 homologies
            for h in range(0, topone.MAXDIM + 1):
                # pairs per homologies
                for i in range(0, len(countries_homologies[cn][gt][h])):
                    row = pd.DataFrame([{ "country": countries[cn],
                                          "country_code": codes[cn],
                                          "homology": h,
                                          "gene_type": gene_types[gt],
                                          "b_time": countries_homologies[cn][gt][h][i][0],
                                          "d_time": countries_homologies[cn][gt][h][i][1]
                                        }])
                    df = pd.concat([df, row], ignore_index=True)
    df.to_csv(f"outputs/country_homologies_{lineage}.csv", index=False)

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
    if args.input != None and args.mats == None and args.seqs == None:
        if args.input == "xbc.1":
            get_recom_dataframe("xbc.1")
        elif args.input == "xbe":
            get_recom_dataframe("xbe")
        elif args.input == "xbz":
            get_recom_dataframe("xbz")
        else:
            print("Error: Unknown recombinant lineage. " +
                  "Please choose between the following: xbc.1, xbe, xbz.")
    elif args.input == None and args.mats != None and args.seqs == None:
        get_custom_filepaths("testlin")
    elif args.input == None and args.mats == None and args.seqs != None:
        pass
    else:
        print("Error: Multiple options chosen. Choose one option only.")


if __name__  == "__main__":
    main()