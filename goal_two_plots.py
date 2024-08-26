from topone import TopONE
import argparse
import logging
import random
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="Input Type")
parser.add_argument("--mats", help="Custom Matrices")
parser.add_argument("--seqs", help="Custom Sequences")

args = parser.parse_args()


def get_custom_hdmatrices(lineage, fname="sequence_counts.csv"):
    """
    Gets the Hamming distance matrices per country per gene type
    from the user's custom input.
    --------------------
    Parameters:

    lineage :   String
        A string of the recombinant lineage being studied.
    
    fname   :   String, default = "sequence_counts.csv"
        A string of the filename containing the number of sequences of the
        recombinant and parent lineages per country.

    
    Returns:

    countries_matrices  :   list
        A 4-D list of dimension (NUM_COUNTRIES, NUM_GENE_TYPE, SAMPLES, SAMPLES), that
        shows the Hamming distance matrices per gene type for each country.
    """

    paths, codes, countries = get_custom_filepaths(lineage)
    df = pd.read_csv(fname)
    countries_matrices = []

    topone = TopONE(samples=-1,
                    seqtype="VIR",
                    maxdim=2,
                    simulations=-1,
                    segsites=-1)


    for c in range(0, len(countries)):
        r = df.loc[c == df["country"]]["recom"]
        p1 = df.loc[c == df["country"]]["parent_one"]
        p2 = df.loc[c == df["country"]]["parent_two"]

        c_hds = []
        # TODO: Currently path contains csv files for HDmats, not sequences; to fix for seqs.
        for p in paths[c]:
            ss = topone.get_sample_sequences(p)
            if "_recom_" in p:
                hd = topone.get_hdmatrices(ss)
                c_hds.append(hd)
            elif "_nonrecom_" in p:
                nr_size = (300 - r) // 2

                ss = ss[:nr_size] + ss[-nr_size:]
                hd = topone.get_hdmatrices(ss)
                c_hds.append(hd)
            else:
                x = max(p1, p2)
                y = min(p1, p2)

                if len(ss[:x]) <= nr_size:
                    nr1 = ss[:x]
                else:
                    nr1 = random.sample(ss[:x], nr_size)

                if len(ss[x:(x+y)]) <= nr_size:
                    nr2 = ss[x:(x+y)]
                else:
                    nr2 = random.sample(ss[x:(x+y)], nr_size)

                r = ss[(x+y):]
                ss = nr1 + nr2 + r

                hd = topone.get_hdmatrices(ss)
                c_hds.append(hd)

        
        countries_matrices.append(c_hds)

    return countries_matrices


def get_custom_filepaths(lineage, fname="countries.csv"):
    """
    Gets the path of the Hamming distance matrices in
    CSV file format provided by the user.
    --------------------
    Parameters:

    lineage :   String
        The name of the lineage being studied by the user.

    fname   :   String, default = "countries.csv"
        The filename of a CSV file containing the countries being studied and their
        2-digit codes.

    
    Returns:

    paths       :   list
        A list of paths of the Hamming distance matrices in CSV file format
        from recombinant, nonrecombinant, and mixed sequences.

    codes       :   list
        A list of the 2-character country codes of each country.

    countries   :   list
        A list of the names of each country being studied.
    """
    
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


def get_sample_filepaths(lineage):
    """
    Gets the filepaths used in the study, used as a sample to demonstrate
    the results and how the pipelines in the study work.
    --------------------
    Parameters:

    lineage :   String
        The recombinant lineage being studied.

    
    Returns:
    
    paths       :   list
        A list of paths of the Hamming distance matrices in CSV file format
        from recombinant, nonrecombinant, and mixed sequences.

    codes       :   list
        A list of the 2-character country codes of each country.

    countries   :   list
        A list of the names of each country being studied.
    """

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


def get_recom_dataframe(lineage, topone=None):
    """
    Constructs the dataframe containing the birth and death time pairs of each
    homology group per gene type per country.
    --------------------
    Parameters:

    lineage :   String
        A string of the recombinant lineage being studied.

    topone  :   TopONE, default = None
        An object instance of the TopONE class.
    """
    gene_types = ["recom", "nonrecom", "mixed"]

    if topone is None:
        topone = TopONE(samples=-1,
                        seqtype="VIR",
                        maxdim=2,
                        simulations=-1,
                        segsites=-1)

    paths, codes, countries = get_sample_filepaths(lineage)

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