from topone import TopONE
import argparse
import logging
import random
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--mats", help="Argument of the recombinant lineage name and folder")
parser.add_argument("--seqs", help="Custom Sequences")
parser.add_argument("--verbose", action='store_true', help="Increases Logging of Messages")

args = parser.parse_args()


def get_custom_hdmatrices(lineage, fname="sequence_counts.csv"):
    """
    Produces the Hamming distance matrices per country per gene type
    from the user's custom input, and are stored in the
    inputs/lineage/ directory.
    --------------------
    Parameters:

    lineage :   String
        A string of the recombinant lineage being studied.
    
    fname   :   String, default = "sequence_counts.csv"
        A string of the filename containing the number of sequences of the
        recombinant and parent lineages per country.
    """

    # Get the codes, countries
    df = pd.read_csv(fname)
    countries = df["country"].to_list()
    codes = df["country_code"].to_list()
        
    # Get the sequence filepaths
    seq_paths = []
    for c in codes:
        c_path = []
        c_path.append(f"inputs/{lineage}/{c}_recom_aligned.fasta")
        c_path.append(f"inputs/{lineage}/{c}_nonrecom_aligned.fasta")
        c_path.append(f"inputs/{lineage}/{c}_mixed_aligned.fasta")
        seq_paths.append(c_path)

    topone = TopONE(samples=-1,
                    seqtype="VIR",
                    maxdim=2,
                    simulations=-1,
                    segsites=-1)

    for c in range(0, len(countries)):
        r = df.loc[c == df["country"]]["recom"]
        p1 = df.loc[c == df["country"]]["parent_one"]
        p2 = df.loc[c == df["country"]]["parent_two"]

        for p in seq_paths[c]:
            ss = topone.get_sample_sequences(p)
            if "_recom_" in p:
                hd = topone.get_hdmatrices(ss)
                np.savetxt(f"inputs/{lineage}/{codes[c]}_recom_aligned.csv", hd, delimiter=",")
            elif "_nonrecom_" in p:
                # nr_size = (300 - r) // 2
                nr_size = 8

                ss = ss[:nr_size] + ss[-nr_size:]
                hd = topone.get_hdmatrices(ss)
                np.savetxt(f"inputs/{lineage}/{codes[c]}_nonrecom_aligned.csv", hd, delimiter=",")
            else:
                x = max(df.iloc[c]["recom"], nr_size)
                y = min(df.iloc[c]["recom"], nr_size)

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
                np.savetxt(f"inputs/{lineage}/{codes[c]}_mixed_aligned.csv", hd, delimiter=",")
    
    return topone


def get_custom_filepaths(lineage, fname="sequence_counts.csv"):
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
    codes = df["country_code"].tolist()
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


def get_recom_dataframe(lineage, topone, fname):
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

    if lineage in ["xbc.1", "xbe", "xbz"]:
        paths, codes, countries = get_sample_filepaths(lineage)
    else:
        paths, codes, countries = get_custom_filepaths(lineage, fname)

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
    if args.mats != None and args.seqs == None:
        get_recom_dataframe(args.mats, None, "sequence_counts_sample.csv")
    elif args.mats == None and args.seqs != None:
        t = get_custom_hdmatrices(args.seqs, "sequence_counts_sample.csv")
        get_recom_dataframe(args.seqs, t, "sequence_counts_sample.csv")
    else:
        print("Error: Multiple options chosen. Choose one option only.")


if __name__  == "__main__":
    main()