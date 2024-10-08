from topone import TopONE
import argparse
import logging
import random
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--mats", help="Argument of the recombinant lineage name and folder to get the CSV matrices")
parser.add_argument("--seqs", help="Argument of the recombinant lineage name and folder to get the FASTA sequences")
parser.add_argument("--fname", help="Argument of the filename for the sequence counts")
parser.add_argument("--verbose", action='store_true', help="Increases Logging of Messages")

args = parser.parse_args()

if args.verbose:
    logging.basicConfig(
        level=logging.INFO
    )
else:
    logging.basicConfig(
        level=logging.CRITICAL + 1
    )

log = logging.getLogger("TOPOne")

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

    # Get the countries and country codes from the sequence counts
    log.info("[START] Create Country and Code Lists")
    df = pd.read_csv(fname)
    countries = df["country"].to_list()
    codes = df["country_code"].to_list()
    log.info("[END]   Create Country and Code Lists")
        
    # Get the filepaths to the recom, nonrecom, and mixed sequences
    log.info("[START] Get Sequence Filepaths")
    seq_paths = []
    for c in codes:
        c_path = []
        c_path.append(f"inputs/{lineage}/{c}_recom_aligned.fasta")
        c_path.append(f"inputs/{lineage}/{c}_nonrecom_aligned.fasta")
        c_path.append(f"inputs/{lineage}/{c}_mixed_aligned.fasta")
        seq_paths.append(c_path)
    log.info("[END]   Get Sequence Filepaths")

    # Instantiate the TopONE object
    topone = TopONE(samples=-1,
                    seqtype="VIR",
                    maxdim=2,
                    simulations=-1,
                    segsites=-1)

    log.info("[START] Create Hamming Distance Dataframes")
    for c in range(0, len(countries)):
        r = df.loc[c == df["country"]]["recom"]

        for p in seq_paths[c]:
            # Get the sequences from the FASTA file
            ss = topone.get_sample_sequences(p)

            # Recombinant Sequences
            if "_recom_" in p:
                hd = topone.get_hdmatrices(ss)
                np.savetxt(f"inputs/{lineage}/{codes[c]}_recom_aligned.csv", hd, delimiter=",")
           
            # Nonrecombinant Sequences
            elif "_nonrecom_" in p:
                nr_size = (300 - r) // 2

                if nr_size > len(ss):
                    nr_size = 999
                else:
                    ss = ss[:nr_size] + ss[-nr_size:]
                    
                hd = topone.get_hdmatrices(ss)
                np.savetxt(f"inputs/{lineage}/{codes[c]}_nonrecom_aligned.csv", hd, delimiter=",")
            
            # Mixed Sequences
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
    log.info("[END]   Create Hamming Distance Dataframes")

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

    log.info("[START] Create Country and Code Lists")
    df = pd.read_csv(fname)
    codes = df["country_code"].tolist()
    countries = df["country"].tolist()
    log.info("[END]   Create Country and Code Lists")

    log.info("[START] Get Sequence Filepaths")
    for c in codes:
        c_path = []
        c_path.append(f"inputs/{lineage}/{c}_recom_aligned.csv")
        c_path.append(f"inputs/{lineage}/{c}_nonrecom_aligned.csv")
        c_path.append(f"inputs/{lineage}/{c}_mixed_aligned.csv")
        paths.append(c_path)
    log.info("[END]   Get Sequence Filepaths")

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

    log.info("[START] Create Country and Code Lists")
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
    log.info("[END]   Create Country and Code Lists")

    log.info("[START] Get Sequence Filepaths")
    for c in codes:
        c_path = []
        c_path.append(f"inputs/{lineage}/{c}_recom_aligned.csv")
        c_path.append(f"inputs/{lineage}/{c}_nonrecom_aligned.csv")
        c_path.append(f"inputs/{lineage}/{c}_mixed_aligned.csv")
        paths.append(c_path)
    log.info("[END]   Get Sequence Filepaths")

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

    # Create List
    # (COUNTRIES, GENETYPES, MAXDIM+1, X)
    countries_homologies = []
    log.info("[START] Get All Homologies from All Countries")
    # COUNTRIES
    for i in range(0, len(paths)):
        c_homs = []
        
        # GENETYPES
        for j in range(0, len(gene_types)):
            hdmat = np.genfromtxt(paths[i][j], delimiter=',')
            
            # (MAXDIM+1, X)
            hom = topone.fit_transform(hdmat)
            c_homs.append(hom)
        countries_homologies.append(c_homs)
    log.info("[END]   Get All Homologies from All Countries")

    cols = np.array(["country", "country_code", "gene_type",
                     "homology", "b_time", "d_time"])
    df = pd.DataFrame(columns=cols)

    
    log.info("[START] Create Dataframe for all Homologies")
    # COUNTRIES
    for cn in range(0, len(countries)):
        # GENETYPES
        for gt in range(0, len(gene_types)):
            # MAXDIM+1 Homologies
            for h in range(0, topone.MAXDIM + 1):
                # X Pairs per Homology
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
    log.info("[END]   Create Dataframe for all Homologies")


def main():
    if args.mats != None and args.fname != None and args.seqs == None:
        get_recom_dataframe(args.mats, None, args.fname)
    elif args.mats == None and args.fname != None and args.seqs != None:
        t = get_custom_hdmatrices(args.seqs, args.fname)
        get_recom_dataframe(args.seqs, t, args.fname)
    else:
        if args.mats != None and args.seqs != None:
            print("Error: Multiple input types chosen. Choose one type only.")
        elif args.mats == None and args.seqs == None:
            print("Error: Please choose one input type.")
        
        if args.fname == None:
            print("Error: Please specify a filename for the sequence counts.")


if __name__  == "__main__":
    main()