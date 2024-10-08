from topone import TopONE
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import logging

parser = argparse.ArgumentParser()
parser.add_argument("--get-one-plot", action='store_true', help="Flag to Get One Plot only")
parser.add_argument("-nreps", type=int, nargs=1, help="Argument for the Number of Simulations")
parser.add_argument("-nsite", type=int, nargs=1, help="Argument for the Number of Segregating Sites")
parser.add_argument("-nsam", type=int, nargs=1, help="Argument for the Number of Samples")
parser.add_argument("-t", type=int, nargs=1, help="Argument for the Mutation Rate (Theta)")
parser.add_argument("-r", type=int, nargs=1, help="Argument for the Recombination Rate (Rho)")
parser.add_argument("--get-all-plots", action='store_true', help="Flag to Get All Plots")
parser.add_argument("-ns", type=int, nargs="+", help="Argument for Multiple Number of Samples")
parser.add_argument("-ts", type=int, nargs="+", help="Argument for Multiple Mutation Rates (Theta)")
parser.add_argument("-rs", type=int, nargs="+", help="Argument for Multiple Recombination Rates (Rho)")
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

def get_single_plot(nreps, nsite, nsam, theta, rho):
    """
    Gets the topological quantity plots for a given group of simulations
    that shows the trends of the quantities when adding noise or increasing
    the sparsity in the population.

    Also gets the CSV files containing the topological quantities
    per simulation used in the plots.
    --------------------
    Parameters:

    nreps   : int
        The number of simulations in a file.

    nsite   : int
        The number of segregating sites in a sequence.

    nsam    : int
        The number of samples per simulation.
    
    theta   : int
        The mutation rate applied in the simulations.

    rho     : int
        The recombination rate applied in the simulations.
    """
    
    # Variables for the PDF file of the topological quantity plots
    pdf_fname = 'outputs/goal_one_plots/plots_'+ f"n{nsam}_" + f"t{theta}_" + f"r{rho}" +'.pdf'
    pdf = PdfPages(pdf_fname)

    # Parameters of the simulations
    params = {
                # Number of lines/chromosomes per simulation
                "SAMPLES": nsam,
                # Number of simulations in the file
                "SIMULATIONS": nreps,
                # Chromosomal length
                "SEGSITES": nsite,
                # Filename of the input text file
                "FNAME": 'inputs/simseq/sims_'+ f"n{nsam}_" + f"t{theta}_" + f"r{rho}" +'.txt'
            }


    # Filenames for the variance and sparsity CSV files with the topological quantities
    var_fname = 'outputs/goal_one_dfs/v_'+ f"n{nsam}_" + f"t{theta}_" + f"r{rho}" + ".csv"
    samp_fname = 'outputs/goal_one_dfs/s_'+ f"n{nsam}_" + f"t{theta}_" + f"r{rho}" + ".csv"

    # Instantiate the TopONE object
    topone = TopONE(samples=params["SAMPLES"],
                    seqtype="SIM",
                    maxdim=2,
                    simulations=params["SIMULATIONS"],
                    segsites=params["SEGSITES"])
    
    # Get the sequences for each simulation from the input text file
    # dim: (SIMULATIONS, SAMPLES)
    log.info("[START] Get Sample Sequences")
    seqs = topone.get_sample_sequences(params["FNAME"])
    log.info("[END]   Get Sample Sequences")


    # Compute the Hamming distance matrices per simulation
    # dim: (SIMULATIONS, SAMPLES, SAMPLES)
    log.info("[START] Get Hamming Distance Matrices")
    hdmatrices = topone.get_hdmatrices(seqs)
    log.info("[END]   Get Hamming Distance Matrices")

    # Add stochasticity or noise to the Hamming distance matrices per simulation
    # dim: (21, SIMULATIONS, SAMPLES, SAMPLES)
    log.info("[START] Add Stochasticity to Hamming Distance Matrices")
    stochmatrices = topone.add_stochasticity(hdmatrices)
    log.info("[END]   Add Stochasticity to Hamming Distance Matrices")

    # Take fractions of the population per simulation
    # dim: (10, SIMULATIONS, SAMPLES)
    log.info("[START] Apply Sparsity Sampling to Sequences")
    sparse_samples = topone.sparsity_sampling(seqs)
    log.info("[END]   Apply Sparsity Sampling to Sequences")
    
    # Get the Hamming distance matrices from sparse populations per simulation
    # dim: (10, SIMULATIONS, SAMPLES, SAMPLES)
    log.info("[START] Get Hamming Distance Matrices from Sparse Populations")
    sparse_hdmats = []
    for s in range(len(topone.spalist)):
        sequences = sparse_samples[s]
        hdmats = topone.get_hdmatrices(sequences)
        sparse_hdmats.append(hdmats)
    log.info("[END]   Get Hamming Distance Matrices from Sparse Populations")


    # Get the homologies of noisy Hamming distance matrices
    # dim: (21, SIMULATIONS, MAXDIM+1, X)
    log.info("[START] Get All Homologies from Noisy Matrices")
    varshom = topone.get_all_homologies(stochmatrices)
    log.info("[END]   Get All Homologies from Noisy Matrices")

    # Get the topological quantities from noisy Hamming distance matrices
    # dim: (21, SIMULATIONS, MAXDIM+1) for all 3 quantities
    log.info("[START] Get All Topological Quantities from Noisy Matrices")
    bettis_v, mean_barcode_lengths_v, var_barcode_lengths_v = topone.get_all_topoquants(varshom)
    log.info("[END]   Get All Topological Quantities from Noisy Matrices")

    # Create the topological quantities dataframe from noisy Hamming distance matrices
    # returns a CSV file with all topological quantities per simulation and variance
    log.info("[START] Create the Topological Quantity Dataframe from Noisy Matrices")
    topone.create_topquant_dataframe(bettis_v, mean_barcode_lengths_v, var_barcode_lengths_v, "VARLIST", var_fname)
    log.info("[END]   Create the Topological Quantity Dataframe from Noisy Matrices")
    

    # Get the homologies of Hamming distance matrices from sparse populations
    # dim: (10, SIMULATIONS, MAXDIM+1, X)
    log.info("[START] Get All Homologies from Sparse Populations")
    sparsehom = topone.get_all_homologies(sparse_hdmats)
    log.info("[END]   Get All Homologies from Sparse Populations")
    
    # Get the topological quantities from Hamming distance matrices from sparse populations
    # dim: (10, SIMULATIONS, MAXDIM+1) for all 3 quantities
    log.info("[START] Get All Topological Quantities from Sparse Populations")
    bettis_s, mean_barcode_lengths_s, var_barcode_lengths_s = topone.get_all_topoquants(sparsehom)
    log.info("[END]   Get All Topological Quantities from Sparse Populations")

    # Create the topological quantities dataframe from Hamming distance matrices from sparse populations
    log.info("[START] Create the Topological Quantity Dataframe from Sparse Populations")
    topone.create_topquant_dataframe(bettis_s, mean_barcode_lengths_s, var_barcode_lengths_s, "SPALIST", samp_fname)
    log.info("[END]   Create the Topological Quantity Dataframe from Sparse Populations")


    # Create the topological quantity plots from Hamming distance matrices
    # that are noisy or matrices from sparse populations, enclosed in one PDF file
    log.info("[START] Create the Topological Quantity Plots from Noisy Matrices")
    topone.plot_topoquant_diagram(bettis_v, mean_barcode_lengths_v, var_barcode_lengths_v, "VARLIST", True, pdf)
    log.info("[END]   Create the Topological Quantity Plots from Noisy Matrices")
    log.info("[START] Create the Topological Quantity Plots from Sparse Populations")
    topone.plot_topoquant_diagram(bettis_s, mean_barcode_lengths_s, var_barcode_lengths_s, "SPALIST", True, pdf)
    pdf.close()
    log.info("[END]   Create the Topological Quantity Plots from Sparse Populations")

    log.info(f"[COMPLETE] Topological Plots Complete for {pdf_fname} with nsam={nsam}, theta={theta}, rho={rho}")


def main():
    if args.get_one_plot and not(args.get_all_plots):
        if args.nreps != None and \
           args.nsite != None and \
           args.nsam  != None and \
           args.t     != None and \
           args.r     != None and \
           args.ns    == None and \
           args.ts    == None and \
           args.rs    == None:
            get_single_plot(args.nreps[0], args.nsite[0], args.nsam[0], args.t[0], args.r[0])
        else:
            if args.nreps == None:
                print("Error! Please indicate the argument for the number of simulations.")
            if args.nsite == None:
                print("Error! Please indicate the argument for the number of segregation sites.")
            if args.nsam == None:
                print("Error! Please indicate the argument for the number of samples.")
            if args.t == None:
                print("Error! Please indicate the argument for the mutation rate (theta).")
            if args.r == None:
                print("Error! Please indicate the argument for the recombination rate (rho).")
            if args.ns != None:
                print("Error! Multiple number of samples is available for --get-all-plots only.")
            if args.ts != None:
                print("Error! Multiple mutation rates is available for --get-all-plots only.")
            if args.rs != None:
                print("Error! Multiple recombination rates is available for --get-all-plots only.")
    elif not(args.get_one_plot) and args.get_all_plots:
        if args.ns    != None and \
           args.ts    != None and \
           args.rs    != None and \
           args.nreps != None and \
           args.nsite != None and \
           args.nsam  == None and \
           args.t     == None and \
           args.r     == None:
            for n in args.ns:
                for t in args.ts:
                    for r in args.rs:
                        get_single_plot(args.nreps[0], args.nsite[0], n, t, r)
        else:
            if args.nreps == None:
                print("Error! Please indicate the argument for the number of simulations.")
            if args.nsite == None:
                print("Error! Please indicate the argument for the number of segregation sites.")
            if args.nsam != None:
                print("Error! Single number of sample value is used only for --get-one-plot only.")
            if args.t != None:
                print("Error! Single mutation rate value is used only for --get-one-plot only.")
            if args.r != None:
                print("Error! Single recombination rate value is used only for --get-one-plot only.")
            if args.ns == None:
                print("Error! Please indicate the argument for the different number of samples.")
            if args.ts == None:
                print("Error! Please indicate the argument for the different mutation rates (theta).")
            if args.rs == None:
                print("Error! Please indicate the argument for the different recombination rates (rho).")
    elif args.get_one_plot == args.get_all_plots:
        print("Error: Both options chosen. Please choose one option only.")

if __name__== "__main__":
    main()