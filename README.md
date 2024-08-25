# TopONE (Topological Omics for Non-tree Evolution)

This project aims to understand how recombination events influence various topological quantities of a genomic dataset using topological data analysis (TDA). The project seeks to determine whether topology can detect recombination events. More specifically, this project asks the following questions :

1. Is topology robust to noise and sparse genomic samples?
2. Does recombination change the topology of genomic samples?

The project's repository makes use of scripts and sequence files to answer the questions above. This project was made during the 2024 Bioinformatics Training and Internship Program headed by the Philippine Genome Center.

## Installation
### Population Genetics Software `ms`
Text Here

### Sequence Alignment Software `nextclade`
Text Here

### Required Libraries
Text Here

Command for Execution: `pip install -r requirements.txt`


## Pipelines
### Goal 01: Is topology robust to noise and sparse genomic samples?
<!-- ![goal01_pipeline](./imgs/goal01_pipeline.png) -->
<a href="url"><img src="https://github.com/alfonsomiguelcruz/pgc_btip_topone/blob/main/imgs/goal01_pipeline.png" height="48" width="48" ></a>

Text Here

### Goal 02: Does recombination change the topology of genomic samples?
![goal01_pipeline](./imgs/goal02_pipeline.png)
## Data
### ms
The simulated samples were generated from `ms`. These sequences are in a biallelic format, where the characters in a sequence are `1` or `0` only.

### GISAID
The virus sequence samples, specifically SARS-CoV-2 samples, were taken from the GISAID database. Three recombinant lineages and their parent lineages were taken across different countries. The table below shows the lineages used and the countries whose samples were used:

| Recombinant Lineage | Parent Lineage 1 | Parent Lineage 2 | Countries
| --- | --- | --- | --- |
| XBC.1 | BA.2 | B.1.617.2 | China, Philippines, Singapore, South Korea, United States of America
| XBE | BA.5.2 | BE.4.1 | United Kingdom, United States of America
| XBZ | BA.5.2.1 | EF.1.3 | Austria, Denmark, Germany

## Execution
### Goal 01: Is topology robust to noise and sparse genomic samples?
The script `goal_one_plots.py` uses simulated sequence data to produce the plots necessary to answer this question.

To execute the script, simply execute `python goal_one_plots.py` on your command line with any of the following arguments:

- `--get-one-plot`: gets one PDF file containing two plots: the variance and sparsity plots. Using this argument must also require the use of three additional arguments:
    
    - `-n`: simulated sample size $n$
    - `-t`: mutation rate $\theta$
    - `-r`: recombination rate $\rho$

- `--get-all-plots`: gets all PDF files containing the variance and sparsity plots, across every combination of possible simulated sample sizes, mutation rates, and recombination rates

- `--verbose`: Optional argument that outputs the logs the start and end of each step in the script.


Some sample commands are shown below:

- `python goal_one_plots.py --get-one-plot -n 100 -t 50 -r 72 --verbose`
- `python goal_one_plots.py --get-all-plots`

### Goal 02: Does recombination change the topology of genomic samples?
The script `goal_two_plots.py` uses the GISAID sequence data to produce the data required for the plots and for the statistical analysis.

To execute the script, simply execute `python goal_two_plots.py` on your command line with any of the following arguments:

- `--input`: Uses the sample Hamming distance matrices provided in the repository. The user can choose between `xbc.1`, `xbe`, and `xbz`.
- `--mats` : Uses the Hamming distance matrices provided by the user themselves.
- `--seqs` : Uses the FASTA files containing the aligned sequences provided by the user themselves.


## Results
<!-- Show Results of the Project in Paragraph form only -->
