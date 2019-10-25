# TODO: Figure out how to use one liners to add two arguments in a list
TIME = ["6W", "12M"]
MET = ["tar", "untar"]
TAX = ["16S"]
METHODS = ["rf", "svm", "enet", "spls"]
EXT = ['rds', 'svg']
OUTCOME = ['assc', 'none']
SNR = [0.2, 0.4, 0.6]
PERC_ASSC = [0.05, 0.1, 0.15]
# Rule set to govern the creation of simulated data sets

rule simluation: 
  input: 
    simulations = expand("snakemake_output/simulation/{outcome}_{snr}_{perc_assc}_overall_datasets.rds", outcome = OUTCOME, 
    snr = SNR, perc_assc = PERC_ASSC)

rule generate_positive_simulations:
  input:
    dataset1 = "data/raw/16S_12M_tar.rds",
    dataset2 = "data/raw/16S_6W_tar.rds",
    script = "R/simulation.R"
  output:
    files = "snakemake_output/simulation/{outcome}_{snr}_{perc_assc}_overall_datasets.rds"
  shell:
    """Rscript {input.script} --calibrate_data_1 {input.dataset1} --calibrate_data_2 {input.dataset2} \
    --n_datasets 50 --outcome_type {wildcard.outcome} --snr {wildcard.snr} --perc_assc {wildcard.perc_assc} --sample_size 300 \
    --calibrate_outcome TRUE"""
