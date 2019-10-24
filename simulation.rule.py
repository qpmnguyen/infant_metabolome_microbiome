


rule simluation: 
  input: 
    simulations = ["snakemake_output/simulation/assc_overall_datasets.rds", "snakemake_output/simulation/none_overall_datasets.rds"]

rule generate_positive_simulations:
  input:
    data = 
