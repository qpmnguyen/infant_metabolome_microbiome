TIME = ["6W", "12M"]
MET = ["tar", "untar"]
TAX = ["16S"]
METHODS = ["rf", "svm", "enet", "spls"]

rule all:
    input:
        expand('data/{tax}_{time}_{met}.rds', tax = TAX, time = TIME, met = MET)

rule data_retrieval:
    input:
        dir_file = "data/data_directory.csv",
        script = "R/data_retrieval.R",
        fasttree_dir = "../software/"
    output: 
        out_file = "data/{tax}_{time}_{met}.rds"
    shell:
        "RScript {input.script} --input {input.dir_file} --output {output.out_file} --time {wildcards.time} --metab_type {wildcards.met} --tax_type {wildcards.tax} --fasttree_dir {input.fasttree_dir}"  

#rule prediction_processing:
#    input:
#        "data/{tax}_{time}_{met}.rds"
#    output:
#        "processed/{tax}_{time}_{met}_processed.rds"
#    shell:
#        "echo {output}"

# rule exploratory_processing:
#    input: 
#        "data/{tax}_{time}_{met}.rds"
#    output:
#        "processed/{tax}_{time}_{met}_ordination_processed.rds"
#    shell:
#       "echo {output}"

#rule fitting_prediction_models: 
#    input:
#    output: 

#rule back_transformation:

#rule calculate_performance_metrics:

#rule generate_prediction_result_figures:

#rule generate_ordination_figures:

#rule generate_correlation_figures:

#rule generate_report:

#rule general_plotting:




