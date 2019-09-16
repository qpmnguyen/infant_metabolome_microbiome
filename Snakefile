TIME = ["6W", "12M"]
MET = ["tar", "untar"]
TAX = ["16S"]
METHODS = ["rf", "svm", "enet", "spls"]
EXT = ['rds','svg']

rule all:
    input:
        joint_ordinations = expand('snakemake_output/figures/ordinations/{tax}/{time}/{met}/joint_protest_ordinations.{ext}', 
        tax = TAX, time = TIME, met = MET, ext = EXT),
        met_ordinations = expand("snakemake_output/figures/ordinations/{tax}/{time}/{met}/met_ordinations.{ext}", 
        tax = TAX, time = TIME, met = MET, ext = EXT),
        tax_ordinations = expand("snakemake_output/figures/ordinations/{tax}/{time}/{met}/tax_ordinations.{ext}", 
        tax = TAX, time = TIME, met = MET, ext = EXT),
        eval_plots = expand("snakemake_output/figures/prediction/{tax}_{time}_{met}_{methods}_prediction_plots.{ext}", tax = TAX, 
        time = TIME, met = MET, methods = METHODS, ext = EXT),
        sparse_cca_plots = expand("snakemake_output/figures/correlation/{tax}_{time}_{met}_scca_plots.{ext}", 
        tax = TAX, time = TIME, met = MET, ext = EXT),
        spearman_correlation_plots = expand("snakemake_output/figures/correlation/{tax}_{time}_{met}_scc_plots.{ext}", 
        tax = TAX, time = TIME, met = MET, ext = EXT) 

rule clean_analyses:
    shell:
        """
        rm -rf figures/
        rm -rf analyses/
        """
rule clean_data:
    shell:
        """
        rm -rf data/
        rm -rf temp/
        """
rule data_retrieval:
    input:
        dir_file = "data/data_directory.csv",
        script = "R/data_retrieval.R",
        fasttree_dir = "../software/"
    output: 
        out_file = "data/raw/{tax}_{time}_{met}.rds"
    shell:
        "RScript {input.script} --input {input.dir_file} --output {output.out_file} --time {wildcards.time} --metab_type {wildcards.met} --tax_type {wildcards.tax} --fasttree_dir {input.fasttree_dir}" 

rule data_processing_non_prediction:
    input:
        file = "data/raw/{tax}_{time}_{met}.rds",
        script = "R/data_processing_non_prediction.R"
    output:
        out_file = "data/processed/{tax}_{time}_{met}_processed_noprediction.rds"
    shell:
        "RScript {input.script} --input {input.file} --output {output.out_file}"

include: "ordination.rule.py"
include: "correlation.rule.py"
include: "prediction.rule.py"








