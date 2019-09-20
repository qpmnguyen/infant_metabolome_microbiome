TIME = ["6W", "12M"]
MET = ["tar", "untar"]
TAX = ["16S"]
METHODS = ["rf", "svm", "enet", "spls"]
EXT = ['rds', 'svg']

localrules: prediction_processing, prediction_eval, prediction_plots

rule prediction:
    input:
        eval_plots = expand("snakemake_output/figures/prediction/{tax}_{time}_{met}_{methods}_prediction_plots.{ext}", tax = TAX, 
        time = TIME, met = MET, methods = METHODS, ext = EXT) 

rule clean_prediction:
    shell:
        """
        rm -rf analyses/prediction/
        rm -rf figures/prediction/ 
        """

rule prediction_processing:
    input:
        data = "data/raw/{tax}_{time}_{met}.rds",
        script = "R/data_processing_prediction.R"
    output:
        out_file = "data/processed/{tax}_{time}_{met}_processed_prediction.rds"
    shell:
        "RScript {input.script} --input {input.data} --output {output.out_file} --metab_type {wildcards.met}"

rule cross_validation:
    input:
        data = "data/processed/{tax}_{time}_{met}_processed_prediction.rds",
        script = "R/cross_validation.R",
        model_script = "R/modelfit.R"
    output:
        out_file = "snakemake_output/analyses/prediction/{tax}_{time}_{met}_{methods}_prediction_values.rds"
    shell:
        "RScript {input.script} --model {input.model_script} --input {input.data}"

rule prediction_eval:
    input: 
        data = "snakemake_output/analyses/prediction/{tax}_{time}_{met}_{methods}_prediction_values.rds",
        script = "R/prediction_eval.R"
    output:
        out_file = "snakemake_output/analyses/prediction/{tax}_{time}_{met}_{methods}_prediction_eval.rds"
    shell:
        "RScript {input.script} --input {input.data}"
    

rule prediction_plots:
    input:
        data = "snakemake_output/analyses/prediction/{tax}_{time}_{met}_{methods}_prediction_eval.rds",
        script = "R/prediction_plots.R"
    output:
        out_file = "snakemake_output/figures/prediction/{tax}_{time}_{met}_{methods}_prediction_plots.{ext}"
    shell:
        "RScript {input.script} --input {input.data}"