TIME = ["6W", "12M"]
MET = ["tar", "untar"]
TAX = ["16S"]
METHODS = ["rf", "svm", "enet", "spls"]

rule all:
    input:
        processed = expand('data/processed/{tax}_{time}_{met}_processed.rds', tax = TAX, time = TIME, met = MET)


rule data_retrieval:
    input:
        dir_file = "data/data_directory.csv",
        script = "R/data_retrieval.R",
        fasttree_dir = "../software/"
    output: 
        out_file = "data/raw/{tax}_{time}_{met}.rds"
    shell:
        "RScript {input.script} --input {input.dir_file} --output {output.out_file} --time {wildcards.time} --metab_type {wildcards.met} --tax_type {wildcards.tax} --fasttree_dir {input.fasttree_dir}"  

rule prediction_processing:
    input:
        data = "data/raw/{tax}_{time}_{met}.rds",
        script = "R/prediction_processing.R"
    output:
        out_file = "data/processed/{tax}_{time}_{met}_processed.rds"
    shell:
        "RScript {input.script} --input {input.data} --output {output.out_file} --metab.type {wildcards.met}"

rule ordination_analyses:
    input: 
        data = "data/raw/{tax}_{time}_{met}.rds",
        script = "R/ordination_analyses.R"
    output:
        out_file = "analyses/{tax}_{time}_{met}_ordinations.rds"
    shell:
       "echo {output.out_file}"





