TIME = ["6W", "12M"]
MET = ["tar", "untar"]
TAX = ["16S"]

rule data_load:
    input:
        "data/data_directory.csv"
    output: 
        "data/data.rds"
    shell:
        "RScript ./R/data_load_transform.R --input {input} "