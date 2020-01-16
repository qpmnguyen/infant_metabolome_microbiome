TIME = ["6W", "12M"]
MET = ["tar", "untar"]

rule all:
    input:
        spearman_correlation = expand("output/analysis/correlation/{time}_{met}_spearman.rds", time = TIME, met = MET)
        sparse_cca = expand("output/analysis/correlation/{time}_{met}_scca.rds", time = TIME, met = MET)

rule sparse_cca:
    input: 
        data = "//dartfs-hpc/rc/lab/H/HoenA/Lab/QNguyen/ResultsFiles/processed_{time}_{met}_prediction_phyloseq_obj.rds",
        script = "R/sparse_cca.R"
    output: 
        out_file = "output/analysis/correlation/{time}_{met}_scca.rds"
    shell:
        "RScript {input.script} --input {input.data} --output {output.out_file} --n_boot 5000 --n_perm 1000"

rule spearman_correlation:
    input:
        data = "//dartfs-hpc/rc/lab/H/HoenA/Lab/QNguyen/ResultsFiles/processed_{time}_{met}_prediction_phyloseq_obj.rds",
        script = "R/spearman_corr.R"
    output: 
        out_file = "output/analysis/correlation/{time}_{met}_spearman.rds"
    shell:
        "RScript {input.script} --input {input.data} --output {output.out_file} --metric spearman --MHC BH"
