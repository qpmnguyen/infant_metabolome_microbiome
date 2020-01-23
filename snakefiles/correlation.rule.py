import platform 

TIME = ["6W", "12M"]
MET = ["tar", "untar"]

def system_info(plat):
    if plat == "Darwin":
        return(["/Volumes/rc-1/Lab/QNguyen/ResultsFiles/data/processed_{time}_{met}_prediction_phyloseq_obj.rds"])
    elif plat == "Linux":
        return(["/mnt/HoenLab/Lab/QNguyen/ResultsFiles/data/processed_{time}_{met}_prediction_phyloseq_obj.rds"])

rule all:
    input:
        spearman_correlation = expand("output/analyses/correlation/{time}_{met}_spearman.rds", time = TIME, met = MET),
        sparse_cca = expand("output/analyses/correlation/{time}_{met}_scca.rds", time = TIME, met = MET)

rule sparse_cca:
    input: 
        data = system_info(platform.system()),
        script = "R/sparse_cca.R"
    output: 
        out_file = "output/analyses/correlation/{time}_{met}_scca.rds"
    shell:
        "Rscript {input.script} --input {input.data} --output {output.out_file} --n_boot 5000 --n_perm 1000"

rule spearman_correlation:
    input:
        data = system_info(platform.system()),
        script = "R/spearman_corr.R"
    output: 
        out_file = "output/analyses/correlation/{time}_{met}_spearman.rds"
    shell:
        "Rscript {input.script} --input {input.data} --output {output.out_file} --metric spearman --MHC BH"
