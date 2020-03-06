import platform 
import os

TIME = ["6W", "12M"]
MET = ["tar", "untar"]

def system_info():
    if platform.system() == "Darwin":
        if os.path.exists("/Volumes/rc-1/Lab/") == True:
            return(["/Volumes/rc-1/Lab/QNguyen/ResultsFiles/data/processed_{time}_{met}_phyloseq_obj.rds"])
        else:
            return(["/Volumes/rc/Lab/QNguyen/ResultsFiles/data/processed_{time}_{met}_phyloseq_obj.rds"])
    elif platform.system() == "Linux":
        if platform.node() in ["polaris.dartmouth.edu", "andes.dartmouth.edu"]:
            return(["/dartfs-hpc/rc/lab/H/HoenA/Lab/QNguyen/ResultsFiles/data/processed_{time}_{met}_phyloseq_obj.rds"])
        else:
            return(["/mnt/HoenLab/Lab/QNguyen/ResultsFiles/data/processed_{time}_{met}_phyloseq_obj.rds"])
rule all:
    input:
        spearman_correlation = expand("output/analyses/correlation/{time}_{met}_spearman.rds", time = TIME, met = MET),
        sparse_cca = expand("output/analyses/correlation/{time}_{met}_scca.rds", time = TIME, met = MET)

rule sparse_cca:
    input: 
        data = system_info(),
        script = "R/sparse_cca.R"
    output: 
        out_file = "output/analyses/correlation/{time}_{met}_scca.rds"
    shell:
        "Rscript {input.script} --input {input.data} --output {output.out_file} --n_boot 5000 --n_perm 1000"

rule spearman_correlation:
    input:
        data = system_info(),
        script = "R/spearman_corr.R"
    output: 
        out_file = "output/analyses/correlation/{time}_{met}_spearman.rds"
    shell:
        "Rscript {input.script} --input {input.data} --output {output.out_file} --metric spearman --MHC BH"