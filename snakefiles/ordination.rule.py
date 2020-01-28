import platform 

TIME = ["6W", "12M"]
MET = ["tar", "untar"]

def system_info(plat):
    if plat == "Darwin":
        if os.path.exists("/Volumes/rc-1/Lab/") == True:
            return(["/Volumes/rc-1/Lab/QNguyen/ResultsFiles/data/processed_{time}_{met}_phyloseq_obj.rds"])
        else:
            return(["/Volumes/rc/Lab/QNguyen/ResultsFiles/data/processed_{time}_{met}_phyloseq_obj.rds"])
    elif plat == "Linux":
        if platform.node() in ["polaris.dartmouth.edu", "andes.dartmouth.edu"]:
            return(["/dartfs-hpc/rc/lab/H/HoenA/Lab/QNguyen/ResultsFiles/data/processed_{time}_{met}_phyloseq_obj.rds"])
        else:
            return(["/mnt/HoenLab/Lab/QNguyen/ResultsFiles/data/processed_{time}_{met}_phyloseq_obj.rds"])
rule ordination:
    input:
        distance = expand("output/analyses/ordinations/{time}_{met}_distance.rds", time = TIME, met = MET),
        ordinations = expand("output/analyses/ordinations/{time}_{met}_ordinations.rds", time = TIME, met = MET)

rule constructing_distance_matrices:
    input:
        data = system_info(plat = platform.system()),
        script = "R/constructing_distance_matrices.R"
    output:
        out_file = "output/analyses/ordinations/{time}_{met}_distance.rds"
    shell:
        "Rscript {input.script} --input {input.data} --output {output.out_file}"

rule procrustes_analyses: 
    input: 
        data = "output/analyses/ordinations/{time}_{met}_distance.rds",
        script = "R/procrustes_analyses.R"
    output:
        out_file = "output/analyses/ordinations/{time}_{met}_ordinations.rds"
    shell:
        "Rscript {input.script} --input {input.data} --output {output.out_file}"


