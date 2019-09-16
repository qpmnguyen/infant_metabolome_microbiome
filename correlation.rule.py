TIME = ["6W", "12M"]
MET = ["tar", "untar"]
TAX = ["16S"]
METHODS = ["rf", "svm", "enet", "spls"]
EXT = ['rds', 'svg']

rule association:
    input:
        sparse_cca_plots = expand("snakemake_output/figures/correlation/{tax}_{time}_{met}_scc_plots.{ext}", 
        tax = TAX, time = TIME, met = MET, ext = EXT),
        spearman_correlation_plots = expand("snakemake_output/figures/correlation/{tax}_{time}_{met}_scc_plots.{ext}", 
        tax = TAX, time = TIME, met = MET, ext = EXT)


rule sparse_cca:
    input: 
        data = "data/processed/{tax}_{time}_{met}_processed_noprediction.rds"
        script = "R/sparse_cca.R"
    output: 
        out_file = "snakemake_output/analyses/correlation/{tax}_{time}_{met}_scca.rds"

rule spearman_correlation:
    input:
        data = "data/processed/{tax}_{time}_{met}_processed_noprediction.rds"
    output: 
        out_file = "snakemake_output/analyses/correlation/{tax}_{time}_{met}_scc.rds"

rule sparse_cca_plots:
    input: 
        data = "snakemake_output/analyses/correlation/{tax}_{time}_{met}_scca.rds"
    output:
        out_file = "snakemake_output/figures/correlation/{tax}_{time}_{met}_scca_plots.{ext}"

rule spearman_correlation_plots:
    input:
        data = "snakemake_output/analyses/correlation/{tax}_{time}_{met}_scc.rds"
    output:
        out_file = "snakemake_output/figures/correlation/{tax}_{time}_{met}_scc_plots.{ext}"


    
