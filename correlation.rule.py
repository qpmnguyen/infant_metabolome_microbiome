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
        data = "data/processed/{tax}_{time}_{met}_processed_prediction.rds",
        script = "R/sparse_cca.R"
    output: 
        out_file = "snakemake_output/analyses/correlation/{tax}_{time}_{met}_scca.rds"
    shell:
        "RScript {input.script} --input {input.data} --output = {output.out_file}"

rule spearman_correlation:
    input:
        data = "data/processed/{tax}_{time}_{met}_processed_prediction.rds",
        script = "R/spearman_corr.R"
    output: 
        out_file = "snakemake_output/analyses/correlation/{tax}_{time}_{met}_scc.rds"

rule sparse_cca_plots:
    input: 
        data = "snakemake_output/analyses/correlation/{tax}_{time}_{met}_scca.rds",
        script = "R/sparse_cca_plots.R"
    output:
        out_file = "snakemake_output/figures/correlation/{tax}_{time}_{met}_scca_plots.{ext}"
    shell:
        "RScript {input.script} --input {input.data} --metab_type {wildcards.met} --time {wildcards.time} --tax_type {wildcards.tax}"

rule spearman_correlation_plots:
    input:
        data = "snakemake_output/analyses/correlation/{tax}_{time}_{met}_scc.rds",
        script = "R/spearman_corr_plots.R"
    output:
        out_file = "snakemake_output/figures/correlation/{tax}_{time}_{met}_scc_plots.{ext}"
