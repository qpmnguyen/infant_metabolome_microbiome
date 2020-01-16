rule sparse_cca_plots:
    input: 
        data = "snakemake_output/analyses/correlation/{tax}_{time}_{met}_scca.rds",
        script = "R/sparse_cca_plots.R",
        pre_gen = "snakemake_output/analyses/correlation/{tax}_{time}_{met}_scc.rds"
    output:
        out_file = "snakemake_output/figures/correlation/{tax}_{time}_{met}_scca_plots.{ext}"
    shell:
        "RScript {input.script} --input {input.data} --correlation {input.pre_gen} --metab_type {wildcards.met} --time {wildcards.time} --tax_type {wildcards.tax}"

rule spearman_correlation_plots:
    input:
        data = "snakemake_output/analyses/correlation/{tax}_{time}_{met}_scc.rds",
        script = "R/spearman_corr_plots.R"
    output:
        out_file = "snakemake_output/figures/correlation/{tax}_{time}_{met}_scc_plots.{ext}"
    shell:
        "RScript {input.script} --input {input.data} --metab_type {wildcards.met} --time {wildcards.time} --tax_type {wildcards.tax}"
