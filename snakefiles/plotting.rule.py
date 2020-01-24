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

rule ordination_plots: 
    input: 
        data = "snakemake_output/analyses/ordinations/{tax}_{time}_{met}_ordinations.rds",
        script = "R/ordination_plots.R"
    output:
        joint_plot = "snakemake_output/figures/ordinations/{tax}/{time}/{met}/joint_protest_ordinations.{ext}",
        met_plot = "snakemake_output/figures/ordinations/{tax}/{time}/{met}/met_ordinations.{ext}",
        tax_plot = "snakemake_output/figures/ordinations/{tax}/{time}/{met}/tax_ordinations.{ext}"
    shell:
        "Rscript {input.script} --input {input.data} --metab_type {wildcards.met} --time {wildcards.time} --tax_type {wildcards.tax}"
