EXT = ["svg", "png"]
TIME = ["6W", "12M"]
METAB = ["tar", "untar"]
TAX_DIST = ["gunifrac", "euclidean"]
MET_DIST = ["euclidean", "manhattan"]

rule all:
    input: 
        joint_ordination = expand("output/figures/ordinations/{tax_dist}_{met_dist}/{time}_{metab}_joint.{ext}", 
        time = TIME, metab = METAB, tax_dist = TAX_DIST, met_dist = MET_DIST, ext = EXT),
        tax_ordination = expand("output/figures/ordinations/{tax_dist}_{met_dist}/{time}_{metab}_tax.{ext}",
        time = TIME, metab = METAB, tax_dist = TAX_DIST, met_dist = MET_DIST, ext = EXT),
        met_ordination = expand("output/figures/ordinations/{tax_dist}_{met_dist}/{time}_{metab}_met.{ext}",
        time = TIME, metab = METAB, tax_dist = TAX_DIST, met_dist = MET_DIST, ext = EXT)


rule ordination_plots: 
     input: 
         data = "output/analyses/ordinations/{time}_{metab}_ordinations.rds",
         script = "R/plots_ordination.R"
     output:
         joint_plot = "output/figures/ordinations/{tax_dist}_{met_dist}/{time}_{metab}_joint.{ext}",
         met_plot = "output/figures/ordinations/{tax_dist}_{met_dist}/{time}_{metab}_met.{ext}",
         tax_plot = "output/figures/ordinations/{tax_dist}_{met_dist}/{time}_{metab}_tax.{ext}"
     shell:
         "Rscript {input.script} --input {input.data} --tax_dist {wildcards.tax_dist} --met_dist {wildcards.met_dist}"



# rule sparse_cca_plots:
#     input: 
#         data = "snakemake_output/analyses/correlation/{tax}_{time}_{met}_scca.rds",
#         script = "R/sparse_cca_plots.R",
#         pre_gen = "snakemake_output/analyses/correlation/{tax}_{time}_{met}_scc.rds"
#     output:
#         out_file = "snakemake_output/figures/correlation/{tax}_{time}_{met}_scca_plots.{ext}"
#     shell:
#         "RScript {input.script} --input {input.data} --correlation {input.pre_gen} --metab_type {wildcards.met} --time {wildcards.time} --tax_type {wildcards.tax}"
# 
# rule spearman_correlation_plots:
#     input:
#         data = "snakemake_output/analyses/correlation/{tax}_{time}_{met}_scc.rds",
#         script = "R/spearman_corr_plots.R"
#     output:
#         out_file = "snakemake_output/figures/correlation/{tax}_{time}_{met}_scc_plots.{ext}"
#     shell:
#         "RScript {input.script} --input {input.data} --metab_type {wildcards.met} --time {wildcards.time} --tax_type {wildcards.tax}"
# 

