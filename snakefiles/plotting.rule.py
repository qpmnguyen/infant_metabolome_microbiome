EXT = ["svg", "png", "rds"]
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
        scca_plot = expand("output/figures/correlation/{time}_{metab}_scca.{ext}"), 
        spearman_plot = "output/figures/correlation/{time}_{metab}_correlation.{ext}"


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

rule correlation_plots:
    input:
        scca = "output/analyses/correlation/{time}_{metab}_scca.rds",
        correlation = "output/analyses/correlation/{time}_{metab}_spearman.rds",
        script = "R/plots_correlation.R"
    output:
        scca_plot = "output/figures/correlation/{time}_{metab}_scca.{ext}",
        spearman_plot = "output/figures/correlation/{time}_{metab}_correlation.{ext}"
    shell:
        "Rscript {input.script} --scca {input.scca} --correlation {input.correlation}"
