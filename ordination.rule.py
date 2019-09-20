TIME = ["6W", "12M"]
MET = ["tar", "untar"]
TAX = ["16S"]
EXT = ['rds', 'svg']

rule ordination:
    input:
        joint_ordinations = expand('snakemake_output/figures/ordinations/{tax}/{time}/{met}/joint_protest_ordinations.{ext}', 
        tax = TAX, time = TIME, met = MET, ext = EXT),
        met_ordinations = expand("snakemake_output/figures/ordinations/{tax}/{time}/{met}/met_ordinations.{ext}", 
        tax = TAX, time = TIME, met = MET, ext = EXT),
        tax_ordinations = expand("snakemake_output/figures/ordinations/{tax}/{time}/{met}/tax_ordinations.{ext}", 
        tax = TAX, time = TIME, met = MET, ext = EXT)

rule clean_ordination:
    shell: 
        """
        rm -rf figures/ordinations/
        rm -rf analyses/ordinations/ 
        """

rule constructing_distance_matrices:
    input:
        data = "data/processed/{tax}_{time}_{met}_processed_noprediction.rds",
        script = "R/constructing_distance_matrices.R"
    output:
        out_file = "snakemake_output/analyses/ordinations/{tax}_{time}_{met}_distance.rds"
    shell:
        "RScript {input.script} --input {input.data} --output {output.out_file}"

rule procrustes_analyses: 
    input: 
        data = "snakemake_output/analyses/ordinations/{tax}_{time}_{met}_distance.rds",
        script = "R/procrustes_analyses.R"
    output:
        out_file = "snakemake_output/analyses/ordinations/{tax}_{time}_{met}_ordinations.rds"
    shell:
        "RScript {input.script} --input {input.data} --output {output.out_file}"

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
