# Reproducible Analyses for manuscript "Healthy Infant Metabolomes are Robust to Normal Variations in Microbiome"

This repository hosts reproducible workflow for analyses represented in the manuscript "Healthy Infant Metabolomes are Robust to Normal Variations in Microbiome" by Nguyen et al. 2019. 

Generated Rmarkdown and figures are located in the `docs/` folder. All scripts are located in the `R` folder. An `renv.lock` file is generated using the `renv` package to manage R package dependencies. A `Snakefile` is generated so all reports/images can be generated automatically with appropriate log files (located in `logs/` folder) and DAG (located in the `docs\` folder)

You can run these results using `snakemake` and the `Snakefile` with your own data. 
