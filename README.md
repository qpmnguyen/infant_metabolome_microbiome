# Infant Gut Metabolome is Robust to Stochastic Variations in the Microbiome

This repository hosts reproducible scripts for analyses presented in the manuscript workflow for analyses represented in the manuscript "Healthy Infant Metabolomes are Robust to Normal Variations in Microbiome" by Nguyen et al. 2020 currently in prep.  

Scripts are all modularized. Argument help can be accessed by running `Rscript script_name.R --help`. Workflow will be described in detail in paper with the following general schematic. The easiest way is to set up `snakemake` scripts.  

No data is currently provided as the manuscript is being prepared. Dependencies can be accessed using the `renv.lock` file through the `renv` package in R. 
  
```r
install.packages('renv')
renv::restore()
```  
The most reproducible method is to use the attached Dockefile. First create an image:  
```bash
docker build --tag infant_robust:1.0 .
```
Then boot into interactive mode with RStudio at port 8787
```bash
docker run --rm -p 8787:8787 -e PASSWORD=password --name infant_robust infant_robust:1.0
```
