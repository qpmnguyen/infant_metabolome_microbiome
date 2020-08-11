# Infant Gut Metabolome is Robust to Stochastic Variations in the Microbiome

This repository hosts reproducible scripts for analyses presented in the manuscript workflow for analyses represented in the manuscript "Healthy Infant Metabolomes are Robust to Normal Variations in Microbiome" by Nguyen et al. 2020 currently in prep.  

Scripts are organized into modularized functions and fed into workflows designed in `drake`. call the `make_<analysis>.R` file to execute the workflow.    

No data is currently provided as the manuscript is being prepared. Dependencies can be accessed using the `dependencies.csv` file. All `Bioconductor` packages were installed under `Bioconductor` version 3.10. R version is `3.6.3`.  

Currently working on creating a new docker file for increased reproducibility. First create an image:  
```bash
docker build --tag infant_robust:1.0 .
```
Then boot into interactive mode with RStudio at port 8787
```bash
docker run --rm -p 8787:8787 -e PASSWORD=password --name infant_robust infant_robust:1.0
```
