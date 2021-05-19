# Associations between the gut microbiome and metabolome in early life

This repository hosts reproducible scripts for analyses presented in the manuscript workflow for analyses represented in the manuscript "Associations between the gut microbiome and metabolome in early life" by Nguyen et al. 2021+.   

Scripts are organized into modularized functions and fed into workflows designed in `drake`. call the `make(plan_name)` in the `make_<analysis>.R` file to execute the workflow. The only exception are the prediction models and Bayesian evaluations, since they are computational intensive and were ran in parallel on Dartmouth's Discovery cluster.  

For supplemental picrust2 analyses, pipeline was ran using Snakemake, all files and folders in the `picrust2` folder.  

`plots_` scripts hold all plotting functions.      
`analysis_` scripts hold all analysis functions.       
`process_` scripts hold all functions that process outputs from analyses, mostly in service of plotting.    
`plan_` scripts hold all `drake` plans for streamlining the dispatch of complex analysis workflow across multiple targets.    

Raw data deposited on their respective repositories on NCBI SRA and Metabolomics Data Repository. Dependencies can be accessed using the `dependencies.csv` file. All `Bioconductor` packages were installed under `Bioconductor` version 3.10. R version is `3.6.3`.  

Currently working on creating a new docker file for increased reproducibility. 

To use Docker (WIP) first create an image:  
```bash
docker build --tag infant_microbiome_metabolome:1.0 .
```
Then boot into interactive mode with RStudio at port 8787
```bash
docker run --rm -p 8787:8787 -e PASSWORD=password --name infant_robust infant_microbiome_metabolome:1.0
```
