library(optparse)
library(ggplot2)
library(viridis)
library(ggpubr)
library(MLmetrics)
library(cowplot)

options_list <- list(
  make_option("--model", help = "Selected model for evaluation"),
  make_option("--metab_type", help = "Type of metabolomics data"),
  make_option("--time", help = "Data time point to extract"),
  make_option("--tax_type", help = "Data type for taxonomic data, so far only '16S' is supported")
)