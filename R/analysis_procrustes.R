# This script takes a phyloseq object with microbiome data in the otu_table slot and metabolomics data in the sample_data
# slot, performs PCOA with calliez correction and then a procrustes comparison
# Quang Nguyen
# Updated 08/11/2020
library(vegan)
library(ade4)
library(optparse)
library(phyloseq)
library(ape)


procrustes_main <- function(data){
  tax <- data[[1]]
  met <- data[[2]]
  tax_ord <- pcoa(D = tax, correction = "cailliez")
  met_ord <- pcoa(D = met, correction = "cailliez")
  proc_test <- protest(tax_ord$vectors[,c(1,2)], met_ord$vectors[,c(1,2)])
  results <- list(tax_ord = tax_ord, met_ord = met_ord, proc_test = proc_test, label = names(data))
}