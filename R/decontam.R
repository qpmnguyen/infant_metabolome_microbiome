# DECONTAM and preliminary analyses  
# Quang Nguyen
# Updated Jan 20th 2020
library(decontam)
library(phyloseq)


###---Step 1: Load data files and concentration ---###
data_repo <- read.csv(file = "/mnt/HoenLab/Lab/QNguyen/ResultsFiles/data_directory.csv", stringsAsFactors = F)
asv_table <- readRDS(file = data_repo$directory[data_repo$key == "tax6W"])
tax_table <- readRDS(file = data_repo$directory[data_repo$key == "tax6Wkey"])
inventory <- read.csv(file = paste0("/mnt/HoenLab",data_repo$directory[data_repo$key == "tax_inventory"]))
###---Step 2: Matching samples from asv_table to inventory and extract DNA concentrations using mblids---###

