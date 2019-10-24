# Script modified Sep 16th 2019
# Quang Nguyen
# This script will do the preliminary data processing step for non-prediction tasks (aka no clr transformation)

library(vegan)
library(MiSPU)
library(phyloseq)
library(phytools)
library(optparse)
library(plyr)

option_list <- list(
  make_option("--input", help = "Input file for data loading"),
  make_option("--output", help = "Output file for data loading")
)

opt <- parse_args(OptionParser(option_list = option_list))

data <- readRDS(file = opt$input)

data <- filter_taxa(data, function(x) count(x > 0)[2,2] >= 0.1*length(x), TRUE) # asvs as to be in at least 10% of samples 
data <- tax_glom(data, taxrank = "Genus") #aggregate to genus level 
data <- transform_sample_counts(data, function(x) x/ sum(x)) # convert to relative abundance
data <- filter_taxa(data, function(x) mean(x) > 0.01e-2, TRUE) # mean of at least 0.01 

saveRDS(data, file = opt$output)