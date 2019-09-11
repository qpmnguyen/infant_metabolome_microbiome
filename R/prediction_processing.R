# Script written 09/11/19
# Quang Nguyen
# This script is used to process data for prediction models 


library(phyloseq)
library(optparse)
library(compositions)

option_list <- list(
  make_option("--input", help = "Input file for data processing"),
  make_option("--output", help = "Output file for data processing"),
  make_option("--metab_type", help = "Data type for metabolite data, can be 'tar' or 'untar'")
)

opt <- parse_args(OptionParser(option_list = option_list))


data <- readRDS(file = opt$input)

# processing phylogenetic data 
data <- filter_taxa(data, function(x) count(x > 0)[2,2] >= 0.1*length(x), TRUE) # asvs as to be in at least 10% of samples 
otu_table(data) <- otu_table(data) + 1 # adding pseudocount
data <- tax_glom(data, taxrank = "Genus") #aggregate to genus level 
data <- transform_sample_counts(data, function(x) x/ sum(x)) # convert to relative abundance
data <- filter_taxa(data, function(x) mean(x) > 0.01e-2, TRUE) # mean of at least 0.01 
tax <- as(otu_table(data), "matrix")
tax <- unclass(clr(tax)) # centerd log_ratio transformation 

# processing metabolomics data 
met <- as(sample_data(data), "matrix")
if (opt$metab_type == "tar"){
    met <- log(met + 1) # log x+1 transformation
} else if (opt$metab_type == "untar"){
    met <- unclass(acomp(met)) # renormalize to relative abundances
    met <- asin(sqrt(met)) # arcsine square root transformation 
}

output <- list(met = met, tax = tax, tab = tax_table(data), tree = phy_tree(data))
saveRDS(output, file = opt$output)


