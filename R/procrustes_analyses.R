library(vegan)
library(ade4)
library(optparse)

option_list <- list(
  make_option("--input", help = "Input file for data loading"),
  make_option("--output", help = "Output file for data loading")
)

opt <- parse_args(OptionParser(option_list = option_list))

# accessing dataset 
data <- readRDS(file = opt$input)
tax_dist <- data$tax_dist
met_dist <- data$met_dist

# Performing ordinations  
tax_ord <- metaMDS(comm = tax_dist, try = 50, engine = "isoMDS")
met_ord <- metaMDS(comm = as.matrix(met_dist), try = 50, engine = "isoMDS")

# Procrustes analyses  
proc_test <- protest(tax_ord, met_ord)
mant_test <- mantel(tax_dist, as.matrix(met_dist), permutations = 9999)

# saving results  
results <- list(tax_ord = tax_ord, met_ord = met_ord, proc_test = proc_test,
                    mant_test = mant_test)

saveRDS(results, file = opt$output)
