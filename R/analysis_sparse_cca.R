# Sparse canonical correlation analysis 
library(PMA)
library(vegan)
library(phyloseq)
library(optparse)
library(plyr)
library(boot)
library(picante)


option_list <- list(
  make_option("--input", help = "Input file for data loading", type = "character"),
  make_option("--output", help = "Output file for data loading", type = "character"),
  make_option("--n_boot", help = "Number of bootstrap iterations", type = "integer"),
  make_option("--n_perm", help = "Number of permutations", type = "integer")
)

opt <- parse_args(OptionParser(option_list = option_list))

data <- readRDS(file = opt$input)

# some pre-processing to make sure everything is in the same format
tax <- otu_table(data)
met <- sample_data(data)

# adding two data frames together and calculate n_cols 
comb <- cbind(tax, met)
n_cols <- ncol(tax)
t_cols <- ncol(comb)

# bootstrap resamplings of the data 
cca_calculation <- function(data, idx, n_cols_first_dat = n_cols, total_cols = t_cols){
    tax <- data[,1: n_cols_first_dat]
    met <- data[,(n_cols_first_dat + 1) : total_cols]
    cross_val <- CCA.permute(x = tax[idx,], z = met[idx,], typex = "standard", typez = "standard", nperms = 25, niter = 5, standardize = T)
    cca_mod <- CCA(x = tax[idx,], z = met[idx,], typex = "standard", typez = "standard")
    return(cca_mod$cors)
}


# permutation test 
perm_test_cca <- function(tax, met,  n_perms){
    null_corr <- c()
    for (i in 1:n_perms){
        rand_tax <- randomizeMatrix(tax, null.model = "richness", iterations = 1000)
        rand_met <- randomizeMatrix(met, null.model = "richness", iterations = 1000)
        cross_val <- CCA.permute(x = rand_tax, z = rand_met, typex = "standard",typez = "standard", nperms = 25, niter = 5, standardize = T)
        cca_mod <- CCA(x = rand_tax, z = rand_met, typex = "standard", typez = "standard")
        null_corr[i] <- cca_mod$cors
    } 
    return(null_corr)
}  

# Running the model 
cv <- CCA.permute(x = tax, z = met, typex = "standard", typez = "standard", nperms = 25, niter = 5, standardize = T)
cca <- CCA(x = tax, z = met, K = 1, penaltyx = cv$bestpenaltyx, penaltyz = cv$bestpenaltyz, niter = 50)

# bootstrapped 
boot <- boot(comb, cca_calculation, R = opt$n_boot)

# permutation test 
perm_test = perm_test_cca(met, tax, n_perms = opt$n_perm)
result <- list(cca = cca, boot = boot, perm_test, tab = tax_table(data))

saveRDS(result, file = opt$output)