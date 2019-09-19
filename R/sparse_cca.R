library(PMA)
library(vegan)
library(phyloseq)
library(optparse)
library(plyr)
library(boot)
library(picante)


option_list <- list(
  make_option("--input", help = "Input file for data loading"),
  make_option("--output", help = "Output file for data loading")
)


data <- readRDS(file = opt$input)

# some pre-processing to make sure everything is in the same format
tax <- as(otu_table(data), "matrix")
met <- as(sample_data(data), "matrix") # this procedure turns everything into character
temp_names <- rownames(met)
met <- apply(met, 2, as.numeric)
rownames(met) <- temp_names
rm(temp_names)

# adding two data frames together and calculate n_cols 
comb <- cbind(tax, met)
n_cols <- ncol(tax)
total_cols <- ncol(comb)

# running the model  
cv <- CCA.permute(x = tax, z = met, typex = "standard", typez = "standard", nperms = 50, niter = 5, standardize = T)
cca <- CCA(x = tax, z = met, K = 1, penaltyx = cv$bestpenaltyx, penaltyz = cv$bestpenaltyz, niter = 50)

# bootstrap resamplings of the data 
cca_calculation <- function(data, idx, n_cols_first_dat, total_cols){
    tax <- data[,1: n_cols_first_dat]
    met <- data[,(n_cols_first_dat + 1) : total_cols]
    cross_val <- CCA.permute(x = tax[idx,], z = met[idx,], typex = "standard", typez = "standard", nperms = 50, niter = 5, standardize = T)
    cca_mod <- CCA(x = tax[idx,], z = met[idx,], typex = "standard", typez = "standard")
    return(cca_mod$cors)
}

# bootstrapped 
boot <- boot(comb, cca_calculation, R = 999)

# permutation test 
perm_test_cca <- function(tax, met,  n_perms){
    null_corr <- c()
    for (i in 1:n_perms){
        rand_tax <- randomizeMatrix(tax, null.model = "richness", iterations = 1000)
        rand_met <- randomizeMatrix(met, null.model = "richness", iterations = 1000)
        cross_val <- CCA.permute(x = rand_tax, z = rand_met, typex = "standard",typez = "standard", nperms = 50, niter = 5, standardize = T)
        cca_mod <- CCA(x = rand_tax, z = rand_met, typex = "standard", typez = "standard")
        null_corr[i] <- cca_mod$cors
    } 
    return(null_corr)
}  

perm_test = perm_test_cca(met, tax, n_perms = 999)
result <- list(cca = cca, boot = boot, perm_test)

saveRDS(result, file = opt$output)

