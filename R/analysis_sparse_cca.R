# Sparse canonical correlation analysis 
library(PMA)
library(vegan)
library(phyloseq)
library(plyr)
library(boot)
library(picante)


sparse_cca_main <- function(data){
  tax <- as(otu_table(data),"matrix")
  met <- as(otu_table(data), "matrix")
  comb <- cbind(tax,met)
  n_cols <- ncol(tax)
  t_cols <- ncol(comb)
  
  # running the model 
  cv <- CCA.permute(x = tax, z = met, typex = "standard", typez = "standard", nperms = 25, niter = 5, standardize = T)
  cca <- CCA(x = tax, z = met, K = 1, penaltyx = cv$bestpenaltyx, penaltyz = cv$bestpenaltyz, niter = 50)
  
  # bootstrapped 
  boot <- boot(comb, cca_calculation, R = opt$n_boot)
  
  # permutation test 
  perm_test = perm_test_cca(met, tax, n_perms = opt$n_perm)
  
  # result 
  result <- list(cca = cca, boot = boot, perm = perm_test, 
                 tab = as(tax_table(data), "matrix"))
  return(result)
}


# bootstrap resamplings of the data 
cca_calculation <- function(data, idx, n_cols_first_dat = n_cols, total_cols = t_cols){
  # weird hack to combine data into one big data frame to run this function in the context 
  # of boot 
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
        rand_tax <- randomizeMatrix(tax, null.model = "richness", iterations = 1000) # one randomization
        rand_met <- randomizeMatrix(met, null.model = "richness", iterations = 1000)
        cross_val <- CCA.permute(x = rand_tax, z = rand_met, typex = "standard",typez = "standard", nperms = 25, niter = 5, standardize = T)
        cca_mod <- CCA(x = rand_tax, z = rand_met, typex = "standard", typez = "standard")
        null_corr[i] <- cca_mod$cors
    } 
    return(null_corr)
}  
