# Sparse canonical correlation analysis 
library(PMA)
library(vegan)
library(phyloseq)
library(plyr)
library(boot)
library(picante)
library(rsample)

# Function to fit cca model to an rsplit object and return a tibble of results 
cca_fit <- function(split, ...){
  df <- analysis(split)
  tax_df <- df %>% select(starts_with("SV")) 
  met_df <- df %>% select(!starts_with("SV"))
  perm <- CCA.permute(tax_df, met_df, typex = "standard", typez = "standard", 
                      nperms = 25, niter = 25, standardize = T)
  mod <- CCA(x = tax_df, z = met_df, typex = "standard", typez = "standard", penaltyx = perm$bestpenaltyx, 
             penaltyz = perm$bestpenaltyz, niter = 25, standardize = T)
  result <- mod$cors
  return(result)
  
}

# permutation
perm_test_cca <- function(tax, met, n_perms){
  null_corr <- c()
  for (i in 1:n_perms){
    rand_tax <- randomizeMatrix(tax, null.model = "richness", iterations = 1000) # one randomization
    rand_met <- randomizeMatrix(met, null.model = "richness", iterations = 1000)
    perm <- CCA.permute(x = rand_tax, z = rand_met, typex = "standard",typez = "standard", 
                             nperms = 50, niter = 25, standardize = T)
    cca_mod <- CCA(x = rand_tax, z = rand_met, typex = "standard", typez = "standard", 
                   niter = 25, penaltyx = perm$bestpenaltyx, 
                   penaltyz = perm$bestpenaltyz, standardize = T)
    null_corr[i] <- cca_mod$cors
  } 
  return(null_corr)
}  

sparse_cca_main <- function(data, n_boot, n_perm){
  tax <- as(otu_table(data),"matrix")
  met <- as(sample_data(data), "matrix")
  comb <- cbind(tax,met) %>% as.data.frame()
  
  print("Running the model...")
  # running the model 
  cv <- CCA.permute(x = tax, z = met, typex = "standard", typez = "standard", 
                    nperms = 50, niter = 25, standardize = T)
  cca <- CCA(x = tax, z = met, K = 1, penaltyx = cv$bestpenaltyx, 
             penaltyz = cv$bestpenaltyz, niter = 25)
  
  print("Estimating bootstrapped confidence intervals")
  # bootstrapped 
  
  boot <- rsample::bootstraps(comb, times = n_boot) %>% 
    mutate(cors = map(splits, ~ cca_fit(.x))) %>% unnest(cors) %>% select(c(id,cors))
  
  # permutation test 
  perm_test <- perm_test_cca(met, tax, n_perms = n_perm)
  
  # result 
  result <- list(cca = cca, boot = boot, perm = perm_test, 
                 tab = as(tax_table(data), "matrix"))
  return(result)
}

