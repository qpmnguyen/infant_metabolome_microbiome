# Script to generate positive and negative simulations 
# by Quang Nguyen
# Last edit: 08/05/2019

library(MCTools)
library(phyloseq)
library(SpiecEasi)
library(compositions)
library(picante)
library(caret)
library(MLmetrics)
source("./scripts/utils.R")
# grab appropriate data and perform preliminary processing  
data <- load_main(file = "./data/data_directory.csv")
tax_data <- rbind(data$tax.12M, data$tax.6W)
key <- key_processing(data$tax.key,type = "NA")
phylo_obj <- phyloseq(otu_table(tax_data, taxa_are_rows = F),
                      tax_table(key))
phylo_obj <- tax_glom(phylo_obj, taxrank = "Genus")

phylo_obj <- prune_taxa(taxa_sums(phylo_obj) > 0, phylo_obj)

tax <- as.matrix(otu_table(phylo_obj))

met <- rbind(data$metabo.12M, data$metabo.12M)

#' @title Generate simulated data calibrated from real data set  
#' @description Using simulating mechanisms from SpiecEasi package to generate calibrated 
#'   simualted microbiome data set with a log-normal outcome. Outcome can be non-associated
#'   or associated.  
#' @param n_datasets: Number of datasets 
#' @param outcome_type: Outcome type can be "assc" for associative or "none" for unassociated
#' @param calibrate: calibration count data set. Taxa as columns and samples as rows  
#' @param correlation_type: Correlation type for generating data, the default is "cluster"
#' @param snr: Signal to noise ratio  
#' @param sample_size: Sample Size 
#' @param perc_assc: Proportion associated with outcome for "assc" outcome_type 
#' @param beta_val: The value of beta coefficients, default to 6/sqrt(10) and constant across all 
#'   taxa. 
#' @param calibrate_outcome: Outcome data (features as columns and samples as rows) to calibrate 
#'   unassociated outcome. If calibrate_outcome is NULL, use meanlog = 0 and sdlog = 1. 
generate_data <- function(n_datasets, outcome_type, calibrate, correlation_type = "cluster", snr, 
                          sample_size, perc_assc = 0.05, 
                          beta_val = 6/sqrt(10), calibrate_outcome = NULL){
  # generate graph for correlation structure from SpiecEasi  
  graph <- make_graph(method = correlation_type, D = ncol(calibrate), e = ncol(calibrate))
  Prec  <- graph2prec(graph)  
  Cor   <- cov2cor(prec2cov(Prec))
  # use SpiecEasi functions to generate data  
  results <- list()
  depth <- median(rowSums(calibrate))
  calibrate <- unclass(acomp(calibrate))
  calibrate <- round(calibrate * depth)
  for(i in 1:n_datasets){
    sim_tax <- synth_comm_from_counts(calibrate, mar=2, distr='zinegbin', Sigma=Cor, n=sample_size, retParams = F) # faster than doing it by hand  
    # Generate outcomes in naive situation with all taxa having the same regression coefficients  
    ass_idx <- sample(1:ncol(calibrate), perc_assc * ncol(calibrate), replace = F)
    coef <- rep(0, ncol(calibrate))
    coef[ass_idx] <- beta_val
    norm_sim_tax <- unclass(acomp(sim_tax))
    if (outcome_type == "assc"){
      outcome <- 1 + sim_tax %*% coef # identical intercepts 
      error <- rlnorm(sample_size, meanlog = 0, sdlog = 1) # standard log normal error 
      k <- sqrt(var(outcome)/(snr*var(error))) # snr multiplier
      suppressWarnings(new_outcome <- outcome + k * error)
    }
    else {
      if (is.null(calibrate_outcome)){
        outcome = rlnorm(sample_size, meanlog = 0, sdlog = 1)
        error <- rlnorm(sample_size, meanlog = 0, sdlog = 1) # standard log normal error 
        k <- sqrt(var(outcome)/(snr*var(error))) # snr multiplier
        suppressWarnings(new_outcome <- outcome + k * error)
      } else {
        fit <- fitdist(calibrate_outcome[calibrate_outcome > 0], distr = "lnorm", method = "mle")
        outcome = rlnorm(sample_size, meanlog = fit$estimate[1], sdlog = fit$estimate[2])
        error <- rlnorm(sample_size, meanlog = 0, sdlog = 1) # standard log normal error 
        k <- sqrt(var(outcome)/(snr*var(error))) # snr multiplier
        suppressWarnings(new_outcome <- outcome + k * error)
      }
    }
    results[[i]] <- list(tax = sim_tax, 
                         coef = coef, 
                         norm_sim_tax = norm_sim_tax, 
                         org_outcome = outcome, 
                         met = new_outcome)
  }
  return(results)
}

associated_results <- generate_data(n_datasets = 50, outcome_type = "assc", calibrate = tax, snr = 2, sample_size = 300)
unassociated_results <- generate_data(n_datasets = 50, outcome_type = "none", calibrate = tax, snr = 2, sample_size = 300, 
                                      calibrate_outcome = as.vector(met))
if (!dir.exists("./simulated_data/")){
  dir.create("./simulated_data/")
}
if (!dir.exists("./simulated_data/positive_sim/")){
  dir.create('./simulated_data/positive_sim/')
}
if (!dir.exists("./simulated_data/negative_sim/")){
  dir.create('./simulated_data/negative_sim/')
}

sapply(1:50, function(x){
  obj <- associated_results[[x]]
  path = paste0("./simulated_data/positive_sim/", "pos_sim_dataset",x,".rds")
  saveRDS(obj, file = path)
})

sapply(1:50, function(x){
  obj <- unassociated_results[[x]]
  path = paste0("./simulated_data/negative_sim/", "neg_sim_dataset",x,".rds")
  saveRDS(obj, file = path)
})

### Create new simulations using permutations  
if (!dir.exists('./simulated_data/permute_sim/')){
  dir.create('./simulated_data/permute_sim')
}

for (i in 1:1000){
  rand_tax <- randomizeMatrix(tax, null.model = "richness", iterations = 1000)
  rand_tax <- unclass(acomp(rand_tax))
  rand_tax[rand_tax == 0] <- 1
  rand_tax <- unclass(clr(rand_tax))
  rand_met <- randomizeMatrix(met, null.model = "richness", iterations = 1000)
  export <- list(tax = rand_tax, met = rand_met)
  saveRDS(export, file = paste0("./simulated_data/permute_sim/randiter_",i,".rds"))
}

