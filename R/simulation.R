library(optparse)
library(phyloseq)
library(SpiecEasi)
library(compositions)
library(picante)
library(caret)
library(MLmetrics)

option_list <- list(
  make_option("--calibrate_data_1", help = "data set 1 to calibrate to"),
  make_option("--calibrate_data_2", help = "data set 2 to calibrate to"),
  make_option("--outcome_type", help = "assc for positive simulations, none for unassociated outcomes"),
  make_option("--n_datasets", help = "Number of data sets"),
  make_option("--snr", help = "Signal to noise ratio"),
  make_option("--sample_size", help = "Sample size per data set"),
  make_option("--perc_assc", help = "Proportion of taxa associated with outcome if outcome_type is assc"),
  make_option("--calibrate_outcome", help = "Logical, indiciate whether or not to calibrate unasscoiated outcomes")
)

opt <- parse_args(OptionParser(option_list = option_list))

####-------------------Function definitions--------------------------------------------------------
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

calibrate_data <- readRDS(file = opt$calibrate_data)
tax <- calibrate_data$tax
met <- calibrate_data$met

if (opt$calibrate_outcome == T){
  if (opt$outcome_type == "assc"){
    stop("calibrate outcome can't be TRUE if outcome type is associative")
  } else {
    calibrate_outcome = as.vector(met)
  }
} else {
  calibrate_outcome <- NULL
}

datasets <- generate_data(n_datasets = opt$n_datasets, perc_assc = opt$perc_assc, outcome_type = opt$outcome_type, calibrate = tax, 
                          snr = opt$snr, sample_size = opt$sample_size, calibrate_outcome = calibrate_outcome)

output <- list(datasets = datasets, params = list(
  n_datasets = opt$n_datasets, 
  perc_assc = opt$perc_assc,
  outcome_type = opt$outcome_type, 
  snr = opt$snr, 
  sample_size = opt$sample_size,
  calibrate_outcome = opt$calibrate_outcome
))

saveRDS(output, file = paste0("snakemake_output/simulation/", opt$outcome_type,"_",opt$snr,"_",opt$perc_assc,"_overall_dataset.rds"))
sapply(1:opt$n_datasets, function(x){
  obj <- datasets[[x]]
  path = paste0("snakemake_output/simulation/single_obj/", opt$outcome_type, "/",opt$snr,"/",opt$perc_assc, "/")
  if(!dir.exists(path)){
    dir.create(path, recursive = T)
  }
  saveRDS(obj, file = paste0(path, "simulated_", x, ".rds"))
})
