# Model posterior distribution of test statistic based on model 
library(tidyposterior)
library(tidyverse)
library(rsample)
library(glue)
library(optparse)
library(doParallel)
library(parallel)
library(foreach)

option_list <- list(
  make_option("--eval", help = "Evaluation metric"),
  make_option("--timepoint", help = "Time point of analysis"),
  make_option("--metabtype", help = "tar or untar"),
  make_option("--parallel", help = "Check if something is parallel"),
  make_option("--ncores", type = "integer", default = NULL, help = "Number of cores"),
  make_option("--output_dir", help = "output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))



#' @title Bayesian modelling of evaluation metrics 
#' @param summary_list Summary list after processing raw data
#' @param time Timepoint to subset
#' @param parallel Check whether paralle is true or false
#' @param ncores number of cores to use if parallel
#' @description Using \code{tidyposterior}, first a subsetted dataframe is converted to proper \code{rset}
#'   object, which is then funelled through the appropriate bayesian hierarchical model. 
get_bayesian <- function(summary_list, time, parallel, ncores = NULL){
  met_names <- names(summary_list)
  if (parallel == T){
    print("Processing parallel the bayesian model per metabolite")
    if (is.null(ncores)){
      ncores <- detectCores() - 2
    }
    # registering the cluster
    cluster <- makeForkCluster(ncores)
    registerDoParallel(cluster) 
    print(glue("The number of cores currently using {cores}...", cores = ncores))
    # starting the main loop
    output <- foreach(i = 1:length(met_names)) %dopar% {
      print(paste0("Current at ", met_names[[i]]))
      # this section remakes the list format into the rset format to use the proper tidyposterior setting
      sum <- summary_list[[met_names[i]]]
      sample_obj <- sum %>% select(id, id2, starts_with(time))
      attr(sample_obj, "class") <- c("vfold_cv", "rset", attr(sample_obj, "class"))
      attr(sample_obj, "v") <- 5
      attr(sample_obj, "repeats") <- 100
      mod <- perf_mod(sample_obj, verbose = FALSE) # standard gaussian
      return(mod)
    }
    stopCluster(cluster) # remove cluster
  } else {
    output <- foreach(i = 1:length(met_names)) %do% {
      print(paste0("Current at ", met_names[[i]]))
      # this section remakes the list format into the rset format to use the proper tidyposterior setting
      sum <- summary_list[[met_names[i]]]
      sample_obj <- sum %>% select(id, id2, starts_with(time))
      attr(sample_obj, "class") <- c("vfold_cv", "rset", attr(sample_obj, "class"))
      attr(sample_obj, "v") <- 5
      attr(sample_obj, "repeats") <- 100
      # fitting bayesian model
      mod <- perf_mod(sample_obj, verbose = FALSE) # standard gaussian
      return(mod)
    }
  }
  names(output) <- met_names
  return(output)
}

print(opt$eval)
print(opt$timepoint)
input <- glue("output/analyses/prediction/processed/{mettype}_{eval}_by_met.rds", mettype = opt$metabtype, eval = opt$eval)
summary_list <- readRDS(file = input)
mods <- get_bayesian(summary_list = summary_list, time = opt$time, parallel = opt$parallel, ncores = opt$ncores)
saveRDS(mods, file = glue("{dir}/{time}_{met}_{eval}_bayes.rds", dir = opt$output_dir, time = opt$timepoint, eval = opt$eval, met = opt$metabtype))


