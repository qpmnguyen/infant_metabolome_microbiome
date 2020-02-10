library(tidyposterior)
library(tidyverse)
library(rsample)
library(glue)
library(rlist)

get_bayesian <- function(eval, time, type = "tar"){
  summary_list <- readRDS(file = glue("output/analyses/prediction/processed/{type}_{eval}_by_met.rds", type = type, eval = eval))
  met_names <- names(summary_list)
  output <- list()
  for (i in 1:length(met_names)){
    print(paste0("Current at ", met_names[[i]]))
    sum <- summary_list[[met_names[i]]]
    print(head(sum))
    sample_obj <- sum %>% select(id, id2, starts_with(time))
    attr(sample_obj, "class") <- c("vfold_cv", "rset", attr(sample_obj, "class"))
    attr(sample_obj, "v") <- 5
    attr(sample_obj, "repeats") <- 100
    mod <- perf_mod(sample_obj) # standard gaussian
    output[[i]] <- mod
  }
  names(output) <- met_names
  return(output)
}

eval <- c("r2", "corr", "rmse")
timepoints <- c("6W", "12M")

summary <- list('r2' = list('6W' = list(), '12M' = list()),
                "corr" = list('6W' = list(), '12M' = list()),
                "rmse" = list('6W' = list(), '12M' = list()))

for (i in 1:length(eval)){
  for (j in 1:length(timepoints)){
    print(paste("Eval", eval[i]))
    print(paste("Time", timepoints[j]))
    if (file.exists(glue("/dartfs-hpc/rc/lab/H/HoenA/Lab/QNguyen/ResultsFiles/{time}_tar_{eval}_bayes.rds", time = timepoints[j], eval = eval[i])) == FALSE){
      mods <- get_bayesian(eval = eval[i], time = timepoints[j], type = "tar")
      saveRDS(mods, glue("/dartfs-hpc/rc/lab/H/HoenA/Lab/QNguyen/ResultsFiles/{time}_tar_{eval}_bayes.rds", time = timepoints[j], eval = eval[i]))
    } else {
      mods <- readRDS(file = glue("/dartfs-hpc/rc/lab/H/HoenA/Lab/QNguyen/ResultsFiles/{time}_tar_{eval}_bayes.rds", time = timepoints[j], eval = eval[i]))
    }
    summaries <- lapply(mods, function(x){
      sample <- tidy(x)
      sum <- summary(sample)
      return(sum)
    })
    summary[[eval[i]]][[timepoints[j]]] <- summaries
  }
}

saveRDS(summary, file = "output/analyses/prediction/processed/bayesian_models_summary.rds")