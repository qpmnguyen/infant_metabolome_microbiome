library(glue)
library(tidyposterior)
library(optparse)

option_list <- list(
  make_option("--dir", help = "Directory where all the bayesian models are kept -- don't include backslash!")
)

opt <- parse_args(OptionParser(option_list = option_list))

eval <- c("r2", "corr")
timepoints <- c("6W", "12M")

for (i in 1:length(eval)){
  for (j in 1:length(timepoints)){
    print(paste("Eval", eval[i]))
    print(paste("Time", timepoints[j]))
    if (file.exists(glue("{dir}/{time}_tar_{eval}_bayes.rds", dir = opt$dir, time = timepoints[j], eval = eval[i])) == FALSE){
      print(glue("File does not exist for {time} and {eval}", time = timepoints[j], eval = eval[i]))
      next
    } else {
      print("Existing Bayes file loaded")
      mods <- readRDS(file = glue("{dir}/{time}_tar_{eval}_bayes.rds", dir = opt$dir, time = timepoints[j], eval = eval[i]))
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