# Script to cobble together all evaluation metrics 
library(glue)
library(phyloseq)
library(optparse)
library(MLmetrics)
library(purrr)
library(tidyverse)


option_list <- list(
  make_option("--input", help="Input rds file of metnames"),
  make_option("--parallel", help = "Check if parallel processing is required"),
  make_option("--metabtype", help = "Tar or untar"),
  make_option("--ncores", help = "Number of cores used if parallel is required", default = NULL)
)

opt <- parse_args(OptionParser(option_list = option_list))

filenames <- list.files(glue("output/analyses/prediction/raw/{metabtype}/", metabtype = opt$metabtype))

#' @title Function to get the names of metabolites 
#' @param input the data frame of the appropriate type when using opt$input
get_names <- function(input){
  input_dat <- readRDS(file = opt$input)
  return(as.vector(colnames(sample_data(input_dat))))
}

#' @title A function to get attribute types 
#' @param filename Name of the file after using list.files function
#' @param attr name of attribute of interest
get_attributes <- function(filename, attr){
  split <- strsplit(filename, "_")[[1]]
  time <- split[1]
  met_type <- split[2]
  model <- split[3]
  met_id <- strsplit(tail(split,1),".", fixed = T)[[1]][1]
  if (attr == "time"){
    return(time)
  } else if (attr == "met_type"){
    return(met_type)
  } else if (attr == "model"){
    return(model)
  } else if (attr == "met_id"){
    return(met_id)
  }
}

#' @title Getting rsquared score 
#' @param result_list list of results cobbled together 
#' @description The results_list is a list of each fold within each repeat of the prediction results 
get_r2_score <- function(result_list){
  score <- map_df(result_list, function(.x){
    map_df(.x, function(.y){
      R2_Score(y_pred = .y[,1], y_true = .y[,2])
    })
  })
  return(score)
}

get_rmse_score <- function(result_list){
  score <- map_df(result_list, function(.x){
    map_df(.x, function(.y){
      RMSE(y_pred = .y[,1], y_true = .y[,2])
    })
  })
  return(score)
}

get_corr_score <- function(result_list){
  score <- map_df(result_list, function(.x){
    map_df(.x, function(.y){
      cor(x = .y[,1], y = .y[,2], method = "spearman")
    })
  })
  score[is.na(score)] <- 0 # add 0 to all NA values 
  return(score)
}

get_summary <- function(filenames, type, met_names, eval){
  ids <- sapply(filenames, function(x){
    str <- strsplit(tail(strsplit(x, "_")[[1]],1),".",fixed = T)[[1]][1]
    return(as.numeric(str))
  })
  metab_type <- sapply(filenames, function(x){
    str <- strsplit(x, "_")[[1]][2]
    return(str)
  })
  summary <- list()
  for (i in unique(ids)){
    names <- filenames[which(ids == i & metab_type == type)]
    output <- expand.grid(id = glue("Repeat{i}", i = sprintf("%03d",1:100)), id2 = glue("Fold{j}", j = 1:5), stringsAsFactors = F)
    for (j in 1:length(names)){
      res <- readRDS(file = glue("output/analyses/prediction/raw/{type}/{filename}", filename = names[j]))
      if (eval == "r2"){
        score <- get_r2_score(res)  
      } else if (eval == "corr") {
        score <- get_corr_score(res)
      } else if (eval == "rmse"){
        score <- get_rmse_score(res)
      } else {
        stop("Not considering that evaluation metric")
      }
      score <- score %>% as.data.frame()
      rownames(score) <- glue("Repeat{i}", i = sprintf("%03d",1:100))
      colnames(score) <- glue("Fold{j}", j = 1:5)
      score <- score %>% rownames_to_column() %>% rename(id = rowname) %>% 
        pivot_longer(-id, names_to = c("id2"), values_to = glue("{time}_{model}", time = get_attributes(names[j], "time"), 
                                                                model = get_attributes(names[j],"model")))
      output <- full_join(output, score, by = c("id", "id2"))
    }
    summary[[i]] <- output
  }
  names(summary) <- met_names
  return(summary)
}

met_names <- get_names(opt$input)
for (j in c("corr", "r2")){
  print(glue("Eval {metric}", metric = j))
  print(glue("Met type {met}", met = opt$metabtype))
  summary <- get_summary(filenames = filenames, type = opt$metabtype, met_names = met_names, eval = j)
  saveRDS(summary, file = glue("output/analyses/prediction/processed/{type}_{eval}_by_met.rds", type = opt$metabtype, eval = j))
}




