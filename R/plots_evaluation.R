library(ggpubr)
library(ggplot2)
library(cowplot)
library(glue)
library(phyloseq)
library(viridis)
library(ggplotify)
library(optparse)
library(MLmetrics)
library(tidyposterior)
library(rsample)

filenames <- list.files("output/analyses/prediction/")

# TODO add support for linux machine 
if(Sys.info()["sysname"] == "Darwin"){
  if (file.exists("/Volumes/rc-1/Lab/QNguyen/ResultsFiles/data/processed_6W_tar_phyloseq_obj.rds")){
    data <- readRDS(file = "/Volumes/rc-1/Lab/QNguyen/ResultsFiles/data/processed_6W_tar_phyloseq_obj.rds")
    met_names <- colnames(sample_data(data))
  } else if (file.exists("/Volumes/rc/Lab/QNguyen/ResultsFiles/data/processed_6W_tar_phyloseq_obj.rds")){
    data <- readRDS(file = "/Volumes/rc/Lab/QNguyen/ResultsFiles/data/processed_6W_tar_phyloseq_obj.rds")
    met_names <- colnames(sample_data(data))
  } else {
    stop("Can't find mounted drives")
  }
}

# A function to get attribute types 
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
      res <- readRDS(file = glue("output/analyses/prediction/{filename}", filename = names[j]))
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

folds <- vfold_cv(data = data.frame(rnorm(1000),100,10),v = 5, repeats = 100)
folds
test
full_join(folds, test, by = c("id","id2"))
tar_sum <- get_summary(filenames = filenames, type = "tar", met_names = met_names, eval = "r2")
test <- tar_sum[[1]]
contrast_models(accuracy_dist)
mod <- perf_mod(test)

