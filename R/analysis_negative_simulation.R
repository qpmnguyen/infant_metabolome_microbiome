# Script to perform 500 permutations
library(tidyverse)
library(MLmetrics)
library(caret)
library(glue)
library(phyloseq)
library(purrr)
library(parallel)
library(optparse)
library(doParallel)
library(foreach)
library(picante)
library(rsample)
# defining functions similar to positive simulations 
option_list <- list(
  make_option("--input", help="Input rds file of results"),
  make_option("--nperms", help="Number of permutations"),
  make_option("--parallel", help = "Check if parallel processing is required"),
  make_option("--ncores", help = "Number of cores used if parallel is required", default = NULL)
)

opt <- parse_args(OptionParser(option_list = option_list))

#' Function to fit models to data 
#' @param data The data to be fit
#' @param model types of model, can be spls, enet, svm_rbf, xgboost, rf
#' @param folds number of inner folds for cross validation 
#' @return a fitted model 
model_fit <- function(data, model, nfolds = 10){
  ctrl <- caret::trainControl(method = "cv", 
                              number = nfolds, 
                              search = "grid", 
                              verboseIter = F,
                              allowParallel = T)
  model_translate <- list(enet = "glmnet", rf = "rf", svm_rbf = "svmRadial", gbm = "gbm", xgboost = "xgbTree",
                          spls = "spls")
  if (model == "enet"){
    tune_length <- 10
    tune_grid <- NULL
  } else if (model == "svm_rbf"){
    tune_length <- NULL
    tune_grid <- expand.grid(
      sigma = sqrt(1/(2*2^seq(-5,15,2))),
      C = 2^seq(-15,3,2)
    )
  } else if (model == "xgboost"){
    tune_length <- NULL
    tune_grid <- expand.grid(
      nrounds = seq(from=200, to=1000, by = 100), # number of trees
      eta = c(0.05, 0.1, 0.2, 0.3), # learning rate 
      gamma = 0, # min_split_loss - thresholding for tree pruning
      colsample_bytree = c(0.4, 0.6, 0.8, 1.0), # column subsample 
      subsample = c(0.5, 0.75, 1.0),# row subsample
      max_depth = c(4,6,8,10), # max tree depth,
      min_child_weight = 3
    )
  } else if (model == "rf"){
    tune_length <- NULL
    tune_grid <- data.frame(mtry = sqrt(ncol(data)-1))
  } else if (model == "spls"){
    tune_length <- NULL
    tune_grid <- expand.grid(kappa = 0.5,
                             K = seq(1,10,1),
                             eta = seq(0.1,0.9,0.1))
  }
  fit <- caret::train(met ~., data = data, trControl = ctrl, method = model_translate[[model]],
                      metric = "RMSE", tuneLength = tune_length, tuneGrid = tune_grid, preProcess = c("center", "scale"))
  return(fit)
}

#' Function to perform holdout assessment
#' @param split rsplit type object from initial_split
#' @param model model to be fitted, can be rf, spls, enet, svm_rbf, xgboost
#' @param eval evaluation metric 
holdout_assessment <- function(split, model, eval){
  fit <- model_fit(training(split), model = model)
  print("At evaluation step...")
  predictions <- predict(fit, testing(split)[,-ncol(testing(split))])
  if (eval == "r2"){
    metric <- R2_Score(y_pred = testing(split)$met, y_true = predictions)
  } else {
    metric <- cor(x = testing(split)$met, y = predictions, method = "spearman")
  }
  return(metric)
}

#' Function to perform evaluation 
#' @param data The data to be split on 
#' @return a list of two elements split by two models 
train_test_eval <- function(data){
  split <- initial_split(data, prop = 0.8)
  print(split)
  output <- list(r2 = rep(0,4), corr = rep(0,4))
  for (j in c("r2", "corr")){
    vec <- c()
    for (i in c("enet", "rf", "spls", "svm_rbf")){
      vec <- c(vec, holdout_assessment(split, model = i, eval = j))
    }
    output[[j]] <- vec
  }
  return(output)
}

data <- readRDS(file = opt$input)
tax <- as(otu_table(data), "matrix")
met <- as(sample_data(data), "matrix")

if (opt$parallel == T){
  if (is.null(opt$ncores) == TRUE){ # if no cores are specified 
    ncores <- parallel::detectCores() - 1        
  } else {
    ncores <- opt$ncores
  }
  print(paste0("Using parallel with ncores = ", ncores))
  cluster <- makeForkCluster(ncores)
  registerDoParallel(cluster) # register
  summary <- foreach(i = 1:opt$nperms) %dopar%{
    perm_tax <- randomizeMatrix(tax, null.model = "richness", iterations = 1000)
    rand_met <- met[,sample(1:ncol(met), size = 1, replace = F)]
    perm_data <- as.data.frame(cbind(perm_tax, met = rand_met))
    return(train_test_eval(perm_data))
  }
  stopCluster(cluster) # remove cluster
} else {
  print("Processing non-parallel")
  summary <- foreach(i = 1:opt$nperms) %do%{
    perm_tax <- randomizeMatrix(tax, null.model = "richness", iterations = 1000)
    rand_met <- met[,sample(1:ncol(met), size = 1, replace = F)]
    perm_data <- as.data.frame(cbind(perm_tax, met = rand_met))
    return(train_test_eval(perm_data))
  }
}

print("Finished simulations...")
split <- strsplit(opt$input,".rds")[[1]]
time <- strsplit(split,"_")[[1]][2]
saveRDS(summary, file = glue("output/analyses/simulations/{time}_negative_simulation_results.rds"))

