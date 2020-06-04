# Script to perform 100 bootstrapped evaluations of each positive simulation setting
library(tidyverse)
library(MLmetrics)
library(caret)
library(rsample)
library(glue)
library(phyloseq)
library(purrr)
library(parallel)
library(optparse)
library(doParallel)
library(foreach)

option_list <- list(
  make_option("--input", help="Input rds file of results"),
  make_option("--parallel", help = "Check if parallel processing is required"),
  make_option("--ncores", help = "Number of cores used if parallel is required", default = NULL)
)

opt <- parse_args(OptionParser(option_list = option_list))

data <- readRDS(file = opt$input)
tax <- as(otu_table(data),"matrix")
met <- as(sample_data(data), "matrix")

#' @title Function to generate simulated data sets from bootstrapped resamples of the original data
#' @param split a split of the rsplit object
#' @param prob Controls sparsity - the probability a predictor is 0
#' @param snr Signal to noise ratio
#' @return data frame containing a met column as the outcome and everything else as predictors 
simulate_y <- function(split, prob, snr){
  data <- analysis(split)
  y <- rep(0, nrow(data))
  while (sd(y) == 0){ # remake y until the standard deviation is larger than 0 
    beta <- as.vector(rnorm(ncol(data)) %*% diag(rbinom(ncol(data), size = 1, prob = prob)))
    beta_0 <- 6/sqrt(10)
    y_mean <- beta_0 + as.matrix(data) %*% beta
    y <- y_mean + rnorm(nrow(data), mean = 0, sd = sd(y_mean)* (1/snr))
  }
  output <- as.data.frame(cbind(data, y))
  colnames(output)[ncol(output)] <- "met"
  return(output)
}

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
    metric <- rsq_trad(data = data.frame(truth = testing(split)$met, estimate = predictions), 
                       truth = truth, estimate = estimate)
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

grid <- expand.grid(snr = c(0.5, 0.7, 3, 5),
                    spar = c(0.05, 0.1, 0.5, 0.95))
if (opt$parallel == T){
  if (is.null(opt$ncores) == TRUE){ # if no cores are specified 
    ncores <- parallel::detectCores() - 1        
  } else {
    ncores <- opt$ncores
  }
  print(paste0("Using parallel with ncores = ", ncores))
  cluster <- makeForkCluster(ncores)
  registerDoParallel(cluster) # register
  summary <- foreach(i=1:nrow(grid)) %dopar% {
    bts <- bootstraps(tax, times = 100)
    data <- map(bts$splits, simulate_y, prob = grid$spar[i], snr = grid$snr[i])
    results <- map(data, train_test_eval)
    return(results)
  }
  stopCluster(cluster) # remove cluster
} else {
  print("Processing non-parallel")
  summary <- foreach(i=1:nrow(grid)) %do% {
    bts <- bootstraps(tax, times = 100)
    data <- map(bts$splits, simulate_y, prob = grid$spar[i], snr = grid$snr[i])
    results <- map(data, train_test_eval)
    return(results)
  }
}
print("Finished simulations...")
split <- strsplit(opt$input,".rds")[[1]]
time <- strsplit(split,"_")[[1]][2]
saveRDS(list(grid = grid, summary = summary), file = glue("output/analyses/simulations/{time}_simulation_results.rds"))
