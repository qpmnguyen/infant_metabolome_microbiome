library(optparse)
library(caret)
library(spls)
library(MLmetrics)
library(phyloseq)
# Install 'randomforest', 'spls', 'gbm', 'xgboost', 'kernlab', 'glmnet', 'plyr'

option_list <- list(
  make_option("--input", help="Input rds file of results"),
  make_option("--output_label", help="Output file directory and name"),
  make_option("--method", type = "character", help = "method of modelling"),
  make_option("--infolds", type = "integer", help = "Number of inner folds"),
  make_option("--outfolds", type = "integer", help = "Number of outer folds"),
  make_option("--nrep", type = "integer", help = "Number of repetitions"),
  make_option("--metid", type = "integer", help = "Index of metabolite to fit"),
  make_option("--preprocess", type = "logical", help = "Whether to center and scale the x matrix"),
  make_option("--metabtype", type = "character", help = "Type of metabolite data to inform back transformations")
)

opt <- parse_args(OptionParser(option_list = option_list))

#' @title Fitting a function for models for single targets  
#' @description This function does nested cross validation for each of the selected models for single target only  
#' @param response Vector of responses 
#' @param predictors Matrix of predictors  
#' @param model Model to fit which includes \code{enet}, \code{rf}, \code{svm_rbf}, \code{spls}, \code{gbm}, \code{xgboost}
#' @param in.folds Number of inner folds
#' @param out.folds Number of outer folds  
#' @param parallel Allows for parallel processing with caret for internal cross validation 
modelfit.fn <- function(response, predictors, model, in.folds, out.folds, folds = NULL,
                        parallel = F, n_cores = NULL){
  # certain model translations to work between common knowledge and caret specific syntax
  model_translate <- list(enet = "glmnet", rf = "rf", svm_rbf = "svmRadial", gbm = "gbm", xgboost = "xgbTree")
  # initialize prediction list
  pred_matrix <- list()
  # custom folds or random folds    
  if (is.null(folds) == T){
    folds <- caret::createFolds(1:nrow(predictors),out.folds)  
  }
  # setting parameters for models 
  ctrl <- caret::trainControl(method = "cv", 
                              number = in.folds, 
                              search = "grid", 
                              verboseIter = T,
                              allowParallel = T)
  if (model == "enet"){
    tune_length <- 100
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
    tune_grid <- data.frame(mtry = sqrt(ncol(predictors)))
  } else if (model == "spls"){
    tune_length <- NULL
    tune_grid <- expand.grid(kappa = 0.5,
                             K = seq(1,10,1),
                             eta = seq(0.1,0.9,0.1))
  }
  for(i in 1:length(folds)){
    test <- folds[[i]]
    met.train <- response[-test]
    met.test <- response[test]
    tax.train <- predictors[-test,]
    tax.test <- predictors[test,]
    # fitting models 
    if (parallel == T){
      if (is.null(n_cores) == T){
        stop("If parallel need to determine number of cores")
      }
      cl <- parallel::makeForkCluster(n_cores)
      doParallel::registerDoParallel(cl)
    }
    fit <- caret::train(x = tax.train, y = met.train, trControl = ctrl, method = model_translate[[model]],
                        metric = "RMSE", tuneLength = tune_length, tuneGrid = tune_grid)
    predictions <- predict(fit, tax.test)
    if (exists(cl) == T){
      parallel::stopCluster(cl)
    }
    # parsing outputs
    # back transformations 
    if (opt$metabtype == "tar"){
      predictions <- exp(predictions) - 1
      met.test <- exp(met.test) - 1
    } else if (opt$metabtype == "untar"){
      predictions <- sin(predictions^2)
      met.test <- sin(met.test^2)
    }
    name <- paste0("outer_fold_", i)
    pred_matrix[[name]] <- cbind(predictions, met.test)
    print(paste("Finish fold",i,"..."))
  }
  return(pred_matrix)
}

data <- readRDS(file = opt$input)
tax <- data$tax
met <- data$met[,opt$metid] 

# check if preprocessing
if (opt$preprocess == T){
  tax <- as.matrix(scale(tax, center = T, scale = T))
}

# running results as repeated CV
result <- list()
for (j in 1:opt$nrep){
  rep <- modelfit.fn(response = met, predictors = tax, model = opt$method, 
                        in.folds = opt$infolds, out.folds = opt$outfolds)
  name <- paste0("rep_", j)
  result[[name]] <- rep
}

path <- paste0("snakemake_output/analyses/prediction/", opt$output_label, "_", opt$method, "_", opt$metid, ".rds")
saveRDS(result, file = path)