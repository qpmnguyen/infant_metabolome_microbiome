library(optparse)
library(caret)
library(spls)
library(MLmetrics)
library(phyloseq)
# Install 'randomforest', 'spls', 'gbm', 'xgboost', 'kernlab'

option_list <- list(
  make_option("--input", help="Input rds file of results"),
  make_option("--output_label", help="Output file directory and name"),
  make_option("--method", type = "character", help = "method of modelling"),
  make_option("--infolds", type = "integer", help = "Number of inner folds"),
  make_option("--outfolds", type = "integer", help = "Number of outer folds"),
  make_option("--metid", type = "integer", help = "Index of metabolite to fit"),
  make_option("--preprocess", type = "logical", help = "Whether to center and scale the x matrix"),
  make_option("--metabtype", type = "character", help = "Type of metabolite data to inform back transformations")
)

opt <- parse_args(OptionParser(option_list = option_list))

#' @title Fitting a function for models for single targets  
#' @description This function does nested cross validation for each of the selected models for single target only  
#' @param response Vector of responses 
#' @param predictors Matrix of predictors  
#' @param model Model to fit which includes \code{enet}, \code{rf}, \code{svm}, \code{spls}, \code{gbm}, \code{xgboost}
#' @param in.folds Number of inner folds
#' @param out.folds Number of outer folds  
#' @param control This controls the parameter tuning portion of all models 
modelfit.fn <- function(response, predictors, model, in.folds, out.folds, folds = NULL, control = NULL){
  # certain model translations to work between common knowledge and caret specific syntax
  model_translate <- list(enet = "glmnet", rf = "rf", svm = "svmRadial", gbm = "gbm", xgboost = "xgbTree")
  if (is.null(control) == T){
    print("Using default parameters for training...")
    control <- list(method = "cv", search = "grid", verboseIter = T, tuneLength = 50)
    print(control)
  }  
  # setting control parameters 
  if (is.null(control) == F & length(control) < 4){
    stop("Control list not completely specified, please fill everything in if using custom control parameters")
  }
  # initialize prediction list
  pred_matrix <- list()
  # custom folds or random folds    
  if (is.null(folds) == T){
    folds <- caret::createFolds(1:nrow(predictors),out.folds)  
  }
  ctrl <- caret::trainControl(method = control$method, number = in.folds, search = control$search, verboseIter = control$verboseIter)
  for(i in 1:length(folds)){
    test <- folds[[i]]
    met.train <- response[-test]
    met.test <- response[test]
    tax.train <- predictors[-test,]
    tax.test <- predictors[test,]
    # fitting models 
    if (model %in% c("enet", "svm", "rf", "xgboost")){
      fit <- train(x = tax.train, y = met.train, trControl = ctrl, method = model_translate[[model]], 
                   metric = "RMSE", tuneLength = control$tuneLength)
      predictions <- predict(fit, tax.test)
    } else if (model == "gbm"){
      fit <- train(x = tax.train, y = met.train, trControl = ctrl, method = model_translate[[model]], metric = "RMSE",
                   tuneLength = control$tuneLength, bag.fraction = 0.75)
      predictions <- predict(fit, tax.test)
    } else if (model == "spls"){
      mod.cv <- cv.spls(tax.train, met.train, K = c(1:10), eta = seq(0.1,0.9,0.1),
                        scale.x = F, scale.y = F, kappa = 0.5, plot.it = F)
      fit <- spls(tax.train, met.train, K = mod.cv$K.opt, eta = mod.cv$eta.opt,
                  scale.x = F, scale.y  = F, kappa = 0.5) 
      predictions <- spls::predict.spls(fit, newx=tax.test, type = "fit")
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
    pred_matrix[[i]] <- cbind(predictions, met.test)
    print(paste("Finish fold",i,"..."))
  }
  return(pred_matrix)
}

data <- readRDS(file = opt$input)
tax <- data$tax
met <- data$met[,opt$metid] 

if (opt$preprocess == T){
  tax <- as.matrix(scale(tax, center = T, scale = T))
}

result <- modelfit.fn(response = met, predictors = tax, model = opt$method, 
                       in.folds = opt$infolds, out.folds = opt$outfolds)
path <- paste0("snakemake_output/analyses/prediction/", opt$output_label, "_", opt$method, "_", opt$metid, ".rds")
saveRDS(result, file = path)


