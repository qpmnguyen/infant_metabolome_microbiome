library(optparse)
library(caret)
library(spls)
library(MLmetrics)

option_list <- list(
  make_option("--input", help="Input rds file of results"),
  make_option("--method", type = "character", help = "method of modelling"),
  make_option("--infolds", type = "integer", help = "Number of inner folds"),
  make_option("--outfolds", type = "integer", help = "Number of outer folds"),
  make_option("--metid", type = "integer", help = "Index of metabolite to fit"),
  make_option()
)

opt <- parse_args(OptionParser(option_list = option_list))

#' @title Fitting a function for models for single targets  
#' @description This function does nested cross validation for each of the selected models for single target only  
#' @param response Vector of responses 
#' @param predictors Matrix of predictors  
#' @param model Model to fit which includes \code{enet}, \code{rf}, \code{svm}, \code{spls}, \code{gbm}
#' @param in.folds Number of inner folds
#' @param out.folds Number of outer folds  
#' @param control This controls the parameter tuning portion of all models 
modelfit.fn <- function(response, predictors, model, in.folds, out.folds, folds = NULL, control = NULL){
  # certain model translations to work between common knowledge and caret specific syntax
  model_translate = list(enet = "glmnet", rf = "rf", svm = "svmRadial")
  if (is.null(control) == T){
    print("Using default parameters for training...")
    control <- list(method = "cv", search = "grid", verboseIter = T, tuneLength = 30)
    print(control)
  }  
  # setting control parameters 
  if (is.null(control) == F & length(control) < 4){
    stop("Control list not completely specified, please fill everything in if using custom control parameters")
  }
  # initialize starting conditions 
  rmse <- c()
  r2 <- c()
  rrse <- c()
  corr <- c()
  # custom folds or random folds    
  if (is.null(folds) == T){
    folds <- createFolds(1:nrow(predictors),out.folds)  
  }
  ctrl <- trainControl(method = control$method, number = in.folds, search = control$search, verboseIter = control$verboseIter)
  for(i in 1:length(folds)){
      met.train <- response[-test]
      met.test <- response[test]
    tax.train <- predictors[-test,]
    tax.test <- predictors[test,]
    # fitting models 
    if (model %in% c("enet", "svm", "rf")){
      fit <- train(x = tax.train, y = met.train, trControl = ctrl, method = model_translate[[model]], 
                   metric = "RMSE", tuneLength = control$tuneLength)
      predictions <- predict(fit, tax.test)
    } else if (model == "spls"){
      mod.cv <- cv.spls(tax.train, met.train, K = c(1:10), eta = seq(0.1,0.9,0.1),
                        scale.x = F, scale.y = F, kappa = 0.5, plot.it = F)
      fit <- spls(tax.train, met.train, K = mod.cv$K.opt, eta = mod.cv$eta.opt,
                  scale.x = F, scale.y  = F, kappa = 0.5) 
      predictions <- predict.spls(fit, newx=tax.test, type = "fit")
    } 
    # parsing outputs
    # metabolite specific results from these models  
    rmse[i] <- MLmetrics::RMSE(y_pred = as.vector(predictions), y_true = as.vector(met.test))
    r2[i] <- MLmetrics::R2_Score(y_pred = as.vector(predictions), y_true = as.vector(met.test))
    rrse[i] <- MLmetrics::RRSE(y_pred = as.vector(predictions), y_true = as.vector(met.test))
    corr[i] <- cor(x = as.vector(predictions), y = as.vector(met.test), method = "pearson")
    print(paste("Finish fold",i,"..."))
  }
  result <- list(rmse = rmse, r2 = r2, rrse = rrse, corr = corr)
  return(result)
}

data <- readRDS(file = opt$input)
tax <- data$tax
met <- data$met 

saveRDS(file = opt$output, )

