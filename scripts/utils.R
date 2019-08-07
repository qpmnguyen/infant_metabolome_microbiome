#' This function processes the taxonomic key matrix 
#' @param key The key matrix  
#' @param type Method to transform unknown taxonomic assignment. NA turns all unknown into NAs where as
#'    'Unknown' adds unknown as a suffix to all non-assignable taxonomic ranks. 

key_processing <- function(key, type = "Unknown"){
  strings <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__" )
  if (type == "NA"){
    for (i in 1:length(strings)){
      key[which(key == strings[i])] <- NA
    }
  }
  
  if (type == "Unknown"){
    for (i in 1:length(strings)){
      key[which(key == strings[i])] <- paste0(strings[i], "Unknown")
    }
    for (i in 1:length(colnames(key))){
      name <- paste0(tolower(substr(colnames(key)[i],1,1)),"__")
      index <- as.vector(which(is.na(key[,colnames(key)[i]])))
      key[,colnames(key)[i]][index] <- paste0(name,"Unknown")
    }
  }
  return(key)
}

# inverse transformations 
inv_transform <- function(met, type){
  if (type == "concentration"){
    return(exp(met) - 1)
  } else if (type == "binned"){
    return(sin(met)^2)
  }
}



## refactoring modelfitting functions for single targets   

#' @title Fitting a function for models for single targets  
#' @description This function does nested cross validation for each of the selected models for single target only  
#' @param response Vector of responses 
#' @param predictors Matrix of predictors  
#' @param model Model to fit which includes \code{enet}, \code{rf}, \code{svm}, \code{spls}, \code{SICS}
#' @param in.folds Number of inner folds
#' @param out.folds Number of outer folds  
#' @param distance Distance matrix of OTU based on the tree for SICS method 
#' @param control This controls the parameter tuning portion of all models 
modelfit.fn <- function(response, predictors, model, in.folds, out.folds, resp_type, distance = NULL, folds = NULL, control = NULL){
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
  # custom folds or random folds    
  if (is.null(folds) == T){
    folds <- createFolds(1:nrow(predictors),out.folds)  
  }
  
  ctrl <- trainControl(method = control$method, number = in.folds, search = control$search, verboseIter = control$verboseIter)
  
  # initialize result 
  result <- list()
  # nested cross-validation  
  for(i in 1:length(folds)){
    test <- folds[[i]]
    tax.train <- predictors[-test,]
    tax.test <- predictors[test,]
    if (resp_type == "concentration"){
      response_norm <- log(response + 1)
    } else if (resp_type == "binned"){
      response_norm <- asin(sqrt(response))
    }
    met.train <- response_norm[-test] # transformed output 
    met.test <- response[test] # original data 
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
    } else if (model == "SICS"){
      if (is.null(distance) == T){
        stop("SICS requires a distance matrix")
      }
      fit <- SICS(X = tax.train, Y = met.train, D = distance, family = "gaussian", nfolds = in.folds, 
                  lambda2 = c(0,2^(seq(-5, 5, length = control$tuneLength))))
      predictions <- predict(fit, tax.test, family = "gaussian")
    } else if (model == "menet"){
      cv.obj <- cv.glmnet(x = tax.train, y = met.train, nfolds = in.folds, family = "mgaussian")
      fit <- glmnet(x=tax.train, y=met.train, family = "mgaussian", alpha = 0.8, 
                    lambda = cv.obj$lambda.min, standardize = FALSE)
      predictions <- predict(fit, newx = tax.test, s = cv.obj$lambda.min, type = "link")[,,1]
    } else if (model == "mspls"){
      mod.cv <- cv.spls(tax.train, met.train, K = c(1:10), eta = seq(0.1,0.9,0.1),
                        scale.x = F, scale.y = F, kappa = 0.5, plot.it = F)
      
      fit <- spls(tax.train, met.train, K = mod.cv$K.opt, eta = mod.cv$eta.opt,
                  scale.x = F, scale.y  = F, kappa = 0.5) 
      predictions <- predict.spls(fit, newx=tax.test, type = "fit")
    } else if (model == "mrf"){
      predictions <- build_forest_predict(trainX = tax.train, trainY = met.train, 
                                          n_tree = 100,  m_feature = 5, min_leaf = 5,
                                          testX = tax.test)
    }
    predictions <- inv_transform(predictions, type = resp_type)
    pair_res <- cbind(predictions, met.test)
    result[[i]] <- pair_res
  }
  return(result)
}


