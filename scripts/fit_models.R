library(foreach)
library(doParallel)
library(caret)
library(MLmetrics)
library(optparse)

source("./utils.R")

option_list <- list(
  make_option("--input", help="Input rds file of results"),
  make_option("--method", type = "character", help = "method of modelling"),
  make_option("--infolds", type = "integer", help = "Number of inner folds"),
  make_option("--outfolds", type = "integer", help = "Number of outer folds"),
  make_option("--label", type = "character", help = "If not null simulations input this to have proper names", 
              default = NULL),
  make_option("--datatype", type = "character", help = "NMR data type which determines different transformations, can be 'concentrations' or 'binned'"),
  make_option("--distance", help = "Distance file for SICS", default = NULL),
  make_option("--ncores", help = "Number of cores")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (opt$method %in% c("spls", "rf", "enet", "SICS", "svm") == FALSE){ stop("Not a supported model type")} # get models  

registerDoParallel(opt$ncores)

if (!is.null(opt$distance)){
  dist <- readRDS(file = opt$distance)
} else {
  dist <- NULL
}

data <- readRDS(file = opt$input)
tax <- data$tax
met <- data$met
in.folds <- opt$infolds
out.folds <- opt$outfolds

tax <- scale(tax, center = T, scale = T)

if(is.null(ncol(met) == T)){ # if met is one vector 
  results <- modelfit.fn(response = met, predictor = tax, model = opt$method, resp_type = opt$datatype, distance = dist, in.folds = in.folds, out.folds = out.folds)
} else {
  results <- foreach(i = 1:ncol(met)) %dopar% {
    mod <- modelfit.fn(response = met[,i], predictor = tax, model = opt$method, resp_type = opt$datatype, distance = dist, in.folds = opt$infolds, out.folds = opt$outfolds)
    mod
  }
}

filename <- paste0("./",opt$label,"_", opt$method, "_", opt$datatype, ".rds")
saveRDS(result, file = filename)






