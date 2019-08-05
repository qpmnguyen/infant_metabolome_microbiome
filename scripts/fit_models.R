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
  make_option("--processed", type = "logical", help = "Indicate whether data is already normalized"),
  make_option("--distance", help = "Distance file for glmmTree", default = NULL)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (opt$method %in% c("spls", "rf", "enet", "SICS", "svm") == FALSE){ stop("Not a supported model type")} # get models  

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

if (opt$procesesed == FALSE){
  tax[tax == 0] <- 1
  tax <- unclass(acomp(tax))
  tax <- unclass(clr(tax))
} 

tax <- scale(tax, center = T, scale = T)
met <- scale(met, center = T, scale = T)


results <- foreach(i = 1:ncol(met)){
  
}






