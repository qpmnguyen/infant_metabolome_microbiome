library(foreach)
library(doParallel)
library(caret)
library(MLmetrics)
library(optparse)
library(parallel)
source("./utils.R")

option_list <- list(
  make_option("--input", help="Input rds file of results"),
  make_option("--method", type = "character", help = "method of modelling"),
  make_option("--infolds", type = "integer", help = "Number of inner folds"),
  make_option("--outfolds", type = "integer", help = "Number of outer folds"),
  make_option("--label", type = "character", help = "If not null simulations input this to have proper names", 
              default = NULL),
  make_option("--datatype", type = "character", help = "NMR data type which determines different transformations, can be 'concentrations' or 'binned'"),
  make_option("--ncores", help = "Number of cores")
)



opt <- parse_args(OptionParser(option_list = option_list))
# check for consistency of methods 
if (opt$method %in% c("spls", "rf", "enet", "SICS", "svm") == FALSE){ stop("Not a supported model type")} # get models

print(paste("The number of cores:", opt$ncores))

# read data files and assign variables 
data <- readRDS(file = opt$input)
tax <- data$tax
met <- data$met
in.folds <- opt$infolds
out.folds <- opt$outfolds

# convert to singular matrix form if there is only one column
if (is.null(ncol(met)) == T){
  met <- as.matrix(met)
}
# scale and center taxonomic tables  
tax <- scale(tax, center = T, scale = T)

print("Starting loop...")
# main dopar loop
registerDoParallel(cores = opt$ncores)
output <- foreach(i = 1:ncol(met)) %dopar% {
    mod <- modelfit.fn(response = met[,i], predictor = tax, model = opt$method, 
                       resp_type = opt$datatype, distance = NULL, 
                       in.folds = in.folds, out.folds = out.folds)
    print(paste("Loop", i, "Finished..."))
    mod
    
}

print("End loop...")
# save files  
filename <- paste0(opt$label,"_", opt$method, "_", opt$datatype, ".rds")
saveRDS(output, file = filename)





