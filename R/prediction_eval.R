# Script to process prediction results -- allowing interfacing to snakemake as well as custom options 
# Quang Nguyen
# Updated 01/13/2020

library(optparse)
library(reshape2)
library(caret)
library(MLmetrics)

option_list <- list(
  make_option("--model", help = "Selected model for evaluation", default = NULL),
  make_option("--metab_type", help = "Type of metabolomics data", default = NULL),
  make_option("--time", help = "Data time point to extract", default = NULL),
  make_option("--tax_type", help = "Data type for taxonomic data, so far only '16S' is supported", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$model) == T & exists(snakemake) == F){
  stop("Need to specify arguments if snakemake object does not exist")
} else if (is.null(opt$model) != T & exists(snakemake) == F){
  tax_type <- opt$tax_type
  time <- opt$time
  metab_type <- opt$metab_type
  model <- opt$model 
} else if (is.null(opt$model) == T & exists(snakemake) == T){
  tax_type <- snakemake@input[["tax_type"]]
  time <- snakemake@input[["time"]]
  metab_type <- snakemake@input[["metab_type"]]
  model <- snakemake@input[["model"]]
}

names <- paste0(tax_type, "_", time, "_", metab_type, "_", model)
files <- paste0("snakemake_output/analyses/prediction/",list.files("snakemake_output/analyses/prediction", pattern = names))
results <- lapply(files, readRDS)
correlation <- c()
for (i in seq_len(length(results))){
  per_rep <- sapply(results[[i]], function(x){
    cor <- c()
    for (j in seq_len(length(x))){
      cor[j] <- cor(x[[j]][,1],x[[j]][,2])
    }
    return(mean(cor))
  })
  correlation <- rbind(correlation, per_rep)
}
correlation <- t(correlation)

r2 <- c()
for (i in seq_len(length(results))){
  per_rep <- sapply(results[[i]], function(x){
    r2_sc <- c()
    for (j in seq_len(length(x))){
      r2_sc[j] <- MLmetrics::R2_Score(y_pred = x[[j]][,1],y_true = x[[j]][,2])
    }
    return(mean(cor))
  })
  r2 <- rbind(r2, per_rep)
}
r2 <- t(r2)

met_names <- colnames(readRDS(file = paste0("data/processed/", tax_type, "_", time, "_", metab_type, "_processed_prediction.rds"))$met)
names_order <- as.vector(sapply(files, function(x){
  as.numeric(strsplit(strsplit(x, ".rds")[[1]], "_")[[1]][6])
}))

names(results) <- met_names[names_order]
colnames(r2) <- met_names[names_order]
colnames(correlation) <- met_names[names_order]
output <- list(raw = results, correlation = correlation, r2 = r2)

path <- paste0("snakemake_output/analyses/prediction_evaluation/", tax_type, "_", time, "_", metab_type, "_", model, ".rds")
saveRDS(output, file = path)

# TODO test this script out  