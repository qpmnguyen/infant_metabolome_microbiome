library(optparse)
library(reshape2)
library(caret)
library(MLmetrics)
options_list <- list(
  make_option("--model", help = "Selected model for evaluation"),
  make_option("--metab_type", help = "Type of metabolomics data"),
  make_option("--time", help = "Data time point to extract"),
  make_option("--tax_type", help = "Data type for taxonomic data, so far only '16S' is supported")
)

opt <- parse_args(OptionParser(option_list = option_list))

names <- paste0(opt$tax_type, "_", opt$time, "_", opt$metab_type, "_", opt$model)
#TODO add names into files line 
files <- paste0("snakemake_output/analyses/prediction/",list.files("snakemake_output/analyses/prediction", pattern = "16S_12M_tar_svm"))
results <- lapply(files, readRDS)
correlation <- c()
for (i in seq_len(length(results))){
  correlation <- rbind(correlation, sapply(results[[i]], function(x){
    cor(x[,1], x[,2])
  }))
}
correlation <- t(correlation)

r2 <- c()
for (i in seq_len(length(results))){
  r2 <- rbind(r2, sapply(results[[i]], function(x){
    MLmetrics::R2_Score(x[,1], x[,2])
  }))
}
r2 <- t(r2)

met_names <- colnames(readRDS(file = paste0("data/processed/", opt$tax_type, "_", opt$time, "_", opt$metab_type, ".rds"))$met)
names_order <- as.vector(sapply(files, function(x){
  as.numeric(strsplit(strsplit(x, ".rds")[[1]], "_")[[1]][6])
}))

names(result) <- met_names[names_order]
colnames(r2) <- met_names[names_order]
colnames(correlation) <- met_names[names_order]
output <- list(raw = results, correlation = correlation, r2 = r2)

path <- paste0("snakemake_output/analyses/prediction_evaluation/", opt$tax_type, "_", opt$time, "_", opt$metab_type, "_", opt$model, ".rds")
saveRDS(output, file = path)



