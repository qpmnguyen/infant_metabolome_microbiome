library(optparse)
library(ggplot2)
library(viridis)
library(ggpubr)
library(MLmetrics)
library(cowplot)
library(reshape2)
library(dplyr)

options_list <- list(
  make_option("--model", help = "Selected model for evaluation"),
  make_option("--metab_type", help = "Type of metabolomics data"),
  make_option("--time", help = "Data time point to extract"),
  make_option("--tax_type", help = "Data type for taxonomic data, so far only '16S' is supported")
)

names <- paste0("snakemake_output/analyses/prediction_evaluation/", opt$tax_type, "_", opt$time, "_", opt$metab_type, "_", opt$model, ".rds")
data <- readRDS(file = names)

raw <- data$raw
correlation <- data$correlation
r2 <- data$r2

correlation <- melt(correlation)
corr_plt <- ggboxplot(correlation, x = "Var2", y= "value", orientation = "horizontal", notch = F, add = "jitter", fill = "Var2", palette = viridis(36))
corr_plt <- ggpar(corr_plt, xlab = "Metabolites", ylab = "Pearson Correlation", legend = "none")

r2 <- melt(r2)
r2_plt <- ggboxplot(r2, x = 'Var2', y = "value", orientation = "horizontal", add = "jitter", fill  = "Var2", palette = viridis(36))
r2_plt <- ggpar(r2_plt, xlab = "Metabolites", ylab = "Pseudo R-squared", legend = "none")

mean_correlation <- correlation %>% group_by(Var2) %>% summarise(mean = mean(value, trim = 0.1)) 
max_correlation <- as.character(mean_correlation[mean_correlation$mean == max(mean_correlation$mean),]$Var2)
min_correlation <- as.character()


