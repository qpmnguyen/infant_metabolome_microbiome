library(vegan)
library(phyloseq)
library(optparse)
library(Hmisc)

option_list <- list(
  make_option("--input", help = "Input file for data loading"),
  make_option("--output", help = "Output file for data loading")
)

opt <- parse_args(OptionParser(option_list = option_list))

data <- readRDS(file = opt$input)

tax <- data$tax
met <- data$met 

correlation <- cor(tax, met, method = "spearman")
p_mat <- t(apply(tax, 2, function(x){
  p_vals <- c()
  for (i in 1:ncol(met)){
    p_vals[i] <- cor.test(x = x, y = met[,i], method = "spearman")$p.value
  }
  return(p_vals)
}))


adj <- matrix(p.adjust(as.vector(p_mat), method = "bonferroni"), ncol = ncol(p_mat), nrow = nrow(p_mat), byrow = F)

sig <- t(apply(adj, 1, function(x){
  ifelse(x < 0.05, 1,0)
}))

idx_tax <- which(apply(sig,1, function(x){
  all(x == 0)
}) == F)
idx_met <- which(apply(sig,2, function(x){
  all(x == 0)
}) == F)

sig <- ifelse(sig == 1, "*", "")
pheatmap(correlation, color = viridis(20), display_numbers = sig)

# Re-adjust correlation calculation and then convert all plotting functionality into a new script. 
