library(vegan)
library(phyloseq)
library(optparse)
library(Hmisc)

option_list <- list(
  make_option("--input", help = "Input file for data loading"),
  make_option("--output", help = "Output file for data loading"),
  make_option("--metric", help = "Correlation metric"),
  make_option("--MHC", help = "multiple hypothesis correction metric")
)

opt <- parse_args(OptionParser(option_list = option_list))

data <- readRDS(file = opt$input)

tax <- otu_table(tax)
met <- sample_data(met) 

#(met)[36] <- "pi-Methylhistidine" # unicode issues 

# simple correlation
correlation <- cor(tax, met, method = opt$metric)

# pairwise correlation
p_mat <- t(apply(tax, 2, function(x){
  p_vals <- c()
  for (i in 1:ncol(met)){
    p_vals[i] <- cor.test(x = x, y = met[,i], method = opt$metric)$p.value
  }
  return(p_vals)
}))

# BH adjustment 
adj_mat <- matrix(p.adjust(as.vector(p_mat), method = "BH"), ncol = ncol(p_mat), nrow = nrow(p_mat), byrow = F) # MHC adjustment

result <- list(
  cor_mat = correlation,
  p_mat = adj_mat,
  tax_tab = tax_table(data)
)

saveRDS(file = opt$output, object = result)