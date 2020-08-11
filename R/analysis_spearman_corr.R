# This script contains function to do pairwise spearman correlation with BH correction
library(vegan)
library(phyloseq)
library(Hmisc)

spearman_main <- function(data){
  tax <- as(otu_table(data))
  met <- as(sample_data(data))
  correlation <- cor(tax, met, method = "spearman")
  # pairwise testing correlation
  p_mat <- t(apply(tax, 2, function(x){
    p_vals <- c()
    print(length(x))
    for (i in 1:ncol(met)){
      p_vals[i] <- cor.test(x = x, y = met[,i], method = "spearman")$p.value
    }
    return(p_vals)
  }))
  # p-value adjustment using BH
  adj_mat <- matrix(p.adjust(as.vector(p_mat), method = "BH"), 
         ncol = ncol(p_mat), nrow = nrow(p_mat), byrow = F)
  output <- list(
    cor_mat = correlation,
    p_mat = adj_mat,
    tax_tab = as(tax_table(data),"matrix")
  )
  return(output)
}