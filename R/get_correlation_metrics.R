library(phyloseq)
library(glue)
library(optparse)

untar <- readRDS(file = "output/analyses/correlation/6W_untar_spearman.rds")
tar <- readRDS(file = "output/analyses/correlation/6W_tar_spearman.rds")

# Getting number of columns or rows with significant values
get_number_significant <- function(spearman_obj){
  cor_mat <- spearman_obj$cor_mat
  cor_mat[which(spearman_obj$p_mat >= 0.05)] <- NA
  arr_idx <- which(!is.na(cor_mat), arr.ind = T)
  print(glue("Number of rows {n_row} with at least one significant value",
             n_row = length(unique(arr_idx[,1]))))
  print(glue("Number of columns {n_row} with at least one significant value",
             n_row = length(unique(arr_idx[,2]))))
  
}

# print the number of significant correlations as a portion of overall cells 
get_sig_corr <- function(spearman_obj){
  n_cells <- length(as.vector(spearman_obj$cor_mat))
  sig_cells <- length(which(spearman_obj$p_mat < 0.05))
  print(glue("Proportion of cells that are significant {prop}", prop = round(100*(sig_cells/n_cells),2)))
}

get_number_significant(untar)
get_number_significant(tar)

get_sig_corr(untar)
get_sig_corr(tar)
