# This script seeks to perform plotting of univariate and multivariate correlation analyses 
library(ggpubr)
library(ggplot2)
library(cowplot)
library(glue)
library(viridis)
library(ggplotify)
library(pheatmap)


correlation_plot <- function(data){
  # extract correlation matrix 
  cor_mat <- data$cor_mat
  if (ncol(cor_mat) > 36){
    show_col <- F
  } else {
    show_col <- T
  }
  sig <- t(apply(data$p_mat, 1, function(x){
    ifelse(x < 0.05, 1,0)
  }))
  
  # get tax names 
  tax_names <- paste(data$tax_tab[,c("Genus")]) #Family-Genus names
  tax_names <- as.expression(sapply(tax_names, function(x){
    bquote(italic(.(x)) ~ "spp.")
  }))
  # get met names
  met_names <- colnames(cor_mat)
  
  
  
  # set non-significant values to NA
  cor_mat[which(sig == 0, arr.ind = T)] <- NA
  
  # plot heatmap
  corr_heatmap <- pheatmap(
    mat = cor_mat,
    color = viridis(40),
    border_color = NA,
    show_colnames = show_col,
    show_rownames = T,
    labels_row = tax_names,
    labels_col = met_names,
    annotation_names_row = F,
    annotation_names_col = F,
    drop_levels = TRUE,
    fontsize = 14, 
    legend = F, cluster_rows = F, cluster_cols = F,
    na_col = "beige"
  )
  corr_heatmap <- as.ggplot(corr_heatmap)
  return(corr_heatmap)
}

scca_plot <- function(data, spearman_res){
  cor_mat <- spearman_res$cor_mat
  if (ncol(cor_mat) > 36){
    show_col <- F
  } else {
    show_col <- T
  }
  
  cca <- data$cca # main cca object
  sig <- t(apply(spearman_res$p_mat, 1, function(x){
    ifelse(x < 0.05, 1,0)
  }))
  
  tax_idx <- which(cca$u != 0)
  met_idx <- which(cca$v != 0)
  
  cca_mat <- cor_mat[tax_idx, met_idx]
  sig_mat <- sig[tax_idx, met_idx]
  
  tax_names <- paste(data$tab[tax_idx,][,c("Genus")]) #Family-Genus names
  tax_names <- as.expression(sapply(tax_names, function(x){
    bquote(italic(.(x)) ~ "spp.")
  }))
  
  # getting positive/negative annotation from scca object 
  row <- data.frame("sCCA Loading" = as.factor(ifelse(cca$u[tax_idx] > 0, "+","-")), check.names = F)
  rownames(row) <- rownames(cca_mat)
  col <- data.frame("sCCA Loading" = as.factor(ifelse(cca$v[met_idx] > 0, "+","-")), check.names = F)
  rownames(col) <- colnames(cca_mat)
  ann_colors <- list("sCCA Loading" = c(cividis(10)[1], cividis(10)[10]))
  names(ann_colors$"sCCA Loading") <- as.factor(c("+", "-")) 
  
  cca_mat[which(sig_mat == 0, arr.ind = 1)] <- NA
  
  reduced_heatmap <-  pheatmap(
    mat               = cca_mat,
    annotation_row    = row,
    annotation_colors = ann_colors,
    annotation_col    = col,
    color             = viridis(40),
    annotation_names_row = T,
    annotation_names_col = T,
    border_color      = NA,
    cluster_rows = F,
    cluster_cols = F,
    show_colnames     = show_col,
    show_rownames     = T,
    labels_row        = tax_names,
    drop_levels       = TRUE,
    fontsize          = 13,
    cex = 1, 
    legend = T,
    na_col = "beige"
  )
  reduced_heatmap <- as.ggplot(reduced_heatmap) 
  return(reduced_heatmap)
}


combine_plots_correlation <-  function(scca_plot, spearman_plot, scca_obj){
  raw_pval <- length(which(scca_obj$perm >= scca_obj$cca$cors))/length(scca_obj$perm)
  
  if (raw_pval == 0){
    raw_pval <- "< 0.001"
  } 
  title <- glue("SCCA Correlation: {correlation} \n (Bootstrapped 95% CI {upper} - {lower}; Permutational p-value: {pval})", correlation = round(scca_obj$cca$cors,3), 
                  upper = round(quantile(scca_obj$boot$cors, 0.05),3), 
                  lower = round(quantile(scca_obj$boot$cors, 0.95), 3),
                  pval = raw_pval)
  
  spearman_plot <- spearman_plot + theme(plot.title = element_blank())
  
  plt <- ggdraw() + draw_plot(spearman_plot, x = 0, y = 0, width = 0.6, height = 1) +
    draw_plot(scca_plot, x = 0.6, y = 0.2, width = 0.4, height = 0.6) + 
    draw_text(title, x = 0.8, y = 0.85)
  
  return(plt)
}




