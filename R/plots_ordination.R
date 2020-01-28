library(ggpubr)
library(cowplot)
library(pheatmap)
library(viridis)
library(phyloseq)
library(reshape2)
library(ggplot2)

input <- readRDS(file = "output/analyses/ordinations/6W_tar_ordinations.rds")

input <- input$tax_dist_Gunifrac_met_dist_euclidean

tax_ord <- input$tax_ord
met_ord <- input$met_ord
proc_test <- input$proc_test


plotting_ord <- function(ord_obj){
  pts <- as.data.frame(ord_obj$points[,c(1,2)])
  perc <- c(abs(ord_obj$eig[1])/sum(abs(ord_obj$eig)) * 100, abs(ord_obj$eig[2])/sum(abs(ord_obj$eig)) * 100)
  colnames(pts) <- c("MDS1", "MDS2")
  plt <- ggplot(pts, aes(x = MDS1, y = MDS2)) + geom_point(color = viridis(100)[50]) + labs(x = "NMDS1", y = "NMDS2") +
    theme_pubr()
  return(plt)
}
plotting_ord(tax_ord)
