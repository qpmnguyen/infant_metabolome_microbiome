library(ggpubr)
library(cowplot)
library(pheatmap)
library(viridis)
library(phyloseq)
library(reshape2)
library(ggplot2)
library(glue)
library(optparse)

option_list <- list(
  make_option("--input", help = "Input file for data loading"),
  make_option("--tax_dist", help = "dist type for taxonomic data"),
  make_option("--met_dist", help = "dist type for metabnomic data")
)

opt <- parse_args(OptionParser(option_list = option_list))

input <- readRDS(file = opt$input)
time <- tail(strsplit(strsplit(opt$input, "_")[[1]][1],"/")[[1]],1)
metab <- strsplit(opt$input,"_")[[1]][2]

# defining some plotting functions  
plotting_ord <- function(ord_obj){
  pts <- as.data.frame(ord_obj$vectors[,c(1,2)])
  if(is.null(ord_obj$trace.cor) == T){
    perc <- round(as.vector(ord_obj$values$Eigenvalues/ord_obj$trace)[c(1,2)] * 100,2)
  } else {
    perc <- round(as.vector(ord_obj$values$Corr_eig/ord_obj$trace.cor)[c(1,2)] * 100,2)
  }
  print(perc)
  colnames(pts) <- c("MDS1", "MDS2")
  plt <- ggplot(pts, aes(x = MDS1, y = MDS2)) + geom_point(color = viridis(100)[50], size = 2) + 
    labs(x = glue("MDS1 ({perc}%)", perc = perc[1]), y = glue("MDS2 ({perc}%)", perc = perc[2])) +
    theme_pubr()
  return(plt)
}

plotting_proc <- function(proc_test){
  pts <- data.frame(x1 = proc_test$X[,1], x2 = proc_test$X[,2], y1 = proc_test$Yrot[,1], y2 = proc_test$Yrot[,2])
  plt <- ggplot(pts) + geom_point(aes(x = x1, y = x2, color = "Taxonomy"), size = 2.5) +  
    geom_point(aes(x = y1, y = y2, color = "Metabolites"), size = 2.5) + 
    stat_ellipse(mapping = aes(x = x1, y = x2, fill = "Taxonomy"), geom = "polygon", show.legend = F, alpha = 0.2) + 
    stat_ellipse(mapping = aes(x = y1, y = y2, fill = "Metabolites"), geom = "polygon", show.legend = F, alpha = 0.2) + 
    scale_color_manual(name = "Ordination", values = viridis(100)[c(1,50)]) + 
    scale_fill_manual(name = "Ordination", values = viridis(100)[c(1,50)]) + 
    labs(x = "MDS1", y = "MDS2", title = glue("Procrustes analysis of samples ordinated by taxonomy and metabolites \n (p={sig}, ss={ss})", 
                                             sig = proc_test$signif, ss = round(proc_test$ss,2))) + 
    theme_pubr() + theme(plot.title=element_text(hjust = 0.5)) 
    # geom_segment(aes(x = y1, y = y2, xend = x1, yend = x2), arrow = arrow(length = unit(0.2, "cm")), alpha = 0.3)
  return(plt)
}

input <- input[[glue("tax_{tax_dist}_met_{met_dist}", tax_dist = opt$tax_dist, met_dist = opt$met_dist)]]

### perform plotting 
tax_ord <- input$tax_ord
met_ord <- input$met_ord
proc_test <- input$proc_test
plot_tax <- plotting_ord(tax_ord)
plot_met <- plotting_ord(met_ord)
plot_proc <- plotting_proc(proc_test)
for (j in c("svg","png")){
  ggsave(plot = plot_tax, filename = glue("output/figures/ordinations/{tax_dist}_{met_dist}/{time}_{metab}_tax.{ext}", 
                                          time = time, metab = metab, tax_dist = opt$tax_dist, met_dist = opt$met_dist, ext = j), device = j)
  ggsave(plot = plot_met, filename = glue("output/figures/ordinations/{tax_dist}_{met_dist}/{time}_{metab}_met.{ext}", 
                                          time = time, metab = metab, met_dist = opt$met_dist, tax_dist = opt$tax_dist, ext = j), device = j)
  ggsave(plot = plot_proc, filename = glue("output/figures/ordinations/{tax_dist}_{met_dist}/{time}_{metab}_joint.{ext}", 
                                          time = time, metab = metab, tax_dist = opt$tax_dist, met_dist = opt$met_dist, ext = j), device = j)
}



