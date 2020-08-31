library(ggpubr)
library(cowplot)
library(viridis)
library(phyloseq)
library(reshape2)
library(ggplot2)
library(glue)


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
    labs(x = "MDS1", y = "MDS2", title = glue("p-value = {sig} \n Sum of Squares: {ss}", 
                                             sig = proc_test$signif, ss = round(proc_test$ss,2))) + 
    theme_pubr() + theme(plot.title=element_text(hjust = 0.5)) 
    # geom_segment(aes(x = y1, y = y2, xend = x1, yend = x2), arrow = arrow(length = unit(0.2, "cm")), alpha = 0.3)
  return(plt)
}

combine_plot_ord <- function(plot1, plot2, plot3, plot4){
  legend <- get_legend(
    # create some space to the left of the legend
    plot1 + theme(legend.box.margin = margin(0, 0, 0, 12))
  )
  
  plot1 <- plot1 + theme(plot.title = element_blank(), legend.position = "none")
  plot2 <- plot2 + theme(plot.title = element_blank(), legend.position = "none")
  plot3 <- plot3 + theme(plot.title = element_blank(), legend.position = "none")
  plot4 <- plot4 + theme(plot.title = element_blank(), legend.position = "none")
  
  ordination_tar <- ggdraw() + 
    draw_text(plot3$labels$title, x = 0.3, y = 0.97, size = 12) + 
    draw_text("A. Euclidean-Euclidean", fontface = "bold", x = 0, hjust = 0, y = 0.97) + 
    draw_plot(plot3, x = 0, y = 0.53, height = 0.38, width = 0.5) + 
    draw_text(plot4$labels$title, x = 0.8, y = 0.97, size = 12) +
    draw_plot(plot4, x = 0.5, y = 0.53, height = 0.38, width = 0.5) +
    draw_text("B. Gunifrac-Euclidean", fontface = "bold", x = 0, hjust = 0, y = 0.5) + 
    draw_text(plot1$labels$title, x = 0.3, y = 0.5, size = 12, hjust = 0.5) +
    draw_plot(plot1, x = 0, y = 0.07, height = 0.38, width = 0.5) + 
    draw_text(plot2$labels$title, x = 0.8, y = 0.5, size = 12) +
    draw_plot(plot2, x = 0.5, y = 0.07 , height = 0.38, width = 0.5) + 
    draw_text("6 weeks", fontface = "bold", x = 0.3, y = 0.05) + 
    draw_text("12 months", fontface = "bold", x = 0.8, y = 0.05) + 
    draw_plot(legend, x = 0, y = 0.02, vjust = 0.5, height = 0.02, width = 1)

  return(ordination_tar)
}




