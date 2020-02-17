library(gridExtra)
library(ggplot2)
library(ggpubr)
library(ggplotify)
library(cowplot)
library(glue)


# Plotting figure for correlation plots ####
p1_6W <- readRDS(file = "./output/figures/correlation/6W_tar_correlation.rds")
p2_6W <- readRDS(file = "./output/figures/correlation/6W_tar_scca.rds")
p1_12M <- readRDS(file = "./output/figures/correlation/12M_tar_correlation.rds")
p2_12M <- readRDS(file = "./output/figures/correlation/12M_tar_scca.rds")

p2_6W <- p2_6W + theme(plot.title = element_blank())
p2_12M <- p2_12M + theme(plot.title = element_blank())

tar_correlation_plot <- ggdraw() + 
  draw_plot(p1_6W, x = 0, y = 0.5, width = 0.6, height = 0.45) + 
  draw_plot(p2_6W , x = 0.6, y = 0.6, width = 0.4, height = 0.3) + 
  draw_plot(p1_12M, x = 0, y = 0, width = 0.6, height = 0.45) + 
  draw_plot(p2_12M, x = 0.6, y = 0.1,  width = 0.4, height = 0.3) + 
  draw_text("sCCA Correlation (1st variate pair): 0.606 \n (Bootstrapped 95% CI: 0.61 - 0.73 ; Permutation p-value < 0.001)", 
            x = 0.8, y = 0.93, size = 20) +
  draw_text("sCCA Correlation (1st variate pair): 0.52 \n (Bootstrapped 95% CI: 0.491 - 0.646 ; Permutation p-value < 0.001)", 
            x = 0.8, y = 0.43, size = 20) + 
  draw_plot_label(label = c("A. 6 Weeks", "B. 12 Months"), x = c(0,0), y = c(1,0.5), size = 30)

save_plot(plot = tar_correlation_plot, "./docs/publication_figures/scc_plots_tar.png", dpi = 300, base_height = 23, 
          base_asp = 1.1)

# plotting supplemental figure similar to figure 2 but for the untar plots #### 
p1_6W <- readRDS(file = "./output/figures/correlation/6W_untar_correlation.rds")
p2_6W <- readRDS(file = "./output/figures/correlation/6W_untar_scca.rds")
p1_12M <- readRDS(file = "./output/figures/correlation/12M_untar_correlation.rds")
p2_12M <- readRDS(file = "./output/figures/correlation/12M_untar_scca.rds")

p2_6W <- p2_6W + theme(plot.title = element_blank())
p2_12M <- p2_12M + theme(plot.title = element_blank())

(untar_correlation_plot <- ggdraw() + 
  draw_plot(p1_6W, x = 0, y = 0.5, width = 0.6, height = 0.45) + 
  draw_plot(p2_6W , x = 0.6, y = 0.6, width = 0.4, height = 0.3) + 
  draw_plot(p1_12M, x = 0, y = 0, width = 0.6, height = 0.45) + 
  draw_plot(p2_12M, x = 0.6, y = 0.1,  width = 0.4, height = 0.3) + 
  draw_text("sCCA Correlation (1st variate pair): 0.636 \n (Bootstrapped 95% CI: 0.621 - 0.733 ; Permutation p-value < 0.001)", 
            x = 0.8, y = 0.93, size = 20) +
  draw_text("sCCA Correlation (1st variate pair): 0.49 \n (Bootstrapped 95% CI: 0.475 - 0.702 ; Permutation p-value < 0.001)", 
            x = 0.8, y = 0.43, size = 20) + 
  draw_plot_label(label = c("A. 6 Weeks", "B. 12 Months"), x = c(0,0), y = c(1,0.5), size = 30))

save_plot(plot = untar_correlation_plot, "./docs/publication_figures/scc_plots_untar.png", dpi = 300, base_height = 23, base_width = 25, 
          base_asp = 1.1)

# plotting ordination plot tar ####
p1_6W <- readRDS(file = "output/figures/ordinations/gunifrac_euclidean/6W_tar_joint.rds")
p1_12M <- readRDS(file = "output/figures/ordinations/gunifrac_euclidean/12M_tar_joint.rds")
p2_6W <- readRDS(file = "output/figures/ordinations/euclidean_euclidean/6W_tar_joint.rds")
p2_12M <- readRDS(file = "output/figures/ordinations/euclidean_euclidean/12M_tar_joint.rds")
legend <- get_legend(
  # create some space to the left of the legend
  p1_6W + theme(legend.box.margin = margin(0, 0, 0, 12))
)

p1_6W <- p1_6W + theme(plot.title = element_blank(), legend.position = "none")
p2_6W <- p2_6W + theme(plot.title = element_blank(), legend.position = "none")
p1_12M <- p1_12M + theme(plot.title = element_blank(), legend.position = "none")
p2_12M <- p2_12M + theme(plot.title = element_blank(), legend.position = "none")

ordination_tar <- ggdraw() + 
  draw_text("p-value = 0.057 \n Sum of Squares: 0.98", x = 0.3, y = 0.97, size = 12) + 
  draw_text("A. Euclidean-Euclidean", fontface = "bold", x = 0, hjust = 0, y = 0.97) + 
  draw_plot(p2_6W, x = 0, y = 0.53, height = 0.38, width = 0.5) + 
  draw_text("p-value = 0.001 \n Sum of Squares: 0.95", x = 0.8, y = 0.97, size = 12) +
  draw_plot(p2_12M, x = 0.5, y = 0.53, height = 0.38, width = 0.5) +
  draw_text("B. Gunifrac-Euclidean", fontface = "bold", x = 0, hjust = 0, y = 0.5) + 
  draw_text("p-value = 0.376 \n Sum of Squares: 0.99", x = 0.3, y = 0.5, size = 12, hjust = 0.5) +
  draw_plot(p1_6W, x = 0, y = 0.07, height = 0.38, width = 0.5) + 
  draw_text("p-value = 0.069 \n Sum of Squares: 0.99", x = 0.8, y = 0.5, size = 12) +
  draw_plot(p1_12M, x = 0.5, y = 0.07 , height = 0.38, width = 0.5) + 
  draw_text("6 weeks", fontface = "bold", x = 0.3, y = 0.05) + 
  draw_text("12 months", fontface = "bold", x = 0.8, y = 0.05) + 
  draw_plot(legend, x = 0, y = 0.02, vjust = 0.5, height = 0.02, width = 1)

  
save_plot(plot = ordination_tar, "./docs/publication_figures/ordination_plots_tar.png", dpi = 300, 
          base_height = 10, base_width = 12)

# plotting ordination plot untar ####
p1_6W <- readRDS(file = "output/figures/ordinations/gunifrac_euclidean/6W_untar_joint.rds")
p1_12M <- readRDS(file = "output/figures/ordinations/gunifrac_euclidean/12M_untar_joint.rds")
p2_6W <- readRDS(file = "output/figures/ordinations/euclidean_euclidean/6W_untar_joint.rds")
p2_12M <- readRDS(file = "output/figures/ordinations/euclidean_euclidean/12M_untar_joint.rds")
legend <- get_legend(
  # create some space to the left of the legend
  p1_6W + theme(legend.box.margin = margin(0, 0, 0, 12))
)

p1_6W <- p1_6W + theme(plot.title = element_blank(), legend.position = "none")
p2_6W <- p2_6W + theme(plot.title = element_blank(), legend.position = "none")
p1_12M <- p1_12M + theme(plot.title = element_blank(), legend.position = "none")
p2_12M <- p2_12M + theme(plot.title = element_blank(), legend.position = "none")

ordination_untar <- ggdraw() + 
  draw_text("p-value = 0.043 \n Sum of Squares: 0.98", x = 0.3, y = 0.97, size = 12) + 
  draw_text("A. Euclidean-Euclidean", fontface = "bold", x = 0, hjust = 0, y = 0.97) + 
  draw_plot(p2_6W, x = 0, y = 0.53, height = 0.38, width = 0.5) + 
  draw_text("p-value = 0.001 \n Sum of Squares: 0.94", x = 0.8, y = 0.97, size = 12) +
  draw_plot(p2_12M, x = 0.5, y = 0.53, height = 0.38, width = 0.5) +
  draw_text("B. Gunifrac-Euclidean", fontface = "bold", x = 0, hjust = 0, y = 0.5) + 
  draw_text("p-value = 0.001 \n Sum of Squares: 0.94", x = 0.3, y = 0.5, size = 12, hjust = 0.5) +
  draw_plot(p1_6W, x = 0, y = 0.07, height = 0.38, width = 0.5) + 
  draw_text("p-value = 0.006 \n Sum of Squares: 0.94", x = 0.8, y = 0.5, size = 12) +
  draw_plot(p1_12M, x = 0.5, y = 0.07 , height = 0.38, width = 0.5) + 
  draw_text("6 weeks", fontface = "bold", x = 0.3, y = 0.05) + 
  draw_text("12 months", fontface = "bold", x = 0.8, y = 0.05) + 
  draw_plot(legend, x = 0, y = 0.02, vjust = 0.5, height = 0.02, width = 1)

save_plot(plot = ordination_untar, "./docs/publication_figures/ordination_plots_untar.png", dpi = 300, 
          base_height = 10, base_width = 12)

# plotting evaluation of prediction with borda ####
b1 <- readRDS(file = "output/figures/prediction/r2_tar_violinplots.rds")
b2 <- readRDS(file = "output/figures/prediction/corr_tar_violinplots.rds")
bor2 <- readRDS(file = "output/figures/prediction/r2_tar_bordaplots.rds")
bor1 <- readRDS(file = "output/figures/prediction/corr_tar_bordaplots.rds")

model_comparison <- ggdraw() + draw_plot(b1 + theme(plot.title = element_text(hjust = 0.5)), x = 0, y = 0.51, height  = 0.47, width = 0.5) +
  draw_plot(b2 + theme(plot.title = element_text(hjust = 0.5)), x = 0.5, y = 0.51, height = 0.47, width = 0.5) + 
  draw_plot(bor1 + theme(legend.position = "none"), x = 0, y = 0.03, height = 0.47, width = 0.5) + 
  draw_plot(bor2 + theme(legend.position = "none"), x = 0.5, y = 0.03, height = 0.47, width = 0.5) 

save_plot(model_comparison, filename = "docs/publication_figures/model_comparison.png", dpi = 300, base_width = 6, base_height = 6)


# plotting results heatmap 
heatmap1 <- readRDS(file = "output/figures/prediction/r2_tar_heatmap.rds")
heatmap2 <- readRDS(file = "output/figures/prediction/corr_tar_heatmap.rds")
heatmap3 <- readRDS(file = "output/figures/prediction/rmse_tar_heatmap.rds")
heatmap_cobble <- ggdraw() + draw_plot(heatmap1, x = 0, y = 0.66, height = 0.30, width = 1) +
  draw_plot(heatmap2, x = 0, y = 0.33, height = 0.30, width = 1) + 
  draw_plot(heatmap3, x = 0, y = 0.02, height = 0.30, width = 1) + 
  draw_label("A. R-squared", fontface = "bold", x = 0, y = 0.98, hjust = 0, size = 15) +
  draw_label("B. Correlation", fontface = "bold", x = 0, y = 0.64, hjust = 0, size = 15) + 
  draw_label("C. RMSE", fontface = "bold", x = 0, y = 0.32, hjust = 0, size = 15) + 
  draw_text("6 Weeks", x = 0.2, y = 0, vjust = 0, hjust = 0, fontface = "bold", size = 15) + 
  draw_text("12 Months", x = 0.7, y = 0, vjust = 0, hjust = 0, fontface = "bold", size = 15)
save_plot(heatmap_cobble, file = "docs/publication_figures/heatmap_cobble_tar.png", dpi = 300, base_height = 15.5, base_width = 12)

# plotting rankings comparison
p1 <- readRDS(file = "output/figures/prediction/ranking_ordinations_tar_rmse_6W.rds")
p2 <- readRDS(file = "output/figures/prediction/ranking_ordinations_tar_rmse_12M.rds")
p3 <- readRDS(file = "output/figures/prediction/ranking_ordinations_tar_corr_6W.rds")
p4 <- readRDS(file = "output/figures/prediction/ranking_ordinations_tar_corr_12M.rds")

rankings_pcoa <- ggdraw() + draw_plot(p1, x = 0, y = 0.51, height  = 0.47, width = 0.5) +
  draw_plot(p2, x = 0.5, y = 0.51, height = 0.47, width = 0.5) + 
  draw_plot(p3, x = 0, y = 0.03, height = 0.47, width = 0.5) + 
  draw_plot(p4, x = 0.5, y = 0.03, height = 0.47, width = 0.5) + 
  draw_label("A. RMSE", fontface = "bold", x = 0, y = 0.99, hjust = 0, size = 15) + 
  draw_label("B. Correlation", fontface = "bold", x = 0, y = 0.51, hjust = 0, size = 15) + 
  draw_text("6 Weeks", x = 0.2, y = 0.01, vjust = 0, hjust = 0, fontface = "bold", size = 15) + 
  draw_text("12 Months", x = 0.7, y = 0.01, vjust = 0, hjust = 0, fontface = "bold", size = 15)
save_plot(rankings_pcoa, file = "docs/publication_figures/rankings_cobble_tar.png", dpi = 300, base_height = 11, base_width = 12)

# plotting forest plots 
p1 <- readRDS(file = "output/figures/prediction/6W_r2_forestplot.rds")
p2 <- readRDS(file = "output/figures/prediction/12M_r2_forestplot.rds")
p3 <- readRDS(file = "output/figures/prediction/6W_corr_forestplot.rds")
p4 <- readRDS(file = "output/figures/prediction/12M_corr_forestplot.rds")

forest_plot <- ggdraw() + draw_plot(p1 + theme(axis.title.x = element_blank(), panel.spacing = unit(1, "lines")), x = 0, y = 0.51, height  = 0.47, width = 0.5) +
  draw_plot(p2 + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
                       panel.spacing = unit(1, "lines")), x = 0.5, y = 0.51, height = 0.47, width = 0.5) + 
  draw_plot(p3 + theme(panel.spacing = unit(1, "lines")), x = 0, y = 0.03, height = 0.47, width = 0.5) + 
  draw_plot(p4 + theme(axis.title.y = element_blank(), panel.spacing = unit(1, "lines")), x = 0.5, y = 0.03, height = 0.47, width = 0.5) + 
  draw_label("A. R2", fontface = "bold", x = 0, y = 0.99, hjust = 0, size = 15) + 
  draw_label("B. Correlation", fontface = "bold", x = 0, y = 0.51, hjust = 0, size = 15) + 
  draw_text("6 Weeks", x = 0.2, y = 0.01, vjust = 0, hjust = 0, fontface = "bold", size = 15) + 
  draw_text("12 Months", x = 0.7, y = 0.01, vjust = 0, hjust = 0, fontface = "bold", size = 15)
forest_plot
save_plot(forest_plot, file = "docs/publication_figures/forresplot_r2_corr_tar.png", dpi = 300, base_height = 15, base_width = 23)



