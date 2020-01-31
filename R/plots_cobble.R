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
p1_6W <- readRDS(file = "output/figures/ordinations/gunifrac_manhattan/6W_tar_joint.rds")
p2_12M <- readRDS(file = "output/figures/ordinations/gunifrac_manhattan/12M_tar_joint.rds")
legend <- get_legend(
  # create some space to the left of the legend
  p1_6W + theme(legend.box.margin = margin(0, 0, 0, 12))
)

ordination_tar <- ggdraw() + 
  draw_text("p-value = 0.376 \n Sum of Squares: 0.99", x = 0.3, y = 0.85, size = 15) + 
  draw_plot(p1_6W + theme(plot.title = element_blank(), legend.position = "none"), x = 0, y = 0.05, width = 0.5, height = 0.75) + 
  draw_text("p-value = 0.069 \n Sum of Squares: 0.99", x = 0.8, y = 0.85, size = 15) +
  draw_plot(p2_12M + theme(plot.title = element_blank(), legend.position = "none"), x = 0.5, y = 0.05, width = 0.5, height = 0.75) + 
  draw_plot(legend, x = 0, y = 0, height = 0.1) + 
  draw_plot_label(label = c("A. 6 Weeks", "B. 12 Months"), x = c(0,0.5), y = c(0.95, 0.95), size = 25)

save_plot(plot = ordination_tar, "./docs/publication_figures/ordination_plots_tar.png", dpi = 300, base_height = 10, base_width = 12)

# plotting ordination plot untar ####
p1_6W <- readRDS(file = "output/figures/ordinations/gunifrac_manhattan/6W_untar_joint.rds")
p2_12M <- readRDS(file = "output/figures/ordinations/gunifrac_manhattan/12M_untar_joint.rds")
legend <- get_legend(
  # create some space to the left of the legend
  p1_6W + theme(legend.box.margin = margin(0, 0, 0, 12))
)

ordination_untar <- ggdraw() + 
  draw_text("p-value = 0.001 \n Sum of Squares: 0.94", x = 0.3, y = 0.85, size = 15) + 
  draw_plot(p1_6W + theme(plot.title = element_blank(), legend.position = "none"), x = 0, y = 0.05, width = 0.5, height = 0.75) + 
  draw_text("p-value = 0.006 \n Sum of Squares: 0.94", x = 0.8, y = 0.85, size = 15) +
  draw_plot(p2_12M + theme(plot.title = element_blank(), legend.position = "none"), x = 0.5, y = 0.05, width = 0.5, height = 0.75) + 
  draw_plot(legend, x = 0, y = 0, height = 0.1) + 
  draw_plot_label(label = c("A. 6 Weeks", "B. 12 Months"), x = c(0,0.5), y = c(0.95, 0.95), size = 25)

save_plot(plot = ordination_untar, "./docs/publication_figures/ordination_plots_untar.png", dpi = 300, base_height = 10, base_width = 12)
