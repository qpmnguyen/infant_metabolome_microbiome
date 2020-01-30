library(gridExtra)
library(ggplot2)
library(ggpubr)
library(ggplotify)
library(cowplot)
library(glue)



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
  draw_text("SCCA Correlation: 0.606 \n (Bootstrapped 95% CI: 0.61 - 0.73 ; Permutation p-value < 0.001)", 
            x = 0.8, y = 0.93, size = 11) +
  draw_text("SCCA Correlation: 0.52 \n (Bootstrapped 95% CI: 0.491 - 0.646 ; Permutation p-value < 0.001)", 
            x = 0.8, y = 0.43, size = 11) + 
  draw_plot_label(label = c("A. 6 Weeks", "B. 12 Months"), x = c(0,0), y = c(1,0.5))

save_plot(plot = tar_correlation_plot, "./docs/publication_figures/scc_plots_tar.png", dpi = 300, base_height = 23, 
          base_asp = 1.1)
