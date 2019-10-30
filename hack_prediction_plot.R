library(ggplot2)
library(ggpubr)
library(viridis)
library(reshape2)
library(dplyr)
library(cowplot)

one_yr <- readRDS(file = "snakemake_output/analyses/prediction_evaluation/16S_12M_tar_rf.rds")
six_wks <- readRDS(file = "snakemake_output/analyses/prediction_evaluation/16S_6W_tar_rf.rds")

cor_12M <- melt(one_yr$correlation)
cor_12M$Var1 <- rep("12M", nrow(cor_12M))
cor_6W <- melt(six_wks$correlation)
cor_6W$Var1 <- rep("6W", nrow(cor_6W))

combined_corr <- rbind(cor_12M, cor_6W)

cor_plt <- ggboxplot(combined, x = "Var2", y = "value", fill = "Var2", orientation = "horizontal",
                 add = "jitter")
cor_plt <- ggpar(cor_plt, palette = viridis(36), legend = "none", xlab = "Metabolite", ylab = "Pearson Correlation")
cor_plt <- facet(cor_plt, facet.by = "Var1")

### R-squared -------------------------------------------------------------------------------
r2_12M <- melt(one_yr$r2)
r2_12M$Var1 <- rep("12M", nrow(r2_12M))
r2_6W <- melt(six_wks$r2)
r2_6W$Var1 <- rep("6W", nrow(r2_6W))

combined_r2 <- rbind(r2_12M, r2_6W)

r2_plt <- ggboxplot(combined, x = "Var2", y = "value", fill = "Var2", orientation = "horizontal",
                     add = "jitter")
r2_plt <- ggpar(r2_plt, palette = viridis(36), legend = "none", xlab = "Metabolite", ylab = "Pseudo R-squared")
r2_plt <- facet(r2_plt, facet.by = "Var1")


### best and worse correlation -----------------------------------------------
trimmed_mean <- combined_corr %>% group_by(Var1, Var2) %>% summarise(mean = mean(value, trim = 0.1))
max <- trimmed_mean[trimmed_mean$mean == max(trimmed_mean$mean),]
min <- trimmed_mean[trimmed_mean$mean == min(trimmed_mean$mean),]
max_frame <- six_wks$raw[[max$Var2]]
min_frame <- one_yr$raw[[min$Var2]]

max_frame <- lapply(1:length(max_frame), function(i){
  max_frame[[i]] <- as.data.frame(max_frame[[i]])
  max_frame[[i]]$Fold <- rep(paste0("Fold ",i), nrow(max_frame[[i]]))
  return(max_frame[[i]])
})

min_frame <- lapply(1:length(min_frame), function(i){
  min_frame[[i]] <- as.data.frame(min_frame[[i]])
  min_frame[[i]]$Fold <- rep(paste0("Fold ",i), nrow(min_frame[[i]]))
  return(min_frame[[i]])
})

max_frame <- do.call(rbind, max_frame)
min_frame <- do.call(rbind, min_frame)

max_plot <- ggscatter(max_frame, x = "met.test", y = "predictions", color = "Fold", palette = viridis(5), add = "reg.line", conf.int = T)
max_plot <- ggpar(max_plot, xlab = "True values", ylab = "Predictions", main = "Butyrate at 6 weeks - Trimmed mean correlation = 0.519")

min_plot <- ggscatter(min_frame, x = "met.test", y = "predictions", color = "Fold", palette = viridis(5), add = "reg.line", conf.int = T)
min_plot <- ggpar(min_plot, xlab = "True values", ylab = "Predictions", main = "Glycerol at 12 months - Trimmed mean correlation = -0.0252")

leg <- get_legend(min_plot)


grid <- plot_grid(r2_plt, 
          cor_plt, 
          max_plot + theme(legend.position = "none"), 
          min_plot + theme(legend.position = "none"), 
          nrow = 2, labels = c("A", "B", "C", "D"))

grid <- plot_grid(grid, leg, ncol = 1, rel_heights = c(1,.1))

ggsave(grid, file = "snakemake_output/figures/predictions/prediction_result_figures.png", 
       device = "png", width = 15, height = 15)

### Model comparison boxplots ------------------------------------------------------------------
models <- c("svm", "enet", "rf", "")





