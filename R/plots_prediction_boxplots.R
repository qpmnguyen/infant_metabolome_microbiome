library(ggplot2)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(ggplotify)
library(glue)
library(stringr)
library(viridis)
library(gridExtra)
library(pheatmap)
library(factoextra)
library(Rankcluster)

# plotting levels violin plots take advantage of the entire distribution ####
plot_comparison_per_metabolite <- function(eval, time, met){
  results <- readRDS(file = glue("output/analyses/prediction/processed/{met}_{eval}_by_met.rds", eval = eval, met = met))
  grid_plt <- list()
  for (i in 1:length(results)){
    data <- results[[i]] %>% select(starts_with(time), id, id2) %>% group_by(id) %>% 
                                       summarise_each(funs(mean), -id2) %>% 
                                       pivot_longer(-id, names_to = "model", values_to = "RMSE") %>% 
                                       mutate(model = str_remove(model, paste0(time,"_")))
    plt <- ggplot(data,aes(x = model, y = RMSE)) + geom_violin(aes(fill = model)) + geom_boxplot(alpha = 0.3) + coord_flip() + theme_pubr() + scale_fill_viridis_d() 
    if (i != 1){
      plt <- plt + theme(legend.position = "none", axis.title = element_blank(), axis.ticks.y = element_blank())
    } else {
      legend <- get_legend(plt)
      plt <- plt + theme(legend.position = "none", axis.title = element_blank(), axis.ticks.y = element_blank())
    }
    grid_plt[[i]] <- plt
  }
  g1 <- plot_grid(plotlist = grid_plt, labels = names(results), nrow = 6, ncol = 6, label_size = 12)
  perf <- ggdraw() + draw_plot(g1, x = 0.02, y = 0.02, height = 0.92, width = 0.98) +
    draw_text("Model", x = 0.01, y = 0.5, angle = 90, size = 15) + draw_text("RMSE", x = 0.5, y = 0.01, size = 15) + 
    draw_plot(legend, x = 0, y = 0.95, width = 1, height = 0.05)
  saveRDS(perf, file = glue("output/figures/prediction/{time}_{met}_{eval}_model_eval_boxplots.rds", time = time, eval = eval, met = met))
  save_plot(perf, filename = glue("output/figures/prediction/{time}_{met}_{eval}_model_eval_boxplots.png", time = time, eval = eval, met = met), 
            base_width = 15, base_height = 13, base_asp = 1.2)
  return(0)
}

eval_metrics <- c("r2", "rmse", "corr")
timepoints <- c("6W", "12M")
for (k in 1:length(eval_metrics)){
  for (j in 1:length(timepoints)){
    plot_comparison_per_metabolite(eval = eval_metrics[k], time = timepoints[j], met = "tar")
  }
}

# Plots of heatmaps by mean of means ####
barplts <- list()
for (k in 1:length(eval_metrics)){
  summary <- readRDS(file = glue("output/analyses/prediction/processed/tar_{eval}_by_met.rds", eval = eval_metrics[k]))
  eval <- map2_df(.x = 1:length(summary), .y = names(summary), .f = function(.x,.y){
    proc <- summary[[.x]] %>% group_by(id) %>% summarise_each(funs(mean), -id2) %>% summarise_each(funs(mean), -id) %>% 
      pivot_longer(everything()) %>% separate(name, c("time", "model")) %>% mutate(met = rep(.y,8))
  })
  
  plt1 <- as.ggplot(pheatmap(pivot_wider(eval %>% filter(time == "6W") %>% select(everything(),-time), names_from = model) 
                             %>% column_to_rownames("met"), color = viridis(100), display_numbers = T, border_color = NA, legend = F,
                             number_color = "white"))
  plt2 <- as.ggplot(pheatmap(pivot_wider(eval %>% filter(time == "12M") %>% select(everything(),-time), names_from = model) 
                             %>% column_to_rownames("met"), color = viridis(100), display_numbers = T, border_color = NA,
                             number_color = "white"))
  g1 <- plot_grid(plt1, plt2, rel_widths = c(1,1.1))
  saveRDS(g1, file = glue("output/figures/prediction/{eval}_tar_heatmap.rds", eval = eval_metrics[k]))
  save_plot(g1, file = glue("output/figures/prediction/{eval}_tar_heatmap.png", eval = eval_metrics[k]), dpi = 300, 
            base_width = 15, base_height = 13, base_asp = 1.2)
  eval2 <- eval %>% group_by(time, model) %>% mutate(eval = rep(eval_metrics[k]))
  # %>% summarise(mean = mean(value)) %>%
  barplts[[k]] <- eval2
}

# boxplots of average means across all metabolites ####
barplts <- bind_rows(barplts)
metrics <- c("Correlation", "R-squared", "RMSE")
names(metrics) <- c("corr", "r2", "rmse")
(barplots <- ggplot(barplts, aes(x = model, y = value, fill = time)) + geom_violin() + geom_boxplot(alpha = 0.8, aes(color = time), show.legend = F) +
  facet_wrap(.~eval, scales = "free",  labeller = labeller(eval = metrics)) + scale_fill_viridis_d() + theme_pubr() + 
  labs(fill = "Time", y = "Mean across all repeats", x = "Models") + 
  scale_x_discrete(labels = c("ENet", "RF", "SPLS", "SVM-RBF")) + theme(axis.text.x = element_text(angle = 90)))
saveRDS(barplots, file="output/figures/prediction/boxplot_across_all_mets.rds")
ggsave(barplots, file = "output/figures/prediction/boxplot_across_all_mets.png", width = 15, height = 13, dpi = 300)

# plotting all borda and getting the ranking matrices ####
for (m in 1:length(eval_metrics)){
  borda_ranking <- matrix(0,nrow = 2, ncol = 4)
  colnames(borda_ranking) <- c("enet", "rf", "svm", "spls")
  rownames(borda_ranking) <- c('6W', "12M")
  rownames(rankings) <- unique(barplts$met)
  colnames(rankings) <- c("enet", "rf", "svm", "spls")
  for (i in 1:length(timepoints)){
    rankings <- matrix(0, nrow = length(unique(barplts$met)), ncol = 4)
    for (j in 1:length(unique(barplts$met))){
      data <- barplts %>% filter(time == timepoints[i], met == unique(barplts$met)[j], eval == eval_metrics[m])
      if (eval_metrics[m] == "rmse"){
        r <- rank(-data$value)
      } else {
        r <- rank(data$value)
      }
      for (k in 1:length(r)){
        borda_ranking[i,r[k]] <- borda_ranking[i,r[k]] + (5 - k)
        rankings[j,r[k]] <- k
      }
    }
    distKendall()
  }
  borda_ranking <- as.data.frame(borda_ranking) %>% rownames_to_column() %>% mutate(time = rowname) %>% select(-rowname) %>%
    pivot_longer(-time, names_to = "model", values_to = "borda_count")
  #%>% pivot_longer(everything(), -rowname, names_to = "model", values_to = "borda_count")
  borda_plt <- ggplot(borda_ranking, aes(x = model, y = borda_count, fill = time)) + geom_bar(stat = "identity", position = "dodge") + 
    scale_fill_viridis_d() + labs(fill = "Time", x = "Models", y = "Borda Score") +
    scale_x_discrete(labels = c("ENet", "RF", "SPLS", "SVM-RBF"))  + theme_pubr() + theme(axis.text.x = element_text(vjust = 0.65)) 
  ggsave(borda_plt, file = glue("output/figures/prediction/borda_plots_tar_{eval}.png", eval = eval_metrics[m]), width = 10, height = 12, dpi = 300)
  saveRDS(borda_plt, file = glue("output/figures/prediction/borda_plots_tar_{eval}.png", eval = eval_metrics[m]))
}


 