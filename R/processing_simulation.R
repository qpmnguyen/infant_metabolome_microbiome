# Script to plot simple evaluation of positive and negative simulations
# Plotting is incorporated 
library(glue)
library(tidyverse)
library(ggpubr)
library(viridis)
library(cowplot)

pos_results <- readRDS(file = "output/analyses/simulation/12M_positive_simulation_results.rds")
neg_summaries <- readRDS(file = "output/analyses/simulation/12M_negative_simulation_results.rds")
param <- pos_results$grid
pos_summaries <- pos_results$summary

# defining some functions
get_eval <- function(boot_list, eval, type){
  if (eval == "r2"){
    list <- map(boot_list, function(.x){
      if (length(.x$r2) > 4){
        r2 <- as.vector(unlist(.x$r2[c(3,6,9,12)])) # did not think through yardstick properly
      } else {
        r2 <- as.vector(.x$r2)
      }
      r2[which(is.na(r2))] <- 0
      return(r2)
    })
  } else if (eval == "corr"){
    list <- map(boot_list, function(.x){
      corr <- as.vector(.x$corr)
      corr[which(is.na(corr))] <- 0
      return(corr)
    })
  }
  if (type == "pos"){
    names(list) <- paste0("Bootstrap",1:100)
    df <- bind_rows(list) %>% mutate(model = c("enet", "rf", "spls", "svm")) %>% pivot_longer(-model) %>% mutate(eval = rep(eval, 400))
  } else {
    names(list) <- paste0("Permutation",1:1000)
    df <- bind_rows(list) %>% mutate(model = c("enet", "rf", "spls", "svm")) %>% pivot_longer(-model) %>% mutate(eval = rep(eval, 4000))
  }
  return(df)
}


total_pos <- map2_df(pos_summaries, 1:length(pos_summaries), function(.x, .y){
  corr <- get_eval(.x, "corr", "pos")
  r2 <- get_eval(.x, "r2", "pos")
  joined <- full_join(corr, r2) %>% mutate(snr = rep(param$snr[.y],800), spar = rep(param$spar[.y],800))
  return(joined)
})

total_neg <- full_join(get_eval(neg_summaries, "r2", "neg"), get_eval(neg_summaries, "corr", "neg"))

# plotting
total_pos$snr <- as.factor(total_pos$snr)
total_pos$spar <- as.factor(total_pos$spar)
neg_filtered <- total_neg %>% group_by(model, eval) %>% summarise(mean = mean(value[is.infinite(value) != TRUE]), 
                                                                  lower = quantile(value[is.infinite(value) != TRUE], 0.05), 
                                                                  upper = quantile(value[is.infinite(value) != TRUE], 0.95))
label <- c("SCC", "R2")
names(label) <- c("corr", "r2")

pos_plt <- ggplot(total_pos, aes(x = model, y = value, fill = eval)) + geom_violin() + geom_boxplot(alpha = 0.5, aes(col = eval), show.legend = F) + 
  facet_grid(snr ~ spar, scales = "free_y", labeller = label_both) + theme_pubr() + theme(axis.text.x = element_text(angle = 25)) +
  scale_fill_viridis_d(name = "Evaluation Metric", labels = c("SCC", "R2")) + labs(y = "Evaluation", x = "Model") 

neg_plt <- ggplot(neg_filtered, aes(x = model, y = mean, ymin = lower, ymax = upper)) + geom_errorbar(width = 0.2) + 
  geom_pointrange(aes(col = model), show.legend = F) + theme_pubr() + 
  coord_flip() + facet_wrap(eval~., scales = "free_x", labeller = labeller(eval = label)) + xlab("Models") + ylab("Evaluation (95% quantile)")

grid_plot <- plot_grid(pos_plt, neg_plt, labels = c("A. Positive", "B. Negative"), scale = 0.95)
save_plot(grid_plot, file = "docs/publication_figures/simulation_cobbplot.png", base_width = 14, base_height = 10)
