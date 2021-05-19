# This script seeks to perform different types of plotting for bayesian analyses  
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
library(ggrepel)
library(openxlsx)
library(patchwork)
library(ggsci)

summary_tar <- readRDS(file = "output/analyses/prediction/processed/bayesian_models_summary.rds")
summary_untar <- readRDS(file = "output/analyses/prediction/processed/bayesian_models_untar_summary.rds")
eval <- c("r2", "corr")
time <- c("6W", "12M")
get_forestplot <- function(summary, time, eval){
  df_6W <- bind_rows(summary_tar[[eval]][['6W']],.id = "Metabolite") %>% mutate(time = rep('6W', n())) %>% mutate(model = gsub("6W_","", model))
  df_12M <- bind_rows(summary_tar[[eval]][['12M']], .id = "Metabolite") %>% mutate(time = rep('12M', n())) %>% mutate(model = gsub("12M_","", model))
  df <- bind_rows(df_6W, df_12M)
  plot <- ggplot(df, aes(y = mean, x = model, ymin = lower, ymax = upper, color = time)) + 
    geom_point(position = position_dodge(width = 0.9), size = 1) + geom_hline(yintercept = 0, col = "red", linetype = 2) + 
    geom_errorbar(position = "dodge") + coord_flip() + 
    theme_bw(base_size = 20)  + 
    scale_color_npg() + 
    theme(panel.spacing.x = unit(1.5, "lines")) + labs(color = "Time") + 
    facet_wrap(~Metabolite) + labs(x = "Models", y = "Posterior Mean (95% Credible Interval)")
  return(plot)
}
get_violinplot <- function(summary, eval){
  if (eval == "r2"){
    title <- "R-squared"
  } else {
    title <- "Spearman correlation"
  }
  df1 <- bind_rows(summary[[eval]][["6W"]], .id = "Metabolite") %>% mutate(model = gsub(paste0("6W","_"),"",model), time = "6W", metric = eval)
  df2 <- bind_rows(summary[[eval]][["12M"]], .id = "Metabolite") %>% mutate(model = gsub(paste0("12M","_"),"",model), time = "12M", metric = eval)
  df <- rbind(df1, df2)
  plt <- ggplot(df, aes(y = mean, x = model, fill = time)) + geom_violin() + geom_boxplot(alpha = 0.8, aes(col = time), show.legend = F) +  
    scale_fill_viridis_d() + theme_pubr() + labs(y = "Posterior Mean", x = "Model", title = title, fill = "Time") + 
    theme(plot.title = element_text(face = "bold"))
  return(plt)
}
get_bordaplot <- function(summary, eval){
  df1 <- bind_rows(summary[[eval]][["6W"]], .id = "Metabolite") %>% mutate(model = gsub(paste0("6W","_"),"",model), time = "6W", metric = eval)
  df2 <- bind_rows(summary[[eval]][["12M"]], .id = "Metabolite") %>% mutate(model = gsub(paste0("12M","_"),"",model), time = "12M", metric = eval)
  df <- rbind(df1, df2)
  borda_ranking <- matrix(0,nrow = 2, ncol = 4)
  colnames(borda_ranking) <- c("enet", "rf", "svm", "spls")
  rownames(borda_ranking) <- c('6W', "12M")
  for (i in 1:length(c("6W", "12M"))){
    ranking <- rep(0, 4)
    for (j in 1:length(unique(df$Metabolite))){
      tbl <- df %>% filter(time == c("6W", "12M")[i], Metabolite == unique(df$Metabolite)[j])
      r <- rank(tbl$mean)
      for (k in 1:length(r)){
        ranking[r[k]] <- ranking[r[k]] + (5 - k)
      }
    }
    print(ranking)
    borda_ranking[i,] <- ranking
  }
  borda_ranking <- as.data.frame(borda_ranking) %>% rownames_to_column() %>% mutate(time = rowname) %>% select(-rowname) %>%
    pivot_longer(-time, names_to = "model", values_to = "borda_count")
  borda_plt <- ggplot(borda_ranking, aes(x = model, y = borda_count, fill = time)) + geom_bar(stat = "identity", position = "dodge") + 
    scale_fill_npg() + labs(fill = "Time", x = "Model", y = "Borda Score") +
    theme_bw() + theme(axis.text.x = element_text(vjust = 0.65)) 
  return(borda_plt)
}

# get plots 
for (i in eval){
  plt <- get_forestplot(summary = summary_tar, eval = i)
  saveRDS(plt, file = glue("output/figures/prediction/{i}_tar_forestplot.rds"))
}


for (i in eval){
  plt <- readRDS(file = glue("output/figures/prediction/{i}_tar_forestplot.rds"))
  ggsave(plt, filename = glue("output/figures/prediction/{i}_tar_forestplot.png"), dpi = 300,
         width = 14, height = 12)
}


for (k in c("tar", "untar")){
  for (i in eval){
    plt <- get_violinplot(summary = get(glue("summary_{type}", type = k)), i)
    saveRDS(plt, file = glue("output/figures/prediction/{i}_{type}_violinplots.rds", i = i, type = k))
  }
  for (i in eval){
    plt <- get_bordaplot(summary = get(glue("summary_{type}", type = k)), i)
    saveRDS(plt, file = glue("output/figures/prediction/{i}_{type}_bordaplots.rds", i = i, type = k))
  }
}


# export it out as files 
export <- list("6W_targeted_corr" = c(), "6W_targeted_r2"=c(), "12M_targeted_corr" = c(), "12M_targeted_r2" = c(), 
               "6W_untargeted_corr" = c(), "6W_untargeted_r2"=c(), "12M_untargeted_corr" = c(), "12M_untargeted_r2" = c())

for (i in time){
  for(j in eval){
    for (k in c("targeted", "untargeted")){
      if (k == "targeted"){
        sum <- summary_tar
      } else {
        sum <- summary_untar
      }
      df <- bind_rows(sum[[j]][[i]],.id = "Metabolite") %>% 
        mutate(time = rep(i, n())) %>% 
        mutate(model = gsub(paste0(time,"_"),"", model))
      df <- df %>% mutate_if(is.numeric, round, 3) %>% 
        mutate("Posterior mean (95% CI)" = paste0(mean,"(",lower,",",upper,")")) %>% 
        select(-mean,-lower,-upper,-time) %>%
        pivot_wider(names_from = model, values_from = "Posterior mean (95% CI)")
      export[[glue("{time}_{type}_{eval}", time = i, type = k, eval = j)]] <- df
    }
  }
}
export
write.xlsx(export, file = "docs/Supplementary_file_1.xlsx", asTable = T)


# Combine forest plots 
tar_plts <- vector(mode = "list", length = 2)
for (i in seq_along(eval)){
  tar_plts[[i]] <- get_forestplot(summary = summary_tar, time = "6W", eval = eval[i])
}
tar_plts[[1]] <- tar_plts[[1]] + theme(legend.position = "top", axis.title.x = element_blank()) 
tar_plts[[2]] <- tar_plts[[2]] + theme(legend.position = "none")
tar_forestplot <- tar_plts[[1]]/tar_plts[[2]] + plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold"))

ggsave(tar_forestplot, filename = "docs/publication_figures/tar_forestplot.png", dpi = 300, 
       width = 17, height = 18)

untar_plts <- vector(mode = "list", length = 2)
for (i in seq_along(eval)){
  untar_plts[[i]] <- get_forestplot(summary = summary_untar, time = "6W", eval = eval[i])
}
untar_plts[[1]] <- untar_plts[[1]] + theme(legend.position = "top", axis.title.x = element_blank()) 
untar_plts[[2]] <- untar_plts[[2]] + theme(legend.position = "none")
untar_forestplot <- untar_plts[[1]]/untar_plts[[2]] + plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold"))

ggsave(untar_forestplot, filename = "docs/publication_figures/untar_forestplot.png", dpi = 300, 
       width = 17, height = 18)
