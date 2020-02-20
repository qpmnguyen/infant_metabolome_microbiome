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


summary <- readRDS(file = "output/analyses/prediction/processed/bayesian_models_summary.rds")
eval <- c("r2", "corr")
time <- c("6W", "12M")
get_forestplot <- function(time, eval){
  df <- bind_rows(summary[[eval]][[time]], .id = "Metabolite") %>% mutate(model = gsub(paste0(time,"_"),"",model), time = time, metric = eval)
  plot <- ggplot(df, aes(y = model, x = mean, xmin = lower, xmax = upper)) + geom_point() + geom_errorbarh() + 
    facet_wrap(~Metabolite) + theme_pubr() + labs(y = "Models", x = "95% Credible Interval")
  if (eval == "r2" | eval == "corr"){
    plot <- plot + geom_vline(xintercept = 0, col = "red", linetype = 2)
  }
  return(plot)
}
get_violinplot <- function(eval){
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
get_bordaplot <- function(eval){
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
    scale_fill_viridis_d() + labs(fill = "Time", x = "Model", y = "Borda Score") +
    theme_pubr() + theme(axis.text.x = element_text(vjust = 0.65))
  return(borda_plt)
}

# get plots 
for (i in eval){
  for (j in time){
    plt <- get_forestplot(time = j, eval = i)
    saveRDS(plt, file = glue("output/figures/prediction/{j}_{i}_forestplot.rds"))
  }
}
for (i in eval){
  plt <- get_violinplot(i)
  saveRDS(plt, file = glue("output/figures/prediction/{i}_tar_violinplots.rds"))
}
for (i in eval){
  plt <- get_bordaplot(i)
  saveRDS(plt, file = glue("output/figures/prediction/{i}_tar_bordaplots.rds"))
}


