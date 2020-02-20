library(glue)
library(tidyverse)
library(ggpubr)
library(viridis)
library(cowplot)

results <- readRDS(file = "output/analyses/simulation/12M_positive_simulation_results.rds")
param <- results$grid
summaries <- results$summary

# defining some functions
get_r2 <- function(boot_list){
  r2_list <- map(boot_list, function(.x){
    r2 <- as.vector(unlist(.x$r2[c(3,6,9,12)])) # did not think through yardstick properly
    return(r2)
  })
  names(r2_list) <- paste0("Bootstrap",1:100)
  r2_df <- bind_rows(r2_list) %>% pivot_longer(everything()) %>% mutate(eval = rep("r2",100))
}

get_eval <- function(boot_list, eval){
  if (eval == "r2"){
    list <- map(boot_list, function(.x){
      r2 <- as.vector(unlist(.x$r2[c(3,6,9,12)])) # did not think through yardstick properly
      return(r2)
    })
  } else if (eval == "corr"){
    list <- map(boot_list, function(.x){
      corr <- as.vector(.x$corr)
      return(corr)
    })
  }
  names(list) <- paste0("Bootstrap",1:100)
  df <- bind_rows(list) %>% mutate(model = c("enet", "rf", "spls", "svm")) %>% pivot_longer(-model) %>% mutate(eval = rep(eval, 400))
  return(df)
}


total <- map2_df(summaries, 1:length(summaries), function(.x, .y){
  corr <- get_eval(.x, "corr")
  r2 <- get_eval(.x, "r2")
  joined <- full_join(corr, r2) %>% mutate(snr = rep(param$snr[.y],800), spar = rep(param$spar[.y],800))
  return(joined)
})
total$snr <- as.factor(total$snr)
total$spar <- as.factor(total$spar)

plt <- ggplot(total, aes(x = model, y = value, fill = eval)) + geom_violin() + geom_boxplot(alpha = 0.5, aes(col = eval), show.legend = F) + 
  facet_grid(snr ~ spar, scales = "free_y") + theme_cleveland() + 
  scale_fill_viridis_d(name = "Evaluation Metric", labels = c("SCC", "R2")) + labs(y = "Value", x = "Model") 

