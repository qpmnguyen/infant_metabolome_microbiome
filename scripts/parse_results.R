library(ggplot2)
library(dplyr)
library(purrr)
library(MLmetrics)
library(reshape2)
library(gridExtra)


files <- list.files("./data/12M_untar_clr_tax/")
met_names <- readRDS(file = "./data/tarNMR_met_names.rds")
null <- readRDS(file = "./data/null_permutation.rds")
null <- do.call(rbind, null)
colnames(null) <- c("RRSE", "RMSE", "Correlation", "R-squared")

RMSE <- RRSE <- Corr <- R2 <- matrix(ncol = 208, nrow = 10)
for (i in 1:length(files)){
  dir <- paste0("./data/12M_untar_clr_tax/",files[i])
  res <- readRDS(dir)
  idx <- as.numeric(strsplit(strsplit(files[i],"_")[[1]][6], ".", fixed = T)[[1]][1])
  rmse <- map_dbl(res, function(.x){
    MLmetrics::RMSE(y_pred = .x[,1], y_true = .x[,2])
  })
  rrse <- map_dbl(res, function(.x){
    RRSE(y_pred = .x[,1], y_true = .x[,2])
  })
  corr <- map_dbl(res, function(.x){
    cor(x = .x[,1], y = .x[,2], method = "spearman")
  })
  r2 <- map_dbl(res, function(.x){
    R2_Score(y_pred = .x[,1], y_true = .x[,2])
  })
  RMSE[,idx] <- rmse
  RRSE[,idx] <- rrse
  Corr[,idx] <- corr
  R2[,idx] <- r2
}

colnames(RMSE) <- colnames(RRSE) <- colnames(Corr) <- colnames(R2) <- paste0("B", 1:208)

generate_mean_barplot <- function(df, ylab, xlab, type, datatype = NULL){
  df <- melt(df) %>% group_by(Var2) %>% summarise(avg = mean(value, trim = 0.1)) 
  plt <- ggplot(df, aes(x = reorder(Var2, -avg), y = avg)) + 
    geom_bar(stat = "identity", color = "black", fill = "blue", alpha = 0.5) + 
    labs(y = ylab, x = xlab, title = type) + coord_flip() + 
    theme(legend.position = "none") + theme_bw() 
  if (datatype == "binned"){
    plt <- plt + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  } else {
    plt <- plt + theme(axis.text.x = element_text(size = 9)) + geom_hline(yintercept = median(null[,type]), col = "red")
  }
  return(plt)
}

p1 <- generate_mean_barplot(RMSE, xlab = "Metabolite", ylab = "Trimmed Mean RMSE", type = "RMSE", datatype = "binned")
p2 <- generate_mean_barplot(R2, xlab = "Metabolite", ylab = "Trimmed Mean R-squared", type = "R-squared", datatype = "binned")
p3 <- generate_mean_barplot(RRSE, xlab = "Metabolite", ylab = "Trimmed Mean RRSE", type = "RRSE", datatype = "binned") + geom_hline(yintercept = 1, col = "red")
p4 <- generate_mean_barplot(Corr, xlab = "Metabolite", ylab = "Trimmed Mean Spearman Cor.", type = "Correlation", datatype = "binned") + geom_hline(yintercept = 0.3, col = "red")


joined <- grid.arrange(p1,p2,p3,p4, ncol = 2, nrow = 2)
plot(joined)

ggsave(filename = "./docs/rf_untarNMR.png", joined, device = "png", width = 15, height = 15, units = "in")






