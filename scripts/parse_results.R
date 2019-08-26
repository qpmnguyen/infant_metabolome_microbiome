library(ggplot2)
library(dplyr)
library(purrr)
library(MLmetrics)
library(reshape2)
library(gridExtra)
library(kableExtra)
library(ggthemes)
library(gtable)
library(grid)
parse_results <- function(method, time, data_type){
  directory <- paste0("./data/", method, "/", time, "_", data_type, "_clr_tax/")
  files <- list.files(directory)
  met_names <- readRDS(file = "./data/tarNMR_met_names.rds")
  if (data_type == "tar"){
    RMSE <- RRSE <- Corr <- R2 <- matrix(ncol = 36, nrow = 10)
    colnames(RMSE) <- colnames(RRSE) <- colnames(Corr) <- colnames(R2) <- met_names
  } else if (data_type == "untar"){
    RMSE <- RRSE <- Corr <- R2 <- matrix(ncol = 208, nrow = 10)
    colnames(RMSE) <- colnames(RRSE) <- colnames(Corr) <- colnames(R2) <- paste0("B", 1:208)
  }
  for (i in 1:length(files)){
    dir <- paste0(directory,files[i])
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
  return(list(RMSE = RMSE, RRSE = RRSE, Correlation = Corr, R2 = R2))
}

generate_mean_barplot <- function(df, ylab, xlab, type, datatype = NULL){
  null <- readRDS(file = "./data/null_permutation.rds")
  null <- do.call(rbind, null)
  colnames(null) <- c("RRSE", "RMSE", "Correlation", "R2")
  df <- melt(df) %>% group_by(Var2) %>% summarise(avg = mean(value, trim = 0.1)) 
  plt <- ggplot(df, aes(x = reorder(Var2, -avg), y = avg)) + 
    geom_bar(stat = "identity", color = "black", fill = "blue", alpha = 0.5) + 
    labs(y = ylab, x = xlab, title = type) + coord_flip() + 
    theme(legend.position = "none") + theme_bw() 
  if (datatype == "untar"){
    plt <- plt + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  } else if (datatype == "tar"){
    plt <- plt + theme(axis.text.x = element_text(size = 9)) + geom_hline(yintercept = median(null[,type]), col = "red")
  }
  return(plt)
}

results <- parse_results("rf", "12M", "untar") 

# making joint plots  
plots <- map2(names(results), results, function(.x, .y){
  p <- generate_mean_barplot(.y, ylab = "Trimmed 10% Mean", xlab = "Metabolites", type = .x, datatype = "untar")
  return(p)
})


joined <- grid.arrange(grobs = plots, ncol = 2, nrow =2)
#ggsave(filename = "./docs/12M_rf_untarNMR.png", joined, device = "png", width = 15, height = 15, units = "in")


# make tables 
grid <- expand.grid(c("12M", "6W"),c("tar", "untar"))

# overall combined results 
generate_boxplots <- function(list){
  df <- map(list, function(.x){
    apply(.x, 2, mean)
  })
  df <- melt(df)
  ggplot(df, aes(x = L1, y = value)) + geom_boxplot() + coord_cartesian(ylim = c(-2.5,2))
}


 # all grid results
combined <- apply(grid, 1, function(x){
  res <- parse_results('rf', x[1], x[2])
  df <- melt(map(res, function(.x){
    apply(.x, 2, mean)
  }))
  n_metab <- nrow(df)/4
  df['Metabolite'] <- rep(paste0("M",1:n_metab),4)
  df['time'] <- rep(x[1], nrow(df))
  df['data_type'] <- rep(x[2], nrow(df))
  return(df)
})

combined <- as.data.frame(do.call(rbind, combined))

comb_filt <- combined %>% filter(L1 != "RMSE")
plt <- ggplot(comb_filt, aes(x = L1, y = value)) + geom_boxplot(aes(fill = L1)) + 
  facet_grid(rows = vars(time), cols = vars(data_type)) + coord_cartesian(ylim = c(-1,1.5)) + labs(y = "Value", x = "Evaluation Measure") +
  theme(legend.position = "None")
#ggsave(plot = plt, filename = "./docs/comparison_boxplots.png", device = "png", width = 10, height = 9, units = "in")
plt
# making maximum correlation table 
dat <- combined %>% group_by(time, data_type) %>% filter(L1 == "Correlation") %>% summarise(max = max(value), max_metab = Metabolite[which.max(value)])

table <- do.call(rbind,apply(dat, 1, function(x){
  combined %>% group_by(time, data_type, L1) %>% filter(time == x[1], data_type == x[2], Metabolite == x[4])
}))

table$data_type[table$data_type == "tar"] <- "Targeted"
table$data_type[table$data_type == "untar"] <- "Untargeted"

table <- table[,c(4,5,3,2,1)]
colnames(table) <- c("Time", "Metabolite Type", "Metabolite Id", "Evaluation Metric", "Value")


kable(table, "html") %>% kable_styling("striped") %>% save_kable("./docs/test.png")


# generating pair plots of all folds
fold_pair_plots <- function(method, time, data_type){
  directory <- paste0("./data/", method, "/", time, "_", data_type, "_clr_tax/")
  files <- list.files(directory)
  if (data_type == "tar"){
    met_names <- readRDS(file = "./data/tarNMR_met_names.rds")
  } else if (data_type == "untar"){
    met_names <- paste0("B", 1:208)
  }
  plot_list <- list()
  for (i in 1:length(files)){
    dir <- paste0(directory, files[i])
    idx <- as.numeric(strsplit(strsplit(files[i], "_")[[1]][6], ".", fixed = T)[[1]][1])
    res <- readRDS(dir)
    new_res <- map2(1:length(res), res, function(.x, .y){
      folds <- rep(.x, nrow(.y))
      return(cbind(.y,folds))
    })
    new_res <- as.data.frame(do.call(rbind, new_res))
    new_res$folds <- as.factor(new_res$folds)
    plt <- ggplot(new_res, aes(x = predictions, y = met.test)) + geom_point(aes(col = folds)) + 
      geom_smooth(method = "lm", aes(col = folds)) + scale_color_stata() + labs(title = met_names[idx]) + ylab(NULL) + xlab(NULL) +
      theme(title = element_text(size = 9))
    plot_list[[i]] <- plt
  }
  return(plot_list)
}

for (i in 1:nrow(grid)){
  plot_list <- fold_pair_plots(method = "rf", time = grid$Var1[i], data_type = grid$Var2[i])
  legend <- gtable_filter(ggplotGrob(plot_list[[1]]), "guide-box")
  plot_list <- map(plot_list, function(.x){
    .x + theme(legend.position = "None")
  })
  grid_arr <- grid.arrange(arrangeGrob(grobs = plot_list,
                                       left = textGrob("True", rot = 90, vjust = 1),
                                       bottom = textGrob("Predicted", hjust = 0.5)),
                                       legend, widths = unit.c(unit(1, "npc") - legend$widths, legend$widths),
                                       nrow = 1)
  name <- paste0(grid$Var1[i],"_", grid$Var2[i], "_pair_plot.png")
  ggsave(plot = grid_arr, filename = name, device = "png", width = 15, height = 10, units = "in")
}


