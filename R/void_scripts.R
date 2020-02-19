# ggplot nightmare with grids
scca <- readRDS(file = "output/analyses/correlation/12M_tar_scca.rds")
corr <- readRDS(file = "output/analyses/correlation/12M_tar_spearman.rds")
get_plot <- function(dendrogram){
  ddata <- dendro_data(dendrogram, type = "rectangle")
  p <- ggplot(segment(ddata)) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) 
  print(ddata)
  return(p)
}
names(scca) <- c("cca", "boot", "perm", "tab")
cor_mat <- corr$cor_mat
p_mat <- corr$p_mat
tax_dendro <- as.dendrogram(hclust(d = dist(cor_mat)))
met_dendro <- as.dendrogram(hclust(d = dist(t(cor_mat))))

order_met <- order.dendrogram(met_dendro)
order_tax <- order.dendrogram(tax_dendro)
cor.long <- melt(cor_mat, varnames = c("Taxonomy", "Metabolite"))
cor.long$Taxonomy <- factor(x = cor.long$Taxonomy, levels = rownames(cor_mat)[order_tax], ordered = T)
cor.long$Metabolite <- factor(x = cor.long$Metabolite, levels = colnames(cor_mat)[order_met], ordered = T)
heatmap.plot <- ggplot(cor.long, aes(x = Metabolite, y = Taxonomy)) + geom_tile(aes(fill = value)) + 
  scale_fill_viridis(name = "Correlation") + theme_pubr() + theme(legend.position = "none") + theme(axis.title = element_blank())
heatmap.plot
tax_dendro_plot <- get_plot(tax_dendro) + scale_y_reverse() + coord_flip() + 
  theme(axis.line.x = element_blank(), axis.line.y = element_blank())
met_dendro_plot <- get_plot(met_dendro)
met_dendro_plot
grid.newpage()
print(heatmap.plot, vp = viewport(width = 0.8, x = 0.6, y = 0.5, height = 1))
print(tax_dendro_plot, vp = viewport(width = 0.15, x = 0.1, y = 0.51, height = 1.04))


#### boxplot predictions 
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

# FUNCTIONS # ------------------------------------------------------------------------
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

rankings_tot <- list()
# plotting all borda and getting the ranking matrices ####
for (m in 1:length(eval_metrics)){
  borda_ranking <- matrix(0,nrow = 2, ncol = 4)
  colnames(borda_ranking) <- c("enet", "rf", "svm", "spls")
  rownames(borda_ranking) <- c('6W', "12M")
  rankings_tot[[m]] <- list()
  for (i in 1:length(timepoints)){
    rankings <- matrix(0, nrow = length(unique(barplts$met)), ncol = 4)
    rownames(rankings) <- unique(barplts$met)
    colnames(rankings) <- c("enet", "rf", "svm", "spls")
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
    rankings_tot[[m]][[i]] <- rankings
    # getting distances 
    distance <- factoextra::get_dist(rankings, method = "spearman")
    ordinations <- ape::pcoa(distance, correction = "cailliez")
    vec <- as.data.frame(ordinations$vectors)
    ranking_plot <- ggplot(vec, aes(x = Axis.1, y = Axis.2)) + geom_point(color = viridis(100)[40]) + 
      geom_text_repel(aes(label = rownames(rankings))) + theme_pubr() + labs(x = "MDS1", y = "MDS2")
    plot(ranking_plot)
    saveRDS(ranking_plot,
            file = glue("output/figures/prediction/ranking_ordinations_tar_{eval}_{time}.rds", 
                        device = "png", eval = eval_metrics[m], time = timepoints[i]))
    ggsave(ranking_plot, 
           file = glue("output/figures/prediction/ranking_ordinations_tar_{eval}_{time}.png",
                       eval = eval_metrics[m], time = timepoints[i]), device = "png",
           dpi = 300, width = 8, height = 8)
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

rankings_tot[]

# Void simulation scripts 
library(optparse)
library(phyloseq)
library(SpiecEasi)
library(compositions)
library(picante)
library(caret)
library(MLmetrics)

option_list <- list(
  make_option("--calibrate_data", help = "dataset to calibrate to"),
  make_option("--outcome_type", help = "assc for positive simulations, none for unassociated outcomes"),
  make_option("--n_datasets", help = "Number of data sets"),
  make_option("--snr", help = "Signal to noise ratio"),
  make_option("--sample_size", help = "Sample size per data set"),
  make_option("--perc_assc", help = "Proportion of taxa associated with outcome if outcome_type is assc"),
  make_option("--calibrate_outcome", help = "Logical, indiciate whether or not to calibrate unasscoiated outcomes")
)

opt <- parse_args(OptionParser(option_list = option_list))

####-------------------Function definitions--------------------------------------------------------
#' @title Generate simulated data calibrated from real data set  
#' @description Using simulating mechanisms from SpiecEasi package to generate calibrated 
#'   simualted microbiome data set with a log-normal outcome. Outcome can be non-associated
#'   or associated.  
#' @param n_datasets: Number of datasets 
#' @param outcome_type: Outcome type can be "assc" for associative or "none" for unassociated
#' @param calibrate: calibration count data set. Taxa as columns and samples as rows  
#' @param correlation_type: Correlation type for generating data, the default is "cluster"
#' @param snr: Signal to noise ratio  
#' @param sample_size: Sample Size 
#' @param perc_assc: Proportion associated with outcome for "assc" outcome_type 
#' @param beta_val: The value of beta coefficients, default to 6/sqrt(10) and constant across all 
#'   taxa. 
#' @param calibrate_outcome: Outcome data (features as columns and samples as rows) to calibrate 
#'   unassociated outcome. If calibrate_outcome is NULL, use meanlog = 0 and sdlog = 1. 
generate_data <- function(n_datasets, outcome_type, calibrate, correlation_type = "cluster", snr, 
                          sample_size, perc_assc = 0.05, 
                          beta_val = 6/sqrt(10), calibrate_outcome = NULL){
  # generate graph for correlation structure from SpiecEasi  
  graph <- make_graph(method = correlation_type, D = ncol(calibrate), e = ncol(calibrate))
  Prec  <- graph2prec(graph)  
  Cor   <- cov2cor(prec2cov(Prec))
  # use SpiecEasi functions to generate data  
  results <- list()
  depth <- median(rowSums(calibrate))
  calibrate <- unclass(acomp(calibrate))
  calibrate <- round(calibrate * depth)
  for(i in 1:n_datasets){
    sim_tax <- synth_comm_from_counts(calibrate, mar=2, distr='zinegbin', Sigma=Cor, n=sample_size, retParams = F) # faster than doing it by hand  
    # Generate outcomes in naive situation with all taxa having the same regression coefficients  
    ass_idx <- sample(1:ncol(calibrate), perc_assc * ncol(calibrate), replace = F)
    coef <- rep(0, ncol(calibrate))
    coef[ass_idx] <- beta_val
    norm_sim_tax <- unclass(acomp(sim_tax))
    if (outcome_type == "assc"){
      outcome <- 1 + sim_tax %*% coef # identical intercepts 
      error <- rlnorm(sample_size, meanlog = 0, sdlog = 1) # standard log normal error 
      k <- sqrt(var(outcome)/(snr*var(error))) # snr multiplier
      suppressWarnings(new_outcome <- outcome + k * error)
    }
    else {
      if (is.null(calibrate_outcome)){
        outcome = rlnorm(sample_size, meanlog = 0, sdlog = 1)
        error <- rlnorm(sample_size, meanlog = 0, sdlog = 1) # standard log normal error 
        k <- sqrt(var(outcome)/(snr*var(error))) # snr multiplier
        suppressWarnings(new_outcome <- outcome + k * error)
      } else {
        fit <- fitdist(calibrate_outcome[calibrate_outcome > 0], distr = "lnorm", method = "mle")
        outcome = rlnorm(sample_size, meanlog = fit$estimate[1], sdlog = fit$estimate[2])
        error <- rlnorm(sample_size, meanlog = 0, sdlog = 1) # standard log normal error 
        k <- sqrt(var(outcome)/(snr*var(error))) # snr multiplier
        suppressWarnings(new_outcome <- outcome + k * error)
      }
    }
    results[[i]] <- list(tax = sim_tax, 
                         coef = coef, 
                         norm_sim_tax = norm_sim_tax, 
                         org_outcome = outcome, 
                         met = new_outcome)
  }
  return(results)
}

calibrate_data <- readRDS(file = opt$calibrate_data)
tax <- calibrate_data$tax
met <- calibrate_data$met

if (opt$calibrate_outcome == T){
  if (opt$outcome_type == "assc"){
    stop("calibrate outcome can't be TRUE if outcome type is associative")
  } else {
    calibrate_outcome = as.vector(met)
  }
} else {
  calibrate_outcome <- NULL
}

datasets <- generate_data(n_datasets = opt$n_datasets, perc_assc = opt$perc_assc, outcome_type = opt$outcome_type, calibrate = tax, 
                          snr = opt$snr, sample_size = opt$sample_size, calibrate_outcome = calibrate_outcome)

output <- list(datasets = datasets, params = list(
  n_datasets = opt$n_datasets, 
  perc_assc = opt$perc_assc,
  outcome_type = opt$outcome_type, 
  snr = opt$snr, 
  sample_size = opt$sample_size,
  calibrate_outcome = opt$calibrate_outcome
))

saveRDS(output, file = paste0("snakemake_output/simulation/", opt$outcome_type,"_",opt$snr,"_",opt$perc_assc,"_overall_dataset.rds"))
sapply(1:opt$n_datasets, function(x){
  obj <- datasets[[x]]
  path = paste0("snakemake_output/simulation/single_obj/", opt$outcome_type, "/",opt$snr,"/",opt$perc_assc, "/")
  if(!dir.exists(path)){
    dir.create(path, recursive = T)
  }
  saveRDS(obj, file = paste0(path, "simulated_", x, ".rds"))
})




