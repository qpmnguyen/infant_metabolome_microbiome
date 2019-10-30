library(ggplot2)
library(ggpubr)
library(viridis)
library(reshape2)
library(dplyr)
library(cowplot)

# getting everything to the correct format 
one_yr <- readRDS(file = "snakemake_output/analyses/prediction_evaluation/16S_12M_tar_rf.rds")
six_wks <- readRDS(file = "snakemake_output/analyses/prediction_evaluation/16S_6W_tar_rf.rds")

cor_12M <- melt(one_yr$correlation) # extract relvant metrics
cor_12M$Var1 <- rep("12M", nrow(cor_12M)) # add time labels 
cor_12M$Var3 <- rep("Pearson Correlation", nrow(cor_12M)) # add correlation labels 
cor_6W <- melt(six_wks$correlation)
cor_6W$Var1 <- rep("6W", nrow(cor_6W))
cor_6W$Var3 <- rep("Pearson Correlation", nrow(cor_6W)) # add correlation labels 

r2_12M <- melt(one_yr$r2)
r2_12M$Var1 <- rep("12M", nrow(r2_12M))
r2_12M$Var3 <- rep("Pseudo R-squared", nrow(r2_12M)) # add correlation labels 
r2_6W <- melt(six_wks$r2)
r2_6W$Var1 <- rep("6W", nrow(r2_6W))
r2_6W$Var3 <- rep("Pseudo R-squared", nrow(r2_6W)) # add correlation labels 


# Plot by time point --------------------------------------------------
combined_12M <- rbind(cor_12M, r2_12M)
mean_12M_corr <- cor_12M %>% group_by(Var1, Var2) %>% summarise(mean = mean(value, trim = 0.1))

ord_names <- mean_12M_corr$Var2[order(mean_12M_corr$mean, decreasing = F)]
plt_12M <- ggboxplot(combined_12M, x = "Var2", y = "value", orientation = "horizontal", add= "jitter", facet.by = "Var3", 
          fill = viridis(20, alpha = 0.8)[9], order = ord_names) + geom_hline(yintercept = 0, color = 'red')
plt_12M <- ggpar(plt_12M, xlab = "Metabolites", ylab = "12 Months")
plot(plt_12M)

combined_6W <- rbind(cor_6W, r2_6W)
mean_6W_corr <- cor_6W %>% group_by(Var1, Var2) %>% summarise(mean = mean(value, trim = 0.1))

ord_names <- mean_6W_corr$Var2[order(mean_6W_corr$mean, decreasing = F)]
plt_6W <- ggboxplot(combined_6W, x = "Var2", y = "value", orientation = "horizontal", add= "jitter", facet.by = "Var3", 
                     fill = viridis(20, alpha = 0.8)[9], order = ord_names) + geom_hline(yintercept = 0, color = 'red')
plt_6W <- ggpar(plt_6W, xlab = "Metabolites", ylab = "6 Weeks")
plot(plt_6W)

## Combined correlation ------------------------------------------------------------------------
combined_corr <- rbind(cor_12M, cor_6W)
trimmed_mean_corr <- combined_corr %>% group_by(Var1, Var2) %>% summarise(mean = mean(value, trim = 0.1))
trimmed_mean_corr[order(trimmed_mean_corr$mean),]


cor_plt <- ggboxplot(combined_corr, x = "Var2", y = "value", orientation = "horizontal",
                 add = "jitter", fill = viridis(20)[8])
cor_plt <- ggpar(cor_plt, legend = "none", xlab = "Metabolite", ylab = "Pearson Correlation") + geom_hline(yintercept = 0, color = 'red')
cor_plt <- facet(cor_plt, facet.by = "Var1")
plot(cor_plt)

### R-squared -------------------------------------------------------------------------------
combined_r2 <- rbind(r2_12M, r2_6W)

r2_plt <- ggboxplot(combined_r2, x = "Var2", y = "value", fill = viridis(20)[8], orientation = "horizontal",
                     add = "jitter")
r2_plt <- ggpar(r2_plt, legend = "none", xlab = "Metabolite", ylab = "Pseudo R-squared") + geom_hline(yintercept = 0, color = 'red')
r2_plt <- facet(r2_plt, facet.by = "Var1")
plot(r2_plt)

### best and correlation -----------------------------------------------
max <- trimmed_mean_corr[trimmed_mean_corr$mean == max(trimmed_mean_corr$mean),]
max_frame <- six_wks$raw[[max$Var2]]

max_frame <- lapply(1:length(max_frame), function(i){
  max_frame[[i]] <- as.data.frame(max_frame[[i]])
  max_frame[[i]]$Fold <- rep(paste0("Fold ",i), nrow(max_frame[[i]]))
  return(max_frame[[i]])
})

max_frame <- do.call(rbind, max_frame)

max_frame <- max_frame %>% mutate(residuals = predictions - met.test)
max_plot <- ggscatter(max_frame, x = "predictions", y = "met.test", color = "Fold", palette = viridis(5), add = "reg.line", conf.int = T)
max_plot <- ggpar(max_plot, xlab = "Predicted (log10)", ylab = "True Values (log10)", main = "Butyrate at 6 weeks - Trimmed mean correlation = 0.519",
                  xscale = "log10", yscale = "log10") 

### Average correlation between time points --------------------------------
combined <- rbind(combined_12M, combined_6W)
library(nlme)
mod_cor <- lme(value ~ Var1, random = ~1|Var2, data = combined[combined$Var3 == "Pearson Correlation",], method = "REML")
anova.lme(mod_cor)
mod_r2 <- lme(value ~ Var1, random = ~1|Var2, data = combined[combined$Var3 == "Pseudo R-squared",], method = "REML")
anova.lme(mod_r2)
lab <- data.frame(label = c('Nested ANOVA: p = 2e-4', 'Nested ANOVA: p = 0.338'), Var3 = c("Pearson Correlation", "Pseudo R-squared"), pos = c(-0.4, 0.6))
overall_plot <- ggboxplot(combined, x = "Var1", y = "value", facet.by = "Var3", xlab = "Time Points", ylab = FALSE, fill = viridis(20)[8]) + geom_text(data = lab, aes(x = 1.5, y = pos, label = label))


grid <- plot_grid(plt_6W, 
          plt_12M, 
          max_plot, 
          overall_plot, 
          nrow = 2, labels = c("A", "B", "C", "D"))
plot(grid)


ggsave(grid, file = "snakemake_output/figures/predictions/prediction_result_figures.png", 
       device = "png", width = 15, height = 15)

### Model comparison boxplots ------------------------------------------------------------------






