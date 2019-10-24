<<<<<<< HEAD
# Updated 09/12/19
# Quang Nguyen
# Script to generate singular and joint ordination plots after ordnation analyses 

library(optparse)
library(ggplot2)
library(ggpubr)
library(vegan)

option_list <- list(
  make_option("--input", help = "Input file for data loading"),
  make_option("--time", help = "Data time point to extract"),
  make_option("--metab_type", help = "Data type for metabolite data, can be 'tar' or 'untar'"),
  make_option("--tax_type", help = "Data type for taxonomic data, so far only '16S' is supported")
)

opt <- parse_args(OptionParser(option_list = option_list))

output_folder <- paste0("figures/ordinations/",opt$tax_type, "/", opt$time, "/", opt$metab_type, "/")
data <- readRDS(file = opt$input)

tax_ord <- data$tax_ord
met_ord <- data$met_ord
proc_test <- data$proc_test
mant_test <- data$mant_test 

# Single data ordinations
met_pts <- as.data.frame(met_ord$points)
tax_pts <- as.data.frame(tax_ord$points)
colnames(met_pts) <- colnames(tax_pts) <- c("NMDS1", "NMDS2")
met_plt <- ggscatter(data =met_pts, x = "NMDS1", y = "NMDS2", color = "coral", ellipse = T, ellipse.type = "t", star.plot = T, mean.point = T, 
                    ellipse.alpha = 0, mean.point.size = 3)
met_plt <- ggpar(met_plt, title = "Metabolites")
tax_plt <- ggscatter(data = tax_pts, x = "NMDS1", y = "NMDS2", color = "steelblue", ellipse = T, ellipse.type = "t", star.plot = T, mean.point = T, 
                        ellipse.alpha = 0)
tax_plt <- ggpar(tax_plt, title = "Taxonomy")


saveRDS(tax_plt, file = paste0(output_folder, "tax_ordinations.rds"))
ggsave(filename = paste0(output_folder, "tax_ordinations.svg"), plot = tax_plt, device = "svg")
saveRDS(met_plt, file = paste0(output_folder, "met_ordinations.rds"))
ggsave(filename = paste0(output_folder, "met_ordinations.svg"), plot = met_plt, device = "svg")
    
# joint data set ordination 
plot_dat <- data.frame(rbind(proc_test$Yrot, proc_test$X))
plot_dat[,3] <- c(rep("Metabolite",nrow(proc_test$Yrot)), rep("Taxonomy", nrow(proc_test$X)))
colnames(plot_dat) <- c("NMDS1", "NMDS2", "Ordination")
plt <- ggscatter(plot_dat, x= "NMDS1", y = "NMDS2", 
            color = "Ordination",
            palette = "jco", 
            size = 3, 
            star.plot.lwd = 1,
            ellipse = T,
            ellipse.type = "t",
            star.plot = F, mean.point = F,
            repel = T)
plt <- ggpar(plt,subtitle = paste("Procrustes SoS:", round(proc_test$ss,4), "; Sig:", round(proc_test$signif,4)))

saveRDS(plt, file = paste0(output_folder, "joint_protest_ordinations.rds"))
=======
# Updated 09/12/19
# Quang Nguyen
# Script to generate singular and joint ordination plots after ordnation analyses 

library(optparse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(viridis)

option_list <- list(
  make_option("--input", help = "Input file for data loading"),
  make_option("--time", help = "Data time point to extract"),
  make_option("--metab_type", help = "Data type for metabolite data, can be 'tar' or 'untar'"),
  make_option("--tax_type", help = "Data type for taxonomic data, so far only '16S' is supported")
)

opt <- parse_args(OptionParser(option_list = option_list))

output_folder <- paste0("snakemake_output/figures/ordinations/",opt$tax_type, "/", opt$time, "/", opt$metab_type, "/")
data <- readRDS(file = opt$input)

tax_ord <- data$tax_ord
met_ord <- data$met_ord
proc_test <- data$proc_test
mant_test <- data$mant_test 

# Single data ordinations
met_pts <- as.data.frame(met_ord$points)
tax_pts <- as.data.frame(tax_ord$points)
colnames(met_pts) <- colnames(tax_pts) <- c("NMDS1", "NMDS2")
met_plt <- ggscatter(data =met_pts, x = "NMDS1", y = "NMDS2", color = "coral", ellipse = T, ellipse.type = "t", star.plot = T, mean.point = T, 
                    ellipse.alpha = 0, mean.point.size = 3)
met_plt <- ggpar(met_plt, title = "Metabolites")
tax_plt <- ggscatter(data = tax_pts, x = "NMDS1", y = "NMDS2", color = "steelblue", ellipse = T, ellipse.type = "t", star.plot = T, mean.point = T, 
                        ellipse.alpha = 0)
tax_plt <- ggpar(tax_plt, title = "Taxonomy")


saveRDS(tax_plt, file = paste0(output_folder, "tax_ordinations.rds"))
ggsave(filename = paste0(output_folder, "tax_ordinations.svg"), plot = tax_plt, device = "svg")
saveRDS(met_plt, file = paste0(output_folder, "met_ordinations.rds"))
ggsave(filename = paste0(output_folder, "met_ordinations.svg"), plot = met_plt, device = "svg")
    
# joint data set ordination 
plot_dat <- data.frame(rbind(proc_test$Yrot, proc_test$X))
plot_dat[,3] <- c(rep("Metabolite",nrow(proc_test$Yrot)), rep("Taxonomy", nrow(proc_test$X)))
colnames(plot_dat) <- c("NMDS1", "NMDS2", "Ordination")
plt <- ggscatter(plot_dat, x= "NMDS1", y = "NMDS2", 
            color = "Ordination",
            palette = "jco", 
            size = 3, 
            star.plot.lwd = 1,
            ellipse = T,
            ellipse.type = "t",
            star.plot = F, mean.point = F,
            repel = T)
plt <- ggpar(plt,subtitle = paste("Procrustes SoS:", round(proc_test$ss,4), "; Sig:", round(proc_test$signif,4)), palette = viridis(2))

saveRDS(plt, file = paste0(output_folder, "joint_protest_ordinations.rds"))
>>>>>>> 5ececc62da3c5f56fbf3afac6d741fd789227535
ggsave(filename = paste0(output_folder, "joint_protest_ordinations.svg"), plot = plt, device = "svg")