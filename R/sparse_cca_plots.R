library(ggpubr)
library(optparse)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(ggplotify)
library(grid)
library(phyloseq)
# important colors to keep consistency  
# #00AFBB nice blue
# #E7B800 nice yellow
# #FC4E07 nice red 


option_list <- list(
  make_option("--input", help = "Input file for data loading"),
  make_option("--correlation", help = "Pre-generated correlation matrix from spearman correlation step"),
  make_option("--time", help = "Data time point to extract"),
  make_option("--metab_type", help = "Data type for metabolite data, can be 'tar' or 'untar'"),
  make_option("--tax_type", help = "Data type for taxonomic data, so far only '16S' is supported")
)

opt <- parse_args(OptionParser(option_list = option_list))

data <- readRDS(file = opt$input)
#data <- readRDS(file = "snakemake_output/analyses/correlation/16S_12M_tar_scca.rds")

correlation <- readRDS(file = opt$correlation)$cor_mat
#correlation <- readRDS(file = "snakemake_output/analyses/correlation/16S_12M_tar_scc.rds")$cor_mat
# Plotting CCA results 
cca <- data$cca
tax_idx <- which(cca$u != 0)
met_idx <- which(cca$v != 0)
tab <- as(data$tab, "matrix")
#combined <- cbind(tax[,tax_idx], met[,met_idx])

# Calculate correlation matrix 
correlation <- correlation[tax_idx, met_idx]
#correlation <- cor(tax[,tax_idx], met[,met_idx], method = "spearman")
# grab family-genus names and italicize them 

# version which has family and genus name italicized 
#tax_names <- paste(tab[tax_idx,][,c("Family")], tab[tax_idx,][,c("Genus")]) #Family-Genus names
#tax_names <- as.expression(sapply(tax_names, function(x){
#  bquote(italic(.(x)) ~ "spp.")
#}))

tax_names <- paste(tab[tax_idx,][,c("Genus")], "spp.")
family_names <- tab[tax_idx,][,c("Family")]

# plotting 
row <- data.frame("sCCA Loading" = as.factor(ifelse(cca$u[tax_idx] > 0, "+","-")), "Family" = as.factor(family_names), check.names = F)
rownames(row) <- rownames(correlation)
col <- data.frame("sCCA Loading" = as.factor(ifelse(cca$v[met_idx] > 0, "+","-")), check.names = F)
rownames(col) <- colnames(correlation)
ann_colors <- list("sCCA Loading" = c("#E7B800", "#00AFBB"))
names(ann_colors$"sCCA Loading") <- as.factor(c("+", "-")) 
raw_pval <- length(which(data[[3]] >= data$boot$t0))/length(data[[3]])
if (raw_pval == 0){
  title = paste("Sparse CCA Correlation:", round(data$boot$t0,3), 
               "(Bootstrapped 95% CI:", round(quantile(data$boot$t, 0.05),3), "-", round(quantile(data$boot$t, 0.95), 3),
               "; Permutation p-value < 0.0001)")
} else {
  title = paste("Sparse CCA Correlation:", round(data$boot$t0,3), 
               "(Bootstrapped 95% CI:", round(quantile(data$boot$t, 0.05),3), "-", round(quantile(data$boot$t, 0.95), 3),
               "; Permutation p-value :", raw_pval,")")
}
cca_heatmap <- pheatmap(
  mat               = correlation,
  color             = viridis(40),
  border_color      = NA,
  show_colnames     = T,
  show_rownames     = T,
  labels_row        = tax_names,
  annotation_row    = row,
  annotation_col    = col,
  annotation_colors = ann_colors,
  annotation_names_row = T,
  annotation_names_col = T,
  drop_levels       = TRUE,
  fontsize          = 11,
  main = title
)
cca_heatmap <- as.ggplot(cca_heatmap)
output_name = paste0("snakemake_output/figures/correlation/", opt$tax_type, "_", opt$time, "_", opt$metab_type, "_scca_plots")
saveRDS(cca_heatmap, file = paste0(output_name, ".rds"))
ggsave(plot = cca_heatmap, filename = paste0(output_name, ".png"), device = "png", width = 13, height = 16)





