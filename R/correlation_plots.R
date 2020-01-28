library(ggpubr)
library(optparse)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(ggplotify)
library(grid)
library(phyloseq)
library(reshape2)
library(ggdendro)
# important colors to keep consistency  
# #00AFBB nice blue
# #E7B800 nice yellow
# #FC4E07 nice red 


# option_list <- list(
#  make_option("--time", help = "Data time point to extract"),
#  make_option("--metab_type", help = "Data type for metabolite data, can be 'tar' or 'untar'")
# )

# opt <- parse_args(OptionParser(option_list = option_list))

# data <- readRDS(file = opt$input)
#data <- readRDS(file = "snakemake_output/analyses/correlation/16S_6W_tar_scca.rds")

correlation <- readRDS(file = opt$correlation)$cor_mat
#correlation <- readRDS(file = "snakemake_output/analyses/correlation/16S_6W_tar_scc.rds")$cor_mat
# Plotting CCA results 
cca <- data$cca
tax_idx <- which(cca$u != 0)
met_idx <- which(cca$v != 0)
tab <- as(data$tab, "matrix")
#combined <- cbind(tax[,tax_idx], met[,met_idx])
tax_idx
met_idx

dim(correlation)
# Calculate correlation matrix 
correlation <- correlation[tax_idx, met_idx]
#correlation <- cor(tax[,tax_idx], met[,met_idx], method = "spearman")
# grab family-genus names and italicize them 

# version which has family and genus name italicized 
tax_names <- paste(tab[tax_idx,][,c("Genus")]) #Family-Genus names
tax_names <- as.expression(sapply(tax_names, function(x){
  bquote(italic(.(x)) ~ "spp.")
}))

#tax_names <- paste(tab[tax_idx,][,c("Genus")], "spp.")
family_names <- tab[tax_idx,][,c("Family")]

# plotting 
#row <- data.frame("sCCA Loading" = as.factor(ifelse(cca$u[tax_idx] > 0, "+","-")), "Family" = as.factor(family_names), check.names = F)
row <- data.frame("sCCA Loading" = as.factor(ifelse(cca$u[tax_idx] > 0, "+","-")), check.names = F)

rownames(row) <- rownames(correlation)
col <- data.frame("sCCA Loading" = as.factor(ifelse(cca$v[met_idx] > 0, "+","-")), check.names = F)
rownames(col) <- colnames(correlation)
ann_colors <- list("sCCA Loading" = c(cividis(10)[1], cividis(10)[10]))
names(ann_colors$"sCCA Loading") <- as.factor(c("+", "-")) 
raw_pval <- length(which(data[[3]] >= data$boot$t0))/length(data[[3]])

if (raw_pval == 0){
  title = paste("SCCA Correlation:", round(data$boot$t0,3), 
               "(Bootstrapped 95% CI:", round(quantile(data$boot$t, 0.05),3), "-", round(quantile(data$boot$t, 0.95), 3),
               "; Permutation p-value < 0.001)")
} else {
  title = paste("SCCA Correlation:", round(data$boot$t0,3),  
               "(Bootstrapped 95% CI:", round(quantile(data$boot$t, 0.05),3), "-", round(quantile(data$boot$t, 0.95), 3),
               "; Permutation p-value :", raw_pval,")")
}

#ord <- hclust(dist(correlation, method = "euclidean"))
#cor_melt <- data.frame(melt(correlation), stringsAsFactors = F)
#cor_melt$Var1 <- factor(cor_melt$Var1, levels = rownames(correlation)[ord$order])
#dendro <- as.dendrogram(ord)
#dendro_plt <- ggdendrogram(data = dendro, rotate = T)


#ggplot(cor_melt, aes(y = Var1, x = Var2, fill = value)) + geom_tile() + scale_fill_viridis() + theme_cleveland()

# removed main = title
cca_heatmap <- pheatmap(
  mat               = correlation,
  annotation_row    = row,
  annotation_colors = ann_colors,
  annotation_col    = col,
  color             = viridis(40),
  annotation_names_row = T,
  annotation_names_col = T,
  border_color      = NA,
  show_colnames     = T,
  show_rownames     = T,
  labels_row        = tax_names,
  drop_levels       = TRUE,
  fontsize          = 11,
  cex = 1, 
  legend = T
)

cca_heatmap <- as.ggplot(cca_heatmap)
output_name = paste0("snakemake_output/figures/correlation/", opt$tax_type, "_", opt$time, "_", opt$metab_type, "_scca_plots")
saveRDS(cca_heatmap, file = paste0(output_name, ".rds"))
ggsave(plot = cca_heatmap, filename = paste0(output_name, ".png"), device = "png", width = 13, height = 16)
ggsave(plot = cca_heatmap, filename = paste0(output_name, ".svg"), device = "svg", width = 8, height = 8)




