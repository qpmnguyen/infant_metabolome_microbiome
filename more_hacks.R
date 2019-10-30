library(ggplot2)
library(ggpubr)
library(phyloseq)

data <- readRDS(file = "snakemake_output/analyses/correlation/16S_12M_tar_scca.rds")

data2 <- readRDS(file = "snakemake_output/analyses/correlation/16S_12M_tar_scc.rds")

quantile(data$boot$t, c(0.05, 0.95))
length(which(data[[3]] > data$boot$t0))
range(data[[3]], data$boot$t0)
data$boot$t0

correlation <- data2$cor_mat
cca <- data$cca
tab <- data$tab
tax_idx <- which(cca$u != 0)
met_idx <- which(cca$v != 0)
tab <- as(data$tab, "matrix")
#combined <- cbind(tax[,tax_idx], met[,met_idx])




# Calculate correlation matrix 
#correlation <- correlation[tax_idx, met_idx]
correlation <- cor(tax[,tax_idx], met[,met_idx], method = "spearman")
# grab family-genus names and italicize them 
tax_names <- paste(tab[tax_idx,][,c("Family")], tab[tax_idx,][,c("Genus")]) #Family-Genus names
tax_names <- as.expression(sapply(tax_names, function(x){
  bquote(italic(.(x)) ~ "spp.")
}))

# plotting 
row <- data.frame("sCCA Loading" = as.factor(ifelse(cca$u[tax_idx] > 0, "+","-")), check.names = F)
rownames(row) <- rownames(correlation)
col <- data.frame("sCCA Loading" = as.factor(ifelse(cca$v[met_idx] > 0, "+","-")), check.names = F)
rownames(col) <- colnames(correlation)
ann_colors <- list("sCCA Loading" = c("#E7B800", "#00AFBB"))
names(ann_colors$"sCCA Loading") <- as.factor(c("+", "-")) 
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
  annotation_names_row = F,
  annotation_names_col = F,
  drop_levels       = TRUE,
  fontsize          = 11,
  cluster_cols = F,
  cluster_rows = F
)
cca_heatmap <- as.ggplot(cca_heatmap)
cca_heatmap
