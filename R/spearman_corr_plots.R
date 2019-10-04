library(pheatmap)
library(viridis)
library(optparse)
library(ggplotify)
library(ggplot2)
option_list <- list(
  make_option("--input", help = "Input file for data loading"),
  make_option("--time", help = "Data time point to extract"),
  make_option("--metab_type", help = "Data type for metabolite data, can be 'tar' or 'untar'"),
  make_option("--tax_type", help = "Data type for taxonomic data, so far only '16S' is supported")
)

opt <- parse_args(OptionParser(option_list = option_list))

data <- readRDS(file = opt$input)

cor_mat <- data$cor_mat
adj_mat <- data$p_mat
tab <- as(data$tax_tab, "matrix")

# making sure that significant taxa are highlighted 
sig <- t(apply(adj_mat, 1, function(x){
  ifelse(x < 0.05, 1,0)
}))

idx_tax <- which(apply(sig,1, function(x){
  all(x == 0)
}) == F)
idx_met <- which(apply(sig,2, function(x){
  all(x == 0)
}) == F)

sig <- ifelse(sig == 1, "*", "")

tax_names <- paste(tab[,c("Family")], tab[,c("Genus")]) #Family-Genus names
tax_names <- as.expression(sapply(tax_names, function(x){
  bquote(italic(.(x)) ~ "spp.")
}))

met_names <- colnames(cor_mat)
if (opt$metab_type == "tar"){
  met_names[36] <- "pi-Methylhistidine"
  display = T
} else {
  display = F
}


corr_heatmap <- pheatmap(
  mat = cor_mat,
  color = viridis(40),
  border_color = NA,
  show_colnames = display,
  show_rownames = T,
  display_numbers = sig,
  labels_row = tax_names,
  labels_col = met_names,
  annotation_names_row = F,
  annotation_names_col = F,
  drop_levels = TRUE,
  fontsize = 10
)

corr_heatmap <- as.ggplot(corr_heatmap)
output_name = paste0("snakemake_output/figures/correlation/", opt$tax_type, "_", opt$time, "_", opt$metab_type, "_scc_plots")
saveRDS(corr_heatmap, file = paste0(output_name, ".rds"))
ggsave(plot = corr_heatmap, filename = paste0(output_name, ".svg"), device = "svg", width = 10, height = 10)
