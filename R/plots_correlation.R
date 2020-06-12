# This script seeks to perform plotting of univariate and multivariate correlation analyses 
library(ggpubr)
library(ggplot2)
library(cowplot)
library(glue)
library(phyloseq)
library(viridis)
library(ggplotify)
library(pheatmap)
library(optparse)


option_list <- list(
  make_option("--scca", help = "Scca input file"),
  make_option("--correlation", help = "correlation input file")
)

opt <- parse_args(OptionParser(option_list = option_list))
time <- tail(strsplit(strsplit(opt$scca, "_")[[1]][1],"/")[[1]],1)
metab <- strsplit(opt$scca,"_")[[1]][2]


if (metab == "untar"){
  show_col <- F
  wid <- 14
  hei <- 10
} else {
  show_col <- T
  wid <- 10
  hei <- 10
}
### loading some stuff #### 
scca <- readRDS(file = opt$scca)
corr <- readRDS(file = opt$correlation)
names(scca) <- c("cca", "boot", "perm", "tab")
cor_mat <- corr$cor_mat
adj_mat <- corr$p_mat

####---------plotting original heatmap---------------------#####
tab <- as(corr$tax_tab, "matrix")
sig <- t(apply(adj_mat, 1, function(x){
  ifelse(x < 0.05, 1,0)
}))

idx_tax <- which(apply(sig,1, function(x){
  all(x == 0)
}) == F)
idx_met <- which(apply(sig,2, function(x){
  all(x == 0)
}) == F)

#sig <- ifelse(sig == 1, "X", "")

tax_names <- paste(tab[,c("Genus")]) #Family-Genus names
tax_names <- as.expression(sapply(tax_names, function(x){
  bquote(italic(.(x)) ~ "spp.")
}))
family_names <- tab[,c("Family")]

met_names <- colnames(cor_mat)

row <- data.frame("Family" = as.factor(family_names), check.names = F)
rownames(row) <- rownames(cor_mat)

cor_mat[which(sig == 0, arr.ind = T)] <- NA 

corr_heatmap <- pheatmap(
  mat = cor_mat,
  color = viridis(40),
  border_color = NA,
  show_colnames = show_col,
  show_rownames = T,
  #display_numbers = sig,
  number_color = "red",
  labels_row = tax_names,
  labels_col = met_names,
  annotation_names_row = F,
  annotation_names_col = F,
  drop_levels = TRUE,
  fontsize = 10, 
  legend = F,cluster_rows = F, cluster_cols = F,
  na_col = "beige"
)

corr_heatmap <- as.ggplot(corr_heatmap)


#####-------plotting cca---------------####
cca <- scca$cca
sig <- t(apply(adj_mat, 1, function(x){
  ifelse(x < 0.05, 1,0)
}))
tax_idx <- which(cca$u != 0)
met_idx <- which(cca$v != 0)
cca_mat <- cor_mat[tax_idx, met_idx]
sig_mat <- sig[tax_idx, met_idx]
#sig_mat <- ifelse(sig_mat == 1, "X", "")
# version which has family and genus name italicized 
tax_names <- paste(tab[tax_idx,][,c("Genus")]) #Family-Genus names
tax_names <- as.expression(sapply(tax_names, function(x){
  bquote(italic(.(x)) ~ "spp.")
}))


#tax_names <- paste(tab[tax_idx,][,c("Genus")], "spp.") # italicized species names 
family_names <- tab[tax_idx,][,c("Family")]

# the following implements family name as well 
# row <- data.frame("sCCA Loading" = as.factor(ifelse(cca$u[tax_idx] > 0, "+","-")), "Family" = as.factor(family_names), check.names = F)
row <- data.frame("sCCA Loading" = as.factor(ifelse(cca$u[tax_idx] > 0, "+","-")), check.names = F)
rownames(row) <- rownames(cca_mat)
col <- data.frame("sCCA Loading" = as.factor(ifelse(cca$v[met_idx] > 0, "+","-")), check.names = F)
rownames(col) <- colnames(cca_mat)
ann_colors <- list("sCCA Loading" = c(cividis(10)[1], cividis(10)[10]))
names(ann_colors$"sCCA Loading") <- as.factor(c("+", "-")) 
raw_pval <- length(which(scca$perm >= scca$boot$t0))/length(scca$perm)

if (raw_pval == 0){
  title = paste("SCCA Correlation:", round(scca$boot$t0,3), 
                "(Bootstrapped 95% CI:", round(quantile(scca$boot$t, 0.05),3), "-", round(quantile(scca$boot$t, 0.95), 3),
                "; Permutation p-value < 0.001)")
} else {
  title = paste("SCCA Correlation:", round(scca$boot$t0,3),  
                "(Bootstrapped 95% CI:", round(quantile(scca$boot$t, 0.05),3), "-", round(quantile(scca$boot$t, 0.95), 3),
                "; Permutation p-value :", raw_pval,")")
}
cca_mat[which(sig_mat == 0, arr.ind = 1)] <- NA
reduced_heatmap <-  pheatmap(
  mat               = cca_mat,
  annotation_row    = row,
  annotation_colors = ann_colors,
  annotation_col    = col,
  color             = viridis(40),
  annotation_names_row = T,
  annotation_names_col = T,
  border_color      = NA,
  #display_numbers   = sig_mat,
  cluster_rows = F,
  cluster_cols = F,
  number_color      = "red",
  show_colnames     = show_col,
  show_rownames     = T,
  labels_row        = tax_names,
  drop_levels       = TRUE,
  fontsize          = 11,
  cex = 1, 
  legend = T,
  na_col = "beige"
)
reduced_heatmap <- as.ggplot(reduced_heatmap)

for (i in c("svg", "png")){
  ggsave(plot = corr_heatmap, filename = glue("output/figures/correlation/{time}_{metab}_correlation.{ext}", 
                                       time = time, metab = metab, ext = i), device = i, dpi = 300, 
                                        width = wid, height = hei,  units = "in")
  ggsave(plot = reduced_heatmap, filename = glue("output/figures/correlation/{time}_{metab}_scca.{ext}",
                                          time = time, metab = metab, ext = i), device = i, dpi = 300, 
                                        width = wid, height = hei, units = "in")
}
saveRDS(corr_heatmap, file = glue("output/figures/correlation/{time}_{metab}_correlation.rds", 
                                  time = time, metab = metab))
saveRDS(reduced_heatmap, file = glue("output/figures/correlation/{time}_{metab}_scca.rds",
                                     time = time, metab = metab))

cat(title, file = glue("output/figures/correlation/{time}_{metab}_labels.txt", time = time, metab = metab))

print(glue("Time is {time} and metab type is {metab}", time = time, metab = metab))




