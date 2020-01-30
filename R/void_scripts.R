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
