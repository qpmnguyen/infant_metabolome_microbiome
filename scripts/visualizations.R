library(ggplot2)
library(phyloseq)
library(reshape2)
library(compositions)
library(MCTools)
library(gridExtra)
library(phytools)
library(MiSPU)
library(ggfortify)
library(vegan)
library(ade4)
source("./scripts/utils.R")

# relative abundance plots 
tar <- readRDS(file = './data/tarNMR_clr_tax.rds')
untar <- readRDS(file = "./data/untarNMR_clr_tax.rds")
data <- load_main(file = "./data/data_directory.csv")
tax <- data$tax.12M
key <- data$tax.key
key <- key_processing(key, type = "NA")
phylo_obj <- phyloseq(otu_table(tax, taxa_are_rows = F),
                      tax_table(key))

phylo_obj <- tax_glom(phylo_obj, taxrank = "Family")
phylo_obj <- prune_taxa(taxa_sums(phylo_obj) > 0, phylo_obj)
phylo_obj <- transform_sample_counts(phylo_obj, function(x) x/sum(x))

p1 <- plot_bar(phylo_obj, fill = "Family", y = "Abundance", title = "Microbial Family") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
ggsave(plot = p1,device = "pdf", filename = "./docs/12M_family_level_RA.pdf")


met <- unclass(acomp(met))
p2 <- ggplot(melt(met), aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity", color = "black") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") + 
  labs(x = "Samples", y = "Abundance", title = "Concentration-fitted Metabolites")
ggsave(plot = p2, device = "pdf", filename="./docs/12M_tarNMR_RA.pdf")

combine_plot <- grid.arrange(p1, p2, nrow = 2, ncol = 1)
ggsave(combine_plot, device = "pdf", filename = "./docs/12M_family_metabolite_RA.pdf")

# ordiation plots 
tax <- data$tax.12M
key <- data$tax.key
key <- key_processing(key, type = "NA")
tree <- read_tree(treefile = "./data/12M_tree.tre")
tree <- midpoint.root(tree)
met <- data$metabo.12M
phylo_obj <- phyloseq(otu_table(tax, taxa_are_rows = F),
                      tax_table(key), 
                      phy_tree(tree))

phylo_obj <- tax_glom(phylo_obj, taxrank = "Genus")
phylo_obj <- prune_taxa(taxa_sums(phylo_obj) > 0, phylo_obj)
otu_table(phylo_obj) <- otu_table(phylo_obj) + 1 # pseudocount
phylo_obj <- transform_sample_counts(phylo_obj, function(x) x/sum(x))

aich_dist <- coda.base::dist(otu_table(phylo_obj), method = "aitchison")

#temp <- unclass(clr(otu_table(phylo_obj)))
#otu_table(phylo_obj) <- otu_table(temp, taxa_are_rows = F)

gunifrac_dist <- MiSPU::GUniFrac(otu.tab = otu_table(phylo_obj), tree = phy_tree(phylo_obj), alpha = 0.5)$GUniF[,,1]
bray_dist <- vegdist(otu_table(phylo_obj), method = "bray")
manhattan_dist <- stats::dist(met, method = "manhattan")

plot(as.matrix(bray_dist), as.matrix(manhattan_dist))
abline(h = quantile(as.matrix(manhattan_dist), 0.25), col = "red")
abline(h = quantile(as.matrix(manhattan_dist), 0.75), col = "red")

length(gunifrac_dist[lower.tri(gunifrac_dist)])





