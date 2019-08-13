# This script is for producing basic visualizations such as relative abundances and ordination 
# Quang Nguyen
# Last updated 08/09/19

library(dplyr)
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
library(RMTstat)
source("./scripts/utils.R")

# relative abundance plots 
#tar <- readRDS(file = './data/tarNMR_clr_tax.rds')
#untar <- readRDS(file = "./data/untarNMR_clr_tax.rds")
data <- load_main(file = "./data/data_directory.csv")
tax <- data$tax.6W
key <- data$tax.key
key <- key_processing(key, type = "NA")
phylo_obj <- phyloseq(otu_table(tax, taxa_are_rows = F),
                      tax_table(key))

phylo_obj <- tax_glom(phylo_obj, taxrank = "Family")
phylo_obj <- prune_taxa(taxa_sums(phylo_obj) > 0, phylo_obj)
phylo_obj <- transform_sample_counts(phylo_obj, function(x) x/sum(x))

p1 <- plot_bar(phylo_obj, fill = "Family", y = "Abundance", title = "Microbial Family") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
ggsave(plot = p1,device = "pdf", filename = "./docs/6W_family_level_RA.pdf")
ggsave(plot = p1,device = "png", filename = "./docs/6W_family_level_RA.png",dpi = 'retina')

met <- data$metabo.6W
met <- unclass(acomp(met))
p2 <- ggplot(melt(met), aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity", color = "black") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") + 
  labs(x = "Samples", y = "Abundance", title = "Concentration-fitted Metabolites")
ggsave(plot = p2, device = "pdf", filename="./docs/6W_tarNMR_RA.pdf")
ggsave(plot = p2, device = "png", filename="./docs/6W_tarNMR_RA.png", dpi = "retina")

combine_plot <- grid.arrange(p1, p2, nrow = 2, ncol = 1)
ggsave(combine_plot, device = "pdf", filename = "./docs/6W_family_metabolite_RA.pdf")
ggsave(combine_plot, device = "png", filename = "./docs/6W_family_metabolite_RA.png", dpi = "retina")
# ordiation plots 

rmt.test <- function(prcomp.output){
  eigenvals <- (prcomp.output$sdev)^2
  n.eigenvec <- nrow(prcomp.output$rotation)
  sample.size <- nrow(prcomp.output$x)
  tws <- c()
  pvals <- c()
  for (i in 1:length(eigenvals)){
    eval <- eigenvals[i] #selecting the eigenvalue
    pdim <- n.eigenvec - i + 1 #number of dimensions of the wishart matrix
    par <- WishartMaxPar(ndf = sample.size, pdim = pdim)
    tws[i] <- (eval - par$centering)/par$scaling
    pvals[i] <- pgamma(tws[i] + 9.84801, shape=46.446, scale=0.186054, lower.tail = F)
  }
  result = list(tws = tws, pvals = pvals)
  return(result)
}

tax <- data$tax.6W
key <- data$tax.key
key <- key_processing(key, type = "NA")
tree <- read_tree(treefile = "./data/6W_tree.tre")
tree <- midpoint.root(tree)
met <- data$metabo.6W
met_untar <-readRDS(file = "./data/6W_untarNMR_clr_tax.rds")$met

phylo_obj <- phyloseq(otu_table(tax, taxa_are_rows = F),
                      tax_table(key), 
                      phy_tree(tree))

#phylo_obj <- tax_glom(phylo_obj, taxrank = "Genus")
phylo_obj <- prune_taxa(taxa_sums(phylo_obj) > 0, phylo_obj)
#otu_table(phylo_obj) <- otu_table(phylo_obj) + 1 # pseudocount
phylo_obj <- transform_sample_counts(phylo_obj, function(x) x/sum(x))

tax <- as(otu_table(phylo_obj), "matrix")

tax_aich_dist <- coda.base::dist(tax, method = "aitchison")
tax_gunifrac_dist <- MiSPU::GUniFrac(otu.tab = tax, tree = phy_tree(phylo_obj), alpha = 0.5)$GUniF[,,1]
tax_bray_dist <- vegdist(tax, method = "bray")

clr_tax <- unclass(clr(as(otu_table(phylo_obj), "matrix")))

tax_euclid_dist <- stats::dist(clr_tax, method = "euclidean")
tax_manhattan_dist <- stats::dist(clr_tax, method = "manhattan")


princomp <- prcomp(met,center = T, scale = T)
rmt_test <- rmt.test(princomp)
#sig_idx <- which(princomp$sdev^2 >= 1)
sig_idx <- which(rmt_test$pvals < 0.05)
n_PCs <- princomp$x[,sig_idx]

met_PCs_manhattan_dist <- stats::dist(n_PCs, method = "manhattan")
met_PCs_euclidean_dist <- stats::dist(n_PCs, method = "euclidean")
met_manhattan_dist <- stats::dist(met, method = "manhattan")
met_euclid_dist <- stats::dist(met, method = "euclidean")
met_cosine_dist <- stats::dist(met, method = "cosine")
untarmet_aich_dist <- coda.base::dist(met_untar, method = "aitchison")
untarmet_manhattahn <- dist(asin(sqrt(met_untar)),  method = "manhattan")

generate_pair_plot <- function(dist_tax, dist_met, rank, xlab, ylab){
  dist_tax <- as.matrix(dist_tax)
  dist_met <- as.matrix(dist_met)
  dist_tax[lower.tri(dist_tax)] <- 0
  dist_met[lower.tri(dist_met)] <- 0
  dist_tax <- melt(dist_tax) %>% filter(value > 0)
  dist_met <- melt(dist_met) %>% filter(value > 0)
  if (rank == T){
    dist_tax$value <- rank(dist_tax$value)
    xlab <- paste("Ranked", xlab, "distance")
  } else {
    xlab <- paste(xlab, "Distance")
  }
  ylab <- paste(ylab, "Distance")
  plot_dat <- data.frame(cbind(dist_tax$value, dist_met$value))
  colnames(plot_dat) <- c("tax", "met")
  plt <- ggplot(plot_dat, aes(x = tax, y = met)) + 
    geom_point() + 
    stat_binhex(bins = 80) + 
    geom_hline(yintercept = quantile(dist_met$value, 0.5), colour ="red") + 
    labs(x = xlab, y = ylab) + theme_bw()
    
  return(plt)
}

(p <- generate_pair_plot(tax_gunifrac_dist, met_PCs_manhattan_dist, rank = F, xlab = "Taxonomic", ylab = "Metabolite") + labs(title = "Targeted"))
(p2 <- generate_pair_plot(tax_gunifrac_dist, untarmet_manhattahn, rank = F, xlab = "Taxonomic", ylab = "Metabolite") + labs(title = "Untargeted"))

(joined <- grid.arrange(grobs = list(p = p, p2 = p2), nrow = 1, col = 2))
ggsave(plot = joined, filename = "./docs/distance_unranked_comparison_12M.png", device = "png", width = 8, height = 5, units = "in")

# Comparing ordinations  
tax_ord <- metaMDS(comm = tax_gunifrac_dist, try = 50, engine = "isoMDS")
met_ord <- metaMDS(comm = as.matrix(met_PCs_manhattan_dist), try = 50, engine = "isoMDS")
met_ord_2 <- metaMDS(comm = as.matrix(untarmet_manhattahn), try = 50, engine = "isoMDS")

proc_test <- protest(tax_ord, met_ord)
proc_test <- protest(tax_ord, met_ord_2)

(mantel_test <- mantel(tax_gunifrac_dist, untarmet_manhattahn))
(mantel_test <- mantel(tax_gunifrac_dist, met_PCs_manhattan_dist))
plot_dat <- data.frame(rbind(proc_test$Yrot, proc_test$X))
plot_dat[,3] <- c(rep("Metabolite",nrow(proc_test$Yrot)), rep("Taxonomy", nrow(proc_test$X)))
(plt <- ggplot(plot_dat, aes(x = NMDS1, y = NMDS2, colour = V3)) + geom_point(size = 2) + scale_color_discrete(name = "Ordination Type") +
    labs(title = "Targeted") + 
  annotate("text", x = 0.6, y = -0.04, label = paste("Procrustes SS:", round(proc_test$ss,4)), color = "red") +
  annotate("text", x = 0.6, y = -0.05, label = paste("Sig:", proc_test$signif), color = "red") + theme_bw())

plot <- grid.arrange(plt, plt2, nrow = 1, ncol = 2)


ggsave("./docs/6W_procrustes_tar_untar.png", device = "png", plot = plot, units = "in", width = 15, height = 5)
