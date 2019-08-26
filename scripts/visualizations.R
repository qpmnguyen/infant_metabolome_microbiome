# This script is for producing ordination plots 
# Quang Nguyen
# Last updated 08/26/19

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
library(ggpubr)
source("./scripts/utils.R")



# ordiation plots 

# function to do testing for significant eigen values  
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

data <- load_main(file = "./data/data_directory.csv")
key <- data$tax.key
key <- key_processing(key, type = "NA")


generate_ordination_plots <- function(time, method, type, data){
  # processing taxonomic data 
  tax <- data[paste0("tax.",time)][[1]]
  tree <- read_tree(treefile = paste0("./data/", time, "_tree.tre"))
  tree <- midpoint.root(tree)
  phylo_obj <- phyloseq(otu_table(tax, taxa_are_rows = F),
                        tax_table(key), 
                        phy_tree(tree))
  phylo_obj <- prune_taxa(taxa_sums(phylo_obj) > 0, phylo_obj)
  phylo_obj <- transform_sample_counts(phylo_obj, function(x) x/sum(x))
  tax <- as(otu_table(phylo_obj), "matrix")
  tax_gunifrac_dist <- MiSPU::GUniFrac(otu.tab = tax, tree = phy_tree(phylo_obj), alpha = 0.5)$GUniF[,,1]
  
  # metabolomics distance
  filename = paste0("./data/",time,"_",method, "NMR_clr_tax.rds")
  met <- readRDS(file = filename)$met
  princomp <- prcomp(met,center = T, scale = T)
  rmt_test <- rmt.test(princomp)
  sig_idx <- which(rmt_test$pvals < 0.05)
  n_PCs <- princomp$x[,sig_idx]
  met_PCs_manhattan_dist <- stats::dist(n_PCs, method = "manhattan")
  
  # Ordinations
  tax_ord <- metaMDS(comm = tax_gunifrac_dist, try = 50, engine = "isoMDS")
  met_ord <- metaMDS(comm = as.matrix(met_PCs_manhattan_dist), try = 50, engine = "isoMDS")
  if(type == "single"){
    met_pts <- as.data.frame(met_ord$points)
    tax_pts <- as.data.frame(tax_ord$points)
    colnames(met_pts) <- colnames(tax_pts) <- c("NMDS1", "NMDS2")
    met_plt <- ggscatter(data =met_pts, x = "NMDS1", y = "NMDS2", color = "coral", ellipse = T, ellipse.type = "t", star.plot = T, mean.point = T, 
                     ellipse.alpha = 0, mean.point.size = 3)
    met_plt <- ggpar(met_plt, title = "Metabolites")
    tax_plt <- ggscatter(data = tax_pts, x = "NMDS1", y = "NMDS2", color = "steelblue", ellipse = T, ellipse.type = "t", star.plot = T, mean.point = T, 
                         ellipse.alpha = 0)
    tax_plt <- ggpar(met_plt, title = "Taxonomy")
    grid <- grid.arrange(met_plt, tax_plt, nrow = 1, ncol = 2)
    ggsave(filename = paste0(time, "_", method, "_single_ordination.png"), plot = plt, width = 15, height = 10, units = "in")
    
  } else if (type == "rot") {
    proc_test <- protest(tax_ord, met_ord)
    plot_dat <- data.frame(rbind(proc_test$Yrot, proc_test$X))
    plot_dat[,3] <- c(rep("Metabolite",nrow(proc_test$Yrot)), rep("Taxonomy", nrow(proc_test$X)))
    colnames(plot_dat) <- c("NMDS1", "NMDS2", "Ordination")
    (plt <- ggscatter(plot_dat, x= "NMDS1", y = "NMDS2", 
              color = "Ordination",
              palette = "jco", 
              size = 1, 
              ellipse = T,
              ellipse.type = "convex",
              star.plot = T, mean.point = T, mean.point.size = 3,
              repel = T))
    plt <- ggpar(plt,subtitle = paste("SoS:", round(proc_test$ss,4), "; Sig:", round(proc_test$signif,4)))
    name <- paste0(time, "_", method, "_protest_ord_plot.png")
    ggsave(filename = name, plot = plt, device = "png", width = 15, height = 10, units = "in")
  }
  # Protest
}


generate_ordination_plots(time = "12M", method = "tar", type = "rot", data = data)
