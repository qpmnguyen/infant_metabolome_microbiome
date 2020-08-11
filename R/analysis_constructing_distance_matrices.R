# Updated 08/11/2020
# Quang Nguyen
# This script is to generate rds files of distance matrices 
# Since gunifrac requires raw data, this script should only take in raw data only 

library(vegan) # Other ecological distances
library(MiSPU) # fast implementation of GUnifrac
library(phyloseq) # Working with phyloseq objects
library(phytools) # mid point root
library(plyr) # some phyloseq transformations require plyr
library(compositions) # for compositional transformations
library(rlang)


source("R/utils.R")

# Getting gunifrac from data 
get_tax_dist <- function(data, dist){
  if (!dist %in% c("gunifrac", "euclidean", "manhattan")){
    stop("Unsupported distance type")
  } else if (dist == "gunifrac"){
    tree <- phytools::midpoint.root(phy_tree(data))
    phy_tree(data) <- tree
    data <- transform_sample_counts(data, function(x) x/sum(x))
    # data <- filter_taxa(data, function(x) mean(x) > 0.005e-2, TRUE)
    tax <- as(otu_table(data), "matrix")
    dist_mat <- as.dist(MiSPU::GUniFrac(otu.tab = tax, tree = phy_tree(data))$d5)
  } else {
    data <- tax_glom(data, taxrank = "Genus") #aggregate to genus level 
    data <- transform_sample_counts(data, function(x) x/ sum(x)) # convert to relative abundance
    data <- filter_taxa(data, function(x) mean(x) > 0.005e-2, TRUE) # mean of at least 0.005 percent (ref. Bokulich 2013)
    tax <- as(otu_table(data), "matrix")
    dist_mat <- dist(unclass(compositions::clr(tax)), method = dist)
  }
  return(dist_mat)
}

get_met_dist <- function(data, mettype, dist){
  if (mettype == "tar"){
    met <- as(sample_data(data), "matrix") # this procedure turns everything into character
    met <- log(met + 1) 
  } else if (mettype == "untar"){
    met <- as(sample_data(data), "matrix")
    met <- unclass(acomp(met))
    met <- asin(sqrt(met))
  }
  if (!dist %in% c("manhattan", "euclidean")){
    stop("Not supported distance metrics")
  } else {
    dist <- stats::dist(met, "manhattan")
  }
}

dist_main <- function(data, mettype, dist_met, dist_tax){
  tax_dist <- get_tax_dist(data, dist_tax)
  met_dist <- get_met_dist(data, dist_met, mettype = mettype)
  output <- list(tax = tax_dist, met = met_dist)
  names_output <- c(glue("tax_{dist_tax}", dist_tax = dist_tax), 
                    glue("met_{dist_met}", dist_met = dist_met))
  names(output) <- names(output)
  return(output)
}

# constructing distance matrices 
 # euclidean distance under clr transformation 

# processing the metabolomics data 

# PCA-based distance  
# princomp <- prcomp(met,center = T, scale = T) # principal component analyses 
# rmt_test <- rmt.test(princomp) # test for significant PCs
# sig_idx <- which(rmt_test$pvals < 0.05) 
# n_PCs <- princomp$x[,sig_idx] # select significant eigen vectors 
# met_PC_dist_euclidean <- stats::dist(n_PCs, method = "euclidean") # construct manhattan distance matrix

