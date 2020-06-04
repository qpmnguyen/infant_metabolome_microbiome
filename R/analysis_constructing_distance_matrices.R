# Updated 09/13/2019
# Quang Nguyen
# This script is to generate rds files of distance matrices 

library(vegan)
library(MiSPU)
library(phyloseq)
library(RMTstat)
library(phytools)
library(optparse)
library(plyr)
library(compositions)

option_list <- list(
  make_option("--input", help = "Input file for data loading - use raw data instead of processed"),
  make_option("--output", help = "Output file for data loading")
)

opt <- parse_args(OptionParser(option_list = option_list))

# supplemental functions 
#' @title Testing for significance of principal components
#' @description Testing for significance of eigenvalues using Tracy-Widom test statistics. Details
#'    see Frost et al. 2015
#' @param prcomp.output Principal component analysis object.  
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


##### main implementation -------------------------------------------------------------------
data <- readRDS(file = opt$input)
TIME <- strsplit(opt$input, "_")[[1]][2]
METAB <- strsplit(opt$input, "_")[[1]][3]
tree <- phytools::midpoint.root(phy_tree(data)) # midpoint rooting the tree 
phyloseq::phy_tree(data) <- tree

#### gunifrac requires that data be kept at ASV level
get_gunifrac <- function(data){
  data <- transform_sample_counts(data, function(x) x/sum(x))
  # data <- filter_taxa(data, function(x) mean(x) > 0.005e-2, TRUE)
  tax <- as(otu_table(data), "matrix")
  tax_gunifrac <- as.dist(MiSPU::GUniFrac(otu.tab = tax, tree = phy_tree(data))$d5)
  return(tax_gunifrac)
}

tax_gunifrac <- get_gunifrac(data)
# TODO check for using euclidean distance vs manhattan distance 

# getting other data sets  
# processing taxonomic data 
data <- tax_glom(data, taxrank = "Genus") #aggregate to genus level 
data <- transform_sample_counts(data, function(x) x/ sum(x)) # convert to relative abundance
data <- filter_taxa(data, function(x) mean(x) > 0.005e-2, TRUE) # mean of at least 0.005 percent (ref. Bokulich 2013)
tax <- as(otu_table(data), "matrix")

# constructing distance matrices 
tax_euclidean <- dist(unclass(compositions::clr(tax)), method = "euclidean") # euclidean distance under clr transformation 

# processing the metabolomics data 
if (METAB == "tar"){
  met <- as(sample_data(data), "matrix") # this procedure turns everything into character
  met <- log(met + 1) 
} else if (METAB == "untar"){
  met <- as(sample_data(data), "matrix")
  met <- unclass(acomp(met))
  met <- asin(sqrt(met))
} else {
  stop("Not supported met type")
}

met_euclidean <- stats::dist(met, "euclidean")
met_manhattan <- stats::dist(met, "manhattan")
# PCA-based distance  
# princomp <- prcomp(met,center = T, scale = T) # principal component analyses 
# rmt_test <- rmt.test(princomp) # test for significant PCs
# sig_idx <- which(rmt_test$pvals < 0.05) 
# n_PCs <- princomp$x[,sig_idx] # select significant eigen vectors 
# met_PC_dist_euclidean <- stats::dist(n_PCs, method = "euclidean") # construct manhattan distance matrix

# saving output  
result <- list(tax_gunifrac = tax_gunifrac, 
               tax_euclidean = tax_euclidean,
               met_euclidean = met_euclidean,
               met_manhattan = met_manhattan)
saveRDS(result, file = opt$output)