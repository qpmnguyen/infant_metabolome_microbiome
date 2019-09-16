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

option_list <- list(
  make_option("--input", help = "Input file for data loading"),
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

# processing the generated data 
tax <- as(otu_table(data), "matrix")
tree <- phytools::midpoint.root(phy_tree(data)) # midpoint rooting the tree 
tax_dist <- MiSPU::GUniFrac(otu.tab = tax, tree = tree, alpha = 0.5)$GUniF[,,1]

# processing the metabolomics data 
met <- as(sample_data(data), "matrix") # this procedure turns everything into character
temp_names <- rownames(met)
met <- apply(met, 2, as.numeric)
rownames(met) <- temp_names
rm(temp_names)

control_idx <- grep("DSS", colnames(met))
met <- met[,-c(1,control_idx)] # removing some batch information 
princomp <- prcomp(met,center = T, scale = T) # principal component analyses 
rmt_test <- rmt.test(princomp) # test for significant PCs
sig_idx <- which(rmt_test$pvals < 0.05) 
n_PCs <- princomp$x[,sig_idx] # select significant eigen vectors 
met_dist <- stats::dist(n_PCs, method = "manhattan") # construct manhattan distance matrix

# saving output  
result <- list(tax_dist = tax_dist, met_dist = met_dist)
saveRDS(result, file = opt$output)