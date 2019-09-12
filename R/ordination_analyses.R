# Updated 09/11/2019
# Quang Nguyen
# This script is to generate rds file of ordnation results for plotting 

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
data <- filter_taxa(data, function(x) count(x > 0)[2,2] >= 0.1*length(x), TRUE) # asvs as to be in at least 10% of samples 
data <- tax_glom(data, taxrank = "Genus") #aggregate to genus level 
data <- transform_sample_counts(data, function(x) x/ sum(x)) # convert to relative abundance
data <- filter_taxa(data, function(x) mean(x) > 0.01e-2, TRUE) # mean of at least 0.01 
tax <- as(otu_table(data), "matrix")
tax_dist <- MiSPU::GUniFrac(otu.tab = tax, tree = phy_tree(data), alpha = 0.5)$GUniF[,,1]

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

# constructing ordinations  
tax_ord <- metaMDS(comm = tax_dist, try = 50, engine = "isoMDS")
met_ord <- metaMDS(comm = as.matrix(met_dist), try = 50, engine = "isoMDS")

results <- list(tax_ord = tax_ord, met_ord = met_ord)

saveRDS(results, file = opt$output)

