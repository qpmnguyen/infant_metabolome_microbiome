library(MCTools)
library(phyloseq)
library(phytools)
source("./scripts/utils.R")
library(compositions)

data <- load_main("./data/data_directory.csv")
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

new_tax <- as(otu_table(phylo_obj), "matrix")



new_tax <- unclass(clr(new_tax))

COV_mult <- function(df){
  S <- cov(df)
  x_bar <- apply(df, 2, mean)
  num <- t(x_bar) %*% S %*% x_bar
  denom <- (t(x_bar) %*% x_bar)^2
  cov <- sqrt(num/denom)
  return(as.numeric(cov))
}

cov_prop <- function(p){
  p_bar <- mean(p)
  sum((p - p_bar)^2)/length(p)
}

cov <- apply(new_tax,2,function(x){
  mean(x)/sd(x)
})

plot(cov)

cov2 <- apply(met,2, function(x){
  mean(x)/sd(x)
})
plot(cov2)


variation <- COV_mult(new_tax)
var_met <- COV_mult(met)

COV_mult(log(met + 1))
hist(log(met + 1))

COV_mult(new_tax)
