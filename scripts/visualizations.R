library(ggplot2)
library(phyloseq)
library(reshape2)
library(compositions)
library(MCTools)

source("./scripts/utils.R")
tar <- readRDS(file = './data/tarNMR_clr_tax.rds')
untar <- readRDS(file = "./data/untarNMR_clr_tax.rds")
data <- load_main(file = "./data/data_directory.csv")
tax <- data$tax.12M
key <- data$tax.key
key <- key_processing(key, type = "NA")
phylo_obj <- phyloseq(otu_table(tax, taxa_are_rows = F),
                      tax_table(key))
phylo_obj <- transform_sample_counts(phylo_obj, function(x) x/sum(x))
phylo_obj <- tax_glom(phylo_obj, taxrank = "Family")
phylo_obj <- prune_taxa(taxa_sums(phylo_obj) > 0, phylo_obj)


ggplot(melt(otu_table(phylo_obj)), aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = 'identity', color = "black") + 
         theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")

tar_met <- tar$met
untar_met <- untar$met


ggsave(plot = p1,device = "pdf", filename = "./family_ra.pdf")


tax <- unclass(clrInv(tax))
tar_met <- unclass(acomp(tar_met))

ggplot(melt(tar_met), aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity", color = "black") + theme(legend.position = "none")
ggplot(melt(untar_met), aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity", color = "black") + theme(legend.position = "none")
ggplot(melt(tax), aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity", color = "black") + theme(legend.position = "none")
