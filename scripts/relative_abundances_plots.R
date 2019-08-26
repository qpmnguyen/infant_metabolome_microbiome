library(MCTools)
library(phyloseq)
library(compositions)


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