# Script to retrieve unq_ids used for data matching
# Retrieving consensus ASV sequences for usage in picrust2
# Quang Nguyen 
# unq ids not shared in repo 
library(phyloseq)
library(tidyverse)
library(biomformat)

# get directory
save_dir <- "//dartfs-hpc/rc/Lab/H/HoenA/Lab/QNguyen/ResultsFiles/data"

# Get raw data 
samples_6w <- readRDS(file = paste0(save_dir, "/raw_6W_tar_phyloseq_obj.rds"))
samples_12m <- readRDS(file = paste0(save_dir, "/raw_12M_tar_phyloseq_obj.rds"))

# Retrieve sample names for quality control statistics
sample_names(samples_6w) %>% 
  write.table(file = "output/samp_ids_6w.csv", row.names = FALSE, col.names = c('ids'))
sample_names(samples_12m) %>% 
  write.table(file = "output/samp_ids_12m.csv", row.names = FALSE, col.names = c('ids'))

# Filter samples and write to biom file
samples_6w <- filter_taxa(samples_6w, function(x) sum(x > 0) >= 0.1*length(x), TRUE)
biom_6w <- make_biom(
  data = otu_table(samples_6w, taxa_are_rows = FALSE)
)
write_biom(biom_6w, biom_file = "output/biom_6w.biom")


samples_12m <- filter_taxa(samples_12m, function(x) sum(x > 0) >= 0.1 * length(x), TRUE)
biom_12m <- make_biom(
  data = otu_table(samples_12m, taxa_are_rows = FALSE)
)
write_biom(biom_12m, biom_file = "output/biom_12m.biom")

# Get sequences used for each data type
seq_12m <- readRDS(file = "//dartfs-hpc/rc/Lab/H/HoenA/Lab/QNguyen/ResultsFiles/seqtab.nochim.colnames_12M_ST.rds")

