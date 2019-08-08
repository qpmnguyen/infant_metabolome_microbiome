library(MCTools)
library(dplyr)
library(compositions)
library(SpiecEasi)
library(phyloseq)
source("./scripts/utils.R")
# data loading
data <- load_main(file = "./data/data_directory.csv")
tax <- data$tax.12M
met <- data$metabo.12M
key <- data$tax.key
key <- key_processing(key, type = "NA")

# grabbing binned data
mblid <- rownames(tax)
nmr_binned <- readRDS("//AFS/northstar/users/m/margk.collab/AnneHoen/Lab/QNguyen/SourceFiles/NMR/nmr_binned_ST_08May2018.rds")
nmr_proc <- nmr_binned %>% filter(TimePeriod == "12M")
idx_match <- match(mblid, nmr_binned$mblid)
nmr_binned <- nmr_binned[idx_match,]

nmr_binned <- nmr_binned[,-c(1,2,3,4,5,6)]
rownames(nmr_binned) <- mblid
nmr_key <- colnames(nmr_binned)
colnames(nmr_binned) <- paste0("B", 1:ncol(nmr_binned))
nmr_binned <- as.matrix(nmr_binned)
nmr_binned <- unclass(acomp(nmr_binned)) # normalize to composition format  

# processing and transformations  
#nmr_binned <- sqrt(asin(nmr_binned))
#met <- log(met + 1)

phylo_obj <- phyloseq(otu_table(tax, taxa_are_rows = F),
                      tax_table(key))
otu_table(phylo_obj) <- otu_table(phylo_obj) + 1
phylo_obj <- tax_glom(phylo_obj, taxrank = "Genus")
phylo_obj <- prune_taxa(taxa_sums(phylo_obj) > 0, phylo_obj)

tax <- as(otu_table(phylo_obj), "matrix")

#saveRDS(list(tax = tax, met = met),file = "./data/raw_data.rds")

tax <- unclass(acomp(tax))
tax <- unclass(compositions::clr(tax))

untar <- list(tax = tax, met = nmr_binned)
tar <- list(tax = tax, met = met)

saveRDS(untar, "12M_untarNMR_clr_tax.rds")
saveRDS(tar, "12M_tarNMR_clr_tax.rds")
