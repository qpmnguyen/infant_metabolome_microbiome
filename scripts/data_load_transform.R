library(MCTools)
library(dplyr)
library(compositions)
library(SpiecEasi)
library(phyloseq)
source("./scripts/utils.R")

grab_data <- function(time_point){
  # data loading  
  data <- load_main(file = "./data/data_directory.csv")
  if (time_point == "12M"){
    tax <- data$tax.12M
    met <- data$metabo.12M
  } else if (time_point == "6W"){
    tax <- data$tax.6W
    met <- data$metabo.6W
  }
  
  key <- data$tax.key
  key <- key_processing(key, type = "NA")
  
  # grabbing binned data
  mblid <- rownames(tax)
  nmr_binned <- readRDS("//AFS/northstar/users/m/margk.collab/AnneHoen/Lab/QNguyen/SourceFiles/NMR/nmr_binned_ST_08May2018.rds")
  nmr_proc <- nmr_binned %>% filter(TimePeriod == time_point)
  idx_match <- match(mblid, nmr_binned$mblid)
  nmr_binned <- nmr_binned[idx_match,]
  
  nmr_binned <- nmr_binned[,-c(1,2,3,4,5,6)]
  rownames(nmr_binned) <- mblid
  nmr_key <- colnames(nmr_binned)
  colnames(nmr_binned) <- paste0("B", 1:ncol(nmr_binned))
  nmr_binned <- as.matrix(nmr_binned)
  nmr_binned <- unclass(acomp(nmr_binned)) # normalize to composition format  
  
  # processing and transformations  
  # nmr_binned <- sqrt(asin(nmr_binned))
  # met <- log(met + 1)
  
  phylo_obj <- phyloseq(otu_table(tax, taxa_are_rows = F),
                        tax_table(key))
  otu_table(phylo_obj) <- otu_table(phylo_obj) + 1 #pseudo count to everything 
  phylo_obj <- tax_glom(phylo_obj, taxrank = "Genus") #aggregate to genus level 
  phylo_obj <- prune_taxa(taxa_sums(phylo_obj) > 0, phylo_obj) #remove taxa with all 0s  
  
  tax <- as(otu_table(phylo_obj), "matrix")
  
  #saveRDS(list(tax = tax, met = met),file = "./data/raw_data.rds")
  # composition normalization and clr transformation
  tax <- unclass(acomp(tax))
  tax <- unclass(compositions::clr(tax))
  
  # targeted and untargeted
  untar <- list(tax = tax, met = nmr_binned)
  tar <- list(tax = tax, met = met)
  
  name_untar <- paste0(time_point, "_untarNMR_clr_tax.rds")
  name_tar <- paste0(time_point, "_tarNMR_clr_tax.rds")
  saveRDS(untar, paste0("./data/", name_untar))
  saveRDS(tar, paste0("./data/", name_tar))
}

grab_data(time_point = "6W")
grab_data(time_point = "12M")


