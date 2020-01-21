# DECONTAM and preliminary analyses  
# Quang Nguyen
# Updated Jan 20th 2020
library(decontam)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(vegan)
library(MiSPU)
library(ape)
library(cowplot)

decontam_fit <- function(time){
  ###---Step 1: Load data files and concentration ---###
  data_repo <- read.csv(file = "/mnt/HoenLab/Lab/QNguyen/ResultsFiles/data_directory.csv", stringsAsFactors = F)
  asv_table <- readRDS(file = data_repo$directory[data_repo$key == paste0("tax",time)])
  tax_table <- readRDS(file = data_repo$directory[data_repo$key == paste0("tax",time,"key")])
  tree <- read.tree(file = paste0("/mnt/HoenLab/Lab/QNguyen/ResultsFiles/",time,"_full_tree.tre"))
  inventory <- read.csv(file = paste0("/mnt/HoenLab",data_repo$directory[data_repo$key == "tax_inventory"]), stringsAsFactors = F)
  ###---Step 2: Matching samples from asv_table to inventory and extract DNA concentrations using mblids---###
  inventory <- inventory %>% arrange(sequence_count) %>% filter(TimePeriod == time & participant_type == "child" & sample_type == "stool")
  matching_idx <- match(inventory$MBL_ID, rownames(asv_table))
  inventory <- inventory[-which(is.na(matching_idx)),]
  asv_table <- asv_table[na.omit(matching_idx),]
  info <- as.data.frame(cbind(dna_conc = inventory$Sample_Concentration, batch = inventory$Batch))
  rownames(info) <- rownames(asv_table)
  df <- phyloseq(otu_table(asv_table,taxa_are_rows = F), tax_table(tax_table), sample_data(info), phy_tree = tree)
  ###---Step 3: Using decontam to remove contaminants using default parameters---###
  contamdf.freq <- isContaminant(df, conc = "dna_conc", batch = "batch")
  output = list(df = df, contamdf.freq = contamdf.freq)
  return(output)
}


diagnostics <- function(df, contamdf.freq, distance_metric, time){
  contam_idx <- which(contamdf.freq$contaminant == T)
  # frequency plot
  freq_plot <- plot_frequency(df, taxa_names(df)[contam_idx], conc = "dna_conc")
  ggsave(freq_plot, file = paste0("docs/decontam_freq_plot", "_", time, ".jpg"), dpi = 300)
  
  # before constructing distance metric convert everything to relative abundances  
  df <- transform_sample_counts(df,function(x) x/sum(x))
  contam_names <- taxa_names(df)[-contam_idx]
  new_df <- prune_taxa(contam_names, df)
  if (distance_metric == "bray"){
    org_dist <- distance(df, method = "bray")
    org_MDS <- ordinate(df, "NMDS", distance = org_dist)
    new_dist <- distance(new_df, method = "bray")
    new_MDS <- ordinate(new_df, method = "NMDS", distance = new_dist)
    p1 <- qplot(x = org_MDS$points[,1], y = org_MDS$points[,2], col = I("steelblue"),
                xlab = "NMDS1", ylab = "NMDS2", main = "Unaltered data")
    p2 <- qplot(x = new_MDS$points[,1], y = new_MDS$points[,2], col = I("steelblue"), 
                xlab = "NMDS1", ylab = "NMDS2", main = "Contaminants removed")
    procrustes <- protest(org_MDS, new_MDS)
    label <- c(rep("original", nrow(otu_table(df))), rep("decontam", nrow(otu_table(df))))
    super_imposed <- data.frame(rbind(procrustes$X, procrustes$Yrot))
    title <- paste0("SoS: ", round(procrustes$ss, 2), " (Null SoS: ", round(procrustes$t0,2),", p = ",procrustes$signif,")")
    p3 <- qplot(x = super_imposed[,1], y = super_imposed[,2], col = label, 
                xlab = "NMDS1", ylab = "NMDS2", main = title)
    comb_plot <- plot_grid(p1,p2,p3, nrow = 2, ncol = 2)
    ggsave(comb_plot, paste0("docs/decontam_ord_plot_", time, "_bray.jpg"))
  } else if (distance_metric == "GUnifrac"){
    invisible()
  }
}
