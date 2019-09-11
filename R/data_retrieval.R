# Updated 09/10/19
# Quang Nguyen
# This script retrieves all data files data files according requirements and export as phyloseq object 

library(purrr)
library(dplyr)
library(phyloseq)
library(optparse)

option_list <- list(
  make_option("--input", help = "Input file for data loading"),
  make_option("--output", help = "Output file for data loading"),
  make_option("--time", help = "Data time point to extract"),
  make_option("--metab_type", help = "Data type for metabolite data, can be 'tar' or 'untar'"),
  make_option("--tax_type", help = "Data type for taxonomic data, so far only '16S' is supported"),
  make_option("--fasttree_dir", help = "Getting where FastTree is installed")
)

opt <- parse_args(OptionParser(option_list = option_list))

#' @title Loading data from NHBCS sources  
#' @description This function takes a .csv file of appropriate directories and load the corresponding types  
#' @param dir_file Directory file loaded onto memory 
#' @param time Time point (6W, 12M)
#' @param metab_type Type of metabolomics data (targeted and untargeted)
#' @param tax_type Type of microbiome data (so far only 16S)
load_data <- function(dir_file, time, metab_type, tax_type){
  print('Loading metabolomic data...')
  if (metab_type == "tar"){
    print(dir_file[dir_file$X == "metabo",]$directory)
    metab <- readRDS(file = as.character(dir_file[dir_file$X == "metabo",]$directory))
  } else if (metab_type == "untar"){
    metab <- readRDS(file = dir_file[dir_file$X == "metabobin",]$directory)
  }
  print("Loading taxonomic data...")
  if (tax_type == "16S"){
    tax <- readRDS(file = dir_file[dir_file$X == "taxall",]$directory)
  } 
  print("Loading ids...")
  if(Sys.info()['sysname'] == "Darwin"){
    crosswalk_path <- dir_file[dir_file$X == "mac_access",]$directory
  } else if (.Platform$OS.type == "windows" | Sys.info['sysname'] == "Linux"){
    crosswalk_path <- dir_file[dir_file$X == "pc_access",]$directory
  }
  # grabbing cross walk ids 
  met_idx <- read.csv(file = paste0(crosswalk_path,"nmr.",time,".paired.mbl.st.csv"), header = F, stringsAsFactors = F)$V1
  tax_idx <- read.csv(file = paste0(crosswalk_path,"mbl.", time, ".paired.nmr.st.csv"), header = F, stringsAsFactors = F)$V1
  print("Matching ids...")
  match_met <- match(met_idx, metab$TubeLabel) # matching crosswalk idx with currently available metabolomics data 
  idx_not_matched <- which(is.na(match_met)) # idx of samples that can't match to a metabolomics profile 
  tax_idx <- tax_idx[-idx_not_matched] # remove those from the taxonomic ids 
  match_met <- match_met[-idx_not_matched]
  match_tax <- match(tax_idx, rownames(tax)) # matching the remainder 
  idx_not_matched <- which(is.na(match_tax)) # check if remaining tax_ids there are those that can't be matched 
  if (length(idx_not_matched > 0)){
    print("Some idx in the cross walk can't be found in the large taxa table")
    match_tax <- match_tax[-idx_not_matched]
    match_met <- match_met[-idx_not_matched]
  }
  print(paste0("Remaining number of samples: ", length(match_tax), "/", 
               length(met_idx), "(", round((length(match_tax)/length(met_idx))* 100,2), 
               "%", ")"))
  tax <- tax[match_tax,]
  metab <- metab[match_met,]
  rownames(metab) <- rownames(tax)
  metab <- metab[,-c(1,2,3,4,5)]
  result <- list(tax = tax, metab = metab)
  return(result)
}

#' @title Loading supplemental data  
#' @description Loading other supplemental data to main met and tax data such as the phylogenetic tree, and the taxa table 
#' @param dir_file Directory .csv file loaded onto memory 
#' @param supp_type The type of files can be gleaned from the dir_file. Defaults to "treeall" and "taxaallkey"
load_supp <- function(dir_file, supp_type = c("treeall", "taxallkey")){
  output <- list()
  for(i in supp_type){
    if (i == "treeall"){
      dir.create("./temp/", showWarnings = F)
      if (file.exists("./temp/full_tree.tre" ==)){
        print("There is no tree file here")
        seqs <- readRDS(file = dir_file[dir_file$X == i,]$directory)
        names(seqs) <- paste0("SV", 1:length(seqs))
        alignment <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(seqs), anchor = NA, verbose = T)
        phagAlign <- phangorn::phyDat(methods::as(alignment, "matrix"), type = "DNA")
        phangorn::write.phyDat(phagAlign, file = "./temp/full_alignment.fasta", format = "fasta")
        command <- paste0(opt$fasttree_dir, " FastTree -gtr -nt ./temp/full_alignment.fasta > ./temp/full_tree.tre")
        print(command)
        system(command)
      }
      output[[i]] <- ape::read.tree(file = paste0("./temp/full_tree.tre"))
    }
    output[[i]] <- readRDS(file = dir_file[dir_file$X == i,]$directory)
  }
  return(output)
}

#####-----------------------------------------------------------------------------------------------------------------------
dir_file <- read.csv(file = opt$input, stringsAsFactors = F)
main_dat <- load_data(dir_file = dir_file, metab_type = opt$metab_type, tax_type = opt$tax_type, time = opt$time)
supp_dat <- load_supp(dir_file = dir_file)

phylo_obj <- phyloseq(otu_table(main_dat$tax, taxa_are_rows = F), sample_data(main_dat$metab), 
                      tax_table(supp_dat$taxallkey))

if (is.null(opt$output) == T){
  dir.create("./data/", showWarnings = FALSE)
  output_dir <- paste0("./data/", time, "_", tax_type, "_", metab_type, ".rds")  
  saveRDS(phylo_obj, file = output_dir)
} else {
  saveRDS(phylo_obj, file = opt$output)
}




