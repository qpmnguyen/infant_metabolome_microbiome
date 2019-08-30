library(purrr)
library(dplyr)



load_data <- function(dir_file, time, metab_type, tax_type){
  print('Loading metabolomic data...')
  if (metab_type == "tar"){
    metab <- readRDS(file = dir_file[dir_file$X == "metabo",]$directory)
  } else if (metab_type == "untar"){
    metab <- readRDS(file = dir_file[dir_file$X == "metabobin",]$directory)
  }
  print("Loading taxonomic data...")
  if (tax_type == "16S"){
    tax <- readRDS(file = dir_file[dir_file$X == "taxall",]$directory)
  } else if (tax_type == "WGS"){
    taxa <- read.table(file = "/Volumes/rc-1/Lab/Share/MBLmetagenomics/metaphlan2output/merged_MC/merged_abundance_table_new.txt", sep = "\t")
    taxa <- t(taxa)
    rownames(taxa) <- taxa[,1]
    colnames(taxa) <- taxa[1,]
    taxa <- taxa[-1,]
    samp_names <- taxa[,"ID"]
    samp_names <- as.vector(sapply(samp_names, function(x){
      strsplit(x, "_")[[1]][1]
    }))
    taxa <- taxa[,-1]
    taxa <- apply(taxa, 2,as.numeric)
    rownames(taxa) <- samp_names
  }
  print("Loading ids...")
  if (time == "12M"){
    met_idx <- read.csv(file = "/Volumes/rc-1/Lab/Share/SourceFiles/IDcrosswalk/nmr.12M.paired.mbl.st.csv", header = F, stringsAsFactors = F)
    tax_idx <- read.csv(file = "/Volumes/rc-1/Lab/Share/SourceFiles/IDcrosswalk/mbl.12M.paired.nmr.st.csv", header = F, stringsAsFactors = F)
  } else if (time == "6W"){
    met_idx <- read.csv(file = "/Volumes/rc-1/lab/Share/SourceFiles/IDcrosswalk/nmr.6W.paired.mbl.st.csv", header = F, stringsAsFactors = F)
    tax_idx <- read.csv(file = "/Volumes/rc-1/lab/Share/SourceFiles/IDcrosswalk/mbl.6W.paired.nmr.st.csv", header = F, stringsAsFactors = F)
  }
  match_met <- match(met_idx$V1, metab$TubeLabel)
  match_met <- match_met[!is.na(match_met)]
  match_tax <- match(tax_idx$V1, rownames(taxa))
  match_tax <- match_tax[!is.na(match_tax)]
  tax <- tax[match_tax,]
  metab <- metab[match_met,]
  metab <- metab[,-c(1,2,3,4,5)]
  result <- list(tax = tax, metab = metab)
  return(result)
}

dir_file <- read.csv(file = "./data/data_directory.csv", stringsAsFactors = F)
load_data(dir_file = dir_file, time = "12M", metab_type = )