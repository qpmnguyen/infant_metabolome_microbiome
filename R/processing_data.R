# setting correct directory for this path and install required packages
library(optparse)
library(phyloseq)
library(compositions)
library(plyr)
library(glue)

option_list <- list(
  make_option("--input_folder", help="Folder where the raw data is stored in phyloseq form. Files should have '_phyloseq_obj.rds' ending and 'raw' beginnings"),
  make_option("--output_folder", default = NULL, help = "Folder where the output file will be stored, defaults to where input files are found")
)

opt <- parse_args(OptionParser(option_list = option_list))


if (is.null(opt$output_folder)){
  opt$output_folder <- opt$input_folder
  print("Outputing results to the same folder as input")
}


#' Function to process data 
#' @param data object loaded onto memory of phyloseq object type
#' @param metabtype "tar" or "untar"
#' @param timepoint "6W" or "12M" 
processing <- function(data, metabtype, timepoint){
  data <- filter_taxa(data, function(x) count(x > 0)[2,2] >= 0.1*length(x), TRUE) # asvs as to be in at least 10% of samples 
  otu_table(data) <- otu_table(data) + 1 # adding pseudo count of 1
  data <- tax_glom(data, taxrank = "Genus") #aggregate to genus level 
  data <- transform_sample_counts(data, function(x) x/ sum(x)) # convert to relative abundance
  data <- filter_taxa(data, function(x) mean(x) > 0.005e-2, TRUE) # mean of at least 0.005 percent (ref. Bokulich 2013)
  tax <- as(otu_table(data), "matrix")
  tax <- unclass(clr(tax)) # centerd log_ratio transformation 
  if (metabtype == "tar"){
    met <- as(sample_data(data), "matrix")
    met <- log(met + 1) # log x+1 transformation
  } else {
    met <- as(sample_data(data), "matrix")
    met <- unclass(acomp(met)) # renormalize to relative abundances
    met <- asin(sqrt(met)) # arcsine square root transformation 
    colnames(met) <- paste0("B",1:ncol(met))
  }
  tree <- phy_tree(data)
  table <- tax_table(data)
  met <- as.data.frame(met)
  new_data <- phyloseq(otu_table(tax, taxa_are_rows = F), 
                       sample_data(met), 
                       tax_table(table), 
                       phy_tree(tree))
  rm(data)
  saveRDS(new_data, file = glue("{folder}/processed_{time}_{metab}_phyloseq_obj.rds", folder = opt$output_folder,
                                time = timepoint, metab = metabtype))
}

TIME <- c("6W", "12M")
METAB <- c("tar", "untar")

for (i in TIME){
  for (j in METAB){
    print(i)
    print(j)
    data <- readRDS(file = glue("{folder}/raw_{time}_{metab}_phyloseq_obj.rds", folder = opt$input_folder, 
                                time = i, metab = j))
    processing(data = data, timepoint = i, metabtype = j)
  }
}



