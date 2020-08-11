# Script used to retrieve data from processed data files with samples matching between two 
# timepoints and data types
library(phyloseq)
library(glue)
library(stringr)

save_dir <- "//dartfs-hpc/rc/Lab/H/HoenA/Lab/QNguyen/ResultsFiles/data"

save_data <- function(mettype = "tar", data_stage){
  data_files <- list.files(path = "//dartfs-hpc/rc/Lab/H/HoenA/Lab/QNguyen/ResultsFiles/data/",
                           full.names = T, pattern = data_stage)
  print(data_files)
  untar_idx <- str_detect(data_files, "untar")
  if(mettype == "tar"){
    data_files <- data_files[!untar_idx]
  } else {
    data_files <- data_files[untar_idx]
  }
  data_6w <- readRDS(file = data_files[str_detect(data_files, "6W")])
  data_12m <- readRDS(file = data_files[str_detect(data_files, "12M")])
  common_samples <- intersect(sample_names(data_6w), sample_names(data_12m))
  data_6w <- subset_samples(data_6w, sample_names(data_6w) %in% common_samples)
  data_12m <- subset_samples(data_12m, sample_names(data_12m) %in% common_samples)
  saveRDS(data_6w, file = glue("{dir}/{stage}_common_{time}_{tartype}_phyloseq_obj.rds", stage = data_stage,
                               dir = save_dir, time = "6W", tartype = mettype))
  saveRDS(data_12m, file = glue("{dir}/{stage}_common_{time}_{tartype}_phyloseq_obj.rds", stage = data_stage,
                               dir = save_dir, time = "12M", tartype = mettype))
  
}

for (i in c("tar", "untar")){
  for (j in c("processed", "raw"))
    save_data(mettype = i, data_stage = j)
}
