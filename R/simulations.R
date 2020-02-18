library(tidyverse)
library(caret)
library(rsample)
library(glue)
library(phyloseq)
library(purrr)

if (Sys.info()['sysname'] == "Darwin"){
  if (file.exists("/Volumes/rc/Lab/")){
    dir <- "/Volumes/rc/Lab/QNguyen/ResultsFiles"
  } else if (file.exists("/Volumes/rc-1/Lab/")){
    dir <- "/Volumes/rc-1/Lab/QNguyen/ResultsFiles"
  }
} else if (Sys.info()['sysname'] == "Linux"){
  if (Sys.info()['nodename'] == "discovery7.dartmouth.edu"){
    dir <- "/dartfs-hpc/rc/lab/H/HoenA/Lab/QNguyen/ResultsFiles"
  } else {
    dir <- "/mnt/HoenLab/Lab/QNguyen/ResultsFiles"
  }
}

data <- readRDS(file = glue("{dir}/data/processed_12M_tar_phyloseq_obj.rds", dir = dir))
tax <- as(otu_table(data),"matrix")
met <- sample_data(data)
bts <- bootstraps(data = tax, times = 500)

dim(rnorm(72) %*% diag(rbinom(72, size = 1, prob = 0.5)))



