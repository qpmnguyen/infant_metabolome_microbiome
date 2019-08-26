# Simple train/test split with 10 fold cross validation for parameter tuning  
library(MCTools)
library()
tar_12M <- readRDS(file = "./data/12M_tarNMR_clr_tax.rds")
tar_6W <- readRDS(file = "./data/6W_tarNMR_clr_tax.rds")
untar_12M <- readRDS(file = "./data/12M_untarNMR_clr_tax.rds")
untar_6W <- readRDS(file = "./data/6W_untarNMR_clr_tax.rds")

