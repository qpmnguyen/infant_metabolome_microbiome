library(tidyverse)
library(phyloseq)
library(table1)
library(kableExtra)
library(glue)
if (Sys.info()['sysname'] == "Darwin"){
  if (file.exists("/Volumes/rc/Lab/") == T){
    dir <- "/Volumes/rc/Lab"
  } else if (file.exists("/Volumes/rc-1/Lab/") == T){
    dir <- "/Volumes/rc-1/Lab"
  }
} 

# get initial reference files 
covariate_list <- read.csv(file = glue("{dir}/QNguyen/SourceFiles/NMR/req23April2019_ah_11nov2019_QN.csv", dir = dir), 
                           stringsAsFactors = F)
feeding_mode <- read.csv(file = glue("{dir}/Share/ResultsFiles/FeedingMode/FeedingModeBySampleAge_2020_02_11.csv", dir = dir), 
                         stringsAsFactors = F)
crosswalk_6W <- read.csv(file = glue("{dir}/Share/SourceFiles/IDcrosswalk/mbl.nmr.6W.paired.st.csv", dir = dir), stringsAsFactors = F, 
                      header = F)
crosswalk_12M <- read.csv(file = glue("{dir}/Share/SourceFiles/IDcrosswalk/mbl.nmr.12M.paired.st.csv", dir = dir), stringsAsFactors = F, 
                         header = F)
# get initial ids to match
six_weeks <- readRDS(file = glue("{dir}/QNguyen/ResultsFiles/data/processed_6W_tar_phyloseq_obj.rds", dir = dir))
ids_6W <- sample_names(six_weeks)
twelve_months <- readRDS(file = glue("{dir}/QNguyen/ResultsFiles/data/processed_12M_tar_phyloseq_obj.rds", dir = dir))
ids_12M <- sample_names(twelve_months)

# get initial covariates
cov_6W <- covariate_list %>% filter(unq_id %in% ids_6W) %>% select(unq_id, imrbthwghtg_all, bbymale, deliverytype, gestage_all, evercigpreg,
                                                      imrantibiotics, bby_race) %>% mutate(time = "6W")
cov_12M <- covariate_list %>% filter(unq_id %in% ids_12M) %>% select(unq_id, imrbthwghtg_all, bbymale, deliverytype, gestage_all, evercigpreg,
                                                                     imrantibiotics, bby_race) %>% mutate(time = "12M")
get_feeding_mode <- function(covariate, crosswalk, feeding){
  matching <- match(covariate$unq_id, crosswalk$V1)
  mblids <- crosswalk$V2[matching]
  fm_idx <- match(mblids, feeding$sampleID)
  return(as.vector(feeding$dFeedingMode[fm_idx]))
}

cov_6W <- cov_6W %>% mutate(feeding = get_feeding_mode(cov_6W, crosswalk_6W, feeding_mode))
cov_12M <- cov_12M %>% mutate(feeding = get_feeding_mode(cov_12M, crosswalk_12M, feeding_mode))

# process covariate table for export as table 1
cov <- rbind(cov_6W, cov_12M)
cov$bbymale <- factor(cov$bbymale, levels = c(1,0), labels = c("Male", "Female"))
cov$deliverytype <- factor(cov$deliverytype, levels = c(1,2), labels = c("Vaginal", "Cesarian"))
cov$time <- factor(cov$time, levels = c("6W", "12M"), labels = c("6W", "12M"))
cov$feeding <- factor(cov$feeding, levels = c(0,1,2,3), labels = c("Unknown", "Exclusively breastfed", "Exclusively formula fed", "Mixed"))
cov$imrantibiotics <- factor(cov$imrantibiotics, levels = c(0,1), labels = c("No", "Yes"))
cov$evercigpreg <- factor(cov$evercigpreg, levels = c(0,1), labels = c("No", "Yes"))
cov$bby_race <- factor(cov$bby_race, labels = c("Asian", "White", "Mixed"))
label(cov$bbymale) <- "Sex"
label(cov$deliverytype) <- "Delivery Mode"
label(cov$imrbthwghtg_all) <- "Birthweight"
label(cov$gestage_all) <- "Gesational Age"
label(cov$feeding) <- "Feeding Mode"
label(cov$imrantibiotics) <- "Infant antibiotics intake during hospitalization"
label(cov$bby_race) <- "Infant Race"
label(cov$evercigpreg) <- "Maternal smoking during pregnancy"
units(cov$imrbthwghtg_all) <- "grams"
units(cov$gestage_all) <- "Weeks"




tbl <- table1(~imrbthwghtg_all + bbymale + feeding + deliverytype + gestage_all + imrantibiotics + evercigpreg + bby_race | time,
       data = cov, overall = FALSE, topclass = "Rtable1-zebra")
tbl
