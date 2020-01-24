library(dplyr)

dir_file <- read.csv(file = "./data/data_directory.csv", stringsAsFactors = F)

met_data <- readRDS(file = dir_file[dir_file$X == "metabo",]$directory)

met_idx <- readRDS(file = "/Volumes/rc/Lab/Share/SourceFiles/IDcrosswalk/nmr.long.rds")
if (met_idx )
met_idx <- readRDS(file = "//dartfs-hpc/rc/lab/H/HoenA/Lab/Share/SourceFiles/IDcrosswalk/nmr.12M.st.rds")
length(met_idx$TimePeriod[met_idx$TimePeriod == "6W"])
length(met_data$TimePeriod[met_data$TimePeriod == "6W"])

met_idx$nmr[met_idx$TimePeriod == "6W"]