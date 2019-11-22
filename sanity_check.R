library(dplyr)

dir_file <- read.csv(file = "./data/data_directory.csv", stringsAsFactors = F)

met_data <- readRDS(file = dir_file[dir_file$X == "metabo",]$directory)

met_idx <- readRDS(file = "/Volumes/rc/Lab/Share/SourceFiles/IDcrosswalk/nmr.long.rds")
length(met_idx$TimePeriod[met_idx$TimePeriod == "6W"])
length(met_data$TimePeriod[met_data$TimePeriod == "6W"])
       