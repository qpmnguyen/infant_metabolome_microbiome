library(drake)
library(tidyverse)
library(glue)
library(rlang)
library(stringr)

source("R/utils.R")
source("R/analysis_constructing_distance_matrices.R")
source("R/analysis_procrustes.R")

plan <- drake::drake_plan(
  data = target(
    readRDS(file = get_dir(time, mettype, data_stage = "raw", sens = T)),
    transform = cross(
      time = c("12M", "6W"),
      mettype = c("tar", "untar")
    )
  ), 
  dist_mat = target(
    dist_main(data = data, mettype = mettype, dist_tax = dist_tax, dist_met = dist_met),
    transform = cross(data, 
                      dist_tax = c("euclidean", "gunifrac"),
                      dist_met = c("euclidean", "manhattan"))
  ),
  procrustes = target(
    procrustes_main(dist_mat),
    transform = map(dist_mat)
  ),
  output = target({
      dir.create("drake_procrustes/")
      saveRDS(procrustes, file = glue("drake_procrustes/{time}_{mettype}_{dist_tax}_{dist_met}.rds",
                                 time = time, mettype = mettype, dist_tax = dist_tax, dist_met = dist_met))
  }, transform = map(procrustes), format = "file")
)
make(plan)
