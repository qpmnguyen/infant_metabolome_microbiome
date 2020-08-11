library(drake)
library(tidyverse)
library(glue)
library(rlang)

source("R/utils.R")
source("R/analysis_spearman_corr.R")
source("R/analysis_sparse_cca.R")

plan <- drake::drake_plan(
  data = target(
    readRDS(file = get_dir(time, mettype, data_stage = "raw", sens = T)),
    transform = cross(
      time = c("12M", "6W"),
      mettype = c("tar", "untar")
    )
  ), 
  spearman = target(
    spearman_main(data),
    transform = map(data)
  ),
  scca = target(
    scca_main(data),
    transform = map(data)
  )
)

vis_drake_graph(plan)
