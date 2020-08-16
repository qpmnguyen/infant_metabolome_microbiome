library(drake)
library(tidyverse)
library(glue)
library(rlang)

source("R/utils.R")
source("R/analysis_spearman_corr.R")
source("R/analysis_sparse_cca.R")

plan <- drake::drake_plan(
  max_expand = 1, 
  data = target(
    readRDS(file = get_dir(time, mettype, data_stage = "processed", sens = T)),
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
    sparse_cca_main(data, n_boot = 5000, n_perm = 1000),
    transform = map(data)
  )
)

drake::predict_runtime(plan)
#vis_drake_graph(plan)
make(plan)
