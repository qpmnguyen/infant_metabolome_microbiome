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
  ),
  output_spearman = target({
    dir.create("drake_procrustes/")
    saveRDS(spearman, file = glue("drake_procrustes/{time}_{mettype}_spearman.rds",
                                    time = time, mettype = mettype))
  }, transform = map(spearman), format = "file"),
  output_scca = target({
    dir.create("drake_procrustes/")
    saveRDS(scca, file = glue("drake_procrustes/{time}_{mettype}_scca.rds",
                                    time = time, mettype = mettype))
  }, transform = map(scca), format = "file")
)


#drake::predict_runtime(plan)
#vis_drake_graph(plan)
make(plan)
