library(drake)

source("R/utils.R")
source("R/plots_correlation.R")
source("R/plots_ordination.R")

correlation_plan <- drake::drake_plan(
  spearman_dat = target(
    readRDS(file = glue("drake_procrustes/{time}_{mettype}_spearman.rds",time = time, mettype = mettype)),
    transform = cross(
      time = c("6W","12M"),
      mettype = c("tar", "untar")    
    )
  ), 
  scca_dat = target(
    readRDS(file = glue("drake_procrustes/{time}_{mettype}_scca.rds",time = time, mettype = mettype)),
    transform = cross(
      time = c("6W","12M"),
      mettype = c("tar", "untar")    
    )
  ),
  spearman = target(
    correlation_plot(data = spearman_dat),
    transform = map(spearman_dat)
  ),
  scca = target(
    scca_plot(data = scca_dat, spearman_res = spearman_dat),
    transform = map(spearman_dat, scca_dat)
  ),
  combined = target(
    combine_plots_correlation(scca_plot = scca, spearman_plot = spearman, scca_obj = scca_dat),
    transform = combine(spearman, scca, scca_dat, .by = c(time, mettype))
  ),
  output_combined_plots = target({
    dir.create("drake/plots/")
    saveRDS(combined, file = glue("drake/plots/{time}_{mettype}_cor_combined_plt.rds",
                              time = time, mettype = mettype))
    ggsave(combined, filename = glue("drake/plots/{time}_{mettype}_cor_combined_plt.png",
                                     time = time, mettype = mettype), dpi = 300, 
           width = 15, height = 12)
  }, transform = map(combined), format = "file")
)

ordination_plan <- drake_plan(
  analysis = target(
    readRDS(file = glue("drake/{time}_{mettype}_{taxdist}_euclidean.rds", 
                        time = time, mettype = mettype, taxdist = taxdist)),
    transform = cross(
      time = c("6W", "12M"),
      mettype = c("tar", "untar"),
      taxdist = c("gunifrac", "euclidean"),
    )
  ),
  procrustes_plot = target(
    plotting_proc(analysis$proc_test),
    transform = map(analysis)
  ),
  combine_plot = target(
    combine_plot_ord(procrustes_plot),
    transform = combine(procrustes_plot, .by = c(mettype))
  ),
  export = target({
      dir.create("drake/plots/")
      saveRDS(combine_plot, file = glue("drake/plots/{mettype}_procrustes_combined_plt.rds",
                                    time = time, mettype = mettype))
      ggsave(combine_plot, filename = glue("drake/plots/{mettype}_procrustes_combined_plt.png",
                                      time = time, mettype = mettype), dpi = 300, width = 15, height = 12)
    }, transform = map(combine_plot), format = "file")
)

drake::vis_drake_graph(correlation_plan)
make(correlation_plan)


drake::vis_drake_graph(ordination_plan)
make(ordination_plan)
