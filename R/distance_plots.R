generate_pair_plot <- function(dist_tax, dist_met, rank, xlab, ylab){
  dist_tax <- as.matrix(dist_tax)
  dist_met <- as.matrix(dist_met)
  dist_tax[lower.tri(dist_tax)] <- 0
  dist_met[lower.tri(dist_met)] <- 0
  dist_tax <- melt(dist_tax) %>% filter(value > 0)
  dist_met <- melt(dist_met) %>% filter(value > 0)
  if (rank == T){
    dist_tax$value <- rank(dist_tax$value)
    xlab <- paste("Ranked", xlab, "distance")
  } else {
    xlab <- paste(xlab, "Distance")
  }
  ylab <- paste(ylab, "Distance")
  plot_dat <- data.frame(cbind(dist_tax$value, dist_met$value))
  colnames(plot_dat) <- c("tax", "met")
  plt <- ggplot(plot_dat, aes(x = tax, y = met)) + 
    geom_point() + 
    stat_binhex(bins = 80) + 
    geom_hline(yintercept = quantile(dist_met$value, 0.5), colour ="red") + 
    labs(x = xlab, y = ylab) + theme_bw()
  
  return(plt)
}

(p <- generate_pair_plot(tax_gunifrac_dist, met_PCs_manhattan_dist, rank = F, xlab = "Taxonomic", ylab = "Metabolite") + labs(title = "Targeted"))
(p2 <- generate_pair_plot(tax_gunifrac_dist, untarmet_manhattahn, rank = F, xlab = "Taxonomic", ylab = "Metabolite") + labs(title = "Untargeted"))
(joined <- grid.arrange(grobs = list(p = p, p2 = p2), nrow = 1, col = 2))
ggsave(plot = joined, filename = "./docs/distance_unranked_comparison_12M.png", device = "png", width = 8, height = 5, units = "in")