library(tidyverse)
library(stringr)
library(ggsci)
library(pheatmap)
library(viridis)
library(ggplotify)
library(glue)
library(patchwork)

result <- readRDS(file = "picrust2/sig_test.rds")
sig <- result %>% filter(qval <= 0.1)
names <- unique(sig$met)
times <- c("6w", "12m")

plt <- sig %>% ggplot(aes(x = str_wrap(description, 20), y = est, col = met)) + 
  geom_point(size = 5, position = position_dodge(width = 1)) + 
  geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 1)) + 
  scale_color_npg() +
  theme_bw(base_size = 20) +
  facet_wrap(~timepoint, scales = "fixed") + 
  coord_flip() + 
  labs(y = "Spearman Correlation", x = "Pathway Description", col = "Metabolites")
plt
ggsave(plt, filename = "docs/publication_figures/picrust_sig.png", dpi = 300, width = 15, height = 12)

ggplot(result, aes(x = met, y = path, fill = est)) + geom_tile()

## Plotting the correlation matrix  



plot_correlation <- function(timepoint){
  cor_obj <- readRDS(file = glue("picrust2/correlation_{timepoint}.rds", timepoint = timepoint))
  r <- cor_obj$r
  adj_p <- cor_obj$p_adj
  adj_p <- ifelse(adj_p <= 0.05, 1, 0)
  
  r <- r * adj_p
  r[r == 0] <- NA
  
  map <- pheatmap(
    mat = t(r),
    color = viridis(40),
    drop_levels = TRUE,
    fontsize_row = 13,
    fontsize_col = 5,
    cluster_rows = F, 
    cluster_cols = F,
    na_col = "beige", show_colnames = FALSE
  )
  map <- ggplotify::as.ggplot(map)
  return(map)
}

htmap_12m <- plot_correlation("12m")
htmap_6w <- plot_correlation("6w")
htmap_12m <- htmap_12m + labs(subtitle = "12-month samples") 
htmap_6w <- htmap_6w + labs(subtitle = "6-week samples")

comb <- htmap_6w/htmap_12m + plot_annotation(tag_levels = "A")
#/ htmap_6w + plot_annotation(tag_levels = "A", tag_suffix = "")

ggsave(comb, filename = "docs/publication_figures/picrust2_heatmap.png", dpi = 300, 
       width = 14, height = 12)


