library(phyloseq)
library(data.table)
library(stringr)
library(tidyverse)
library(glue)
library(patchwork)

dir <- "//dartfs-hpc/rc/Lab/H/HoenA/Lab/QNguyen/"
generate_contrib_plot <- function(timepoint){
  asv_contrib <- fread(glue("picrust2/temp/pathway_out_{timepoint}/path_abun_contrib.tsv.gz", timepoint = timepoint)) %>% 
    as_tibble()
  sig <- readRDS(file = "picrust2/sig_test.rds")
  raw_asv <- readRDS(file = paste0(dir, glue("ResultsFiles/data/raw_{timepoint}_tar_phyloseq_obj.rds", timepoint = timepoint)))
  table <- raw_asv %>% tax_table() %>% as("matrix") %>% as.data.frame() %>% 
    rownames_to_column(var = "taxon")
  
  sig_6w <- sig %>% filter(timepoint == timepoint) %>% filter(qval <= 0.05) %>% 
    rename("function" = "path")
  unq_path <- unique(sig_6w %>% pull(`function`))
  asv_contrib <- left_join(asv_contrib, table, by = "taxon")
  asv_contrib <- asv_contrib %>% filter(`function` %in% unq_path)
  asv_contrib <- left_join(asv_contrib, sig_6w %>% select(`function`, description), by = "function")
  asv_contrib <- asv_contrib %>% filter(!is.na(Genus))
  contrib <- asv_contrib %>% group_by(`function`, description, Genus) %>% 
    summarise(abun = sum(taxon_function_abun)) %>%
    group_by(`function`) %>% 
    mutate(total = sum(abun)) %>% 
    ungroup() %>% mutate(rel_abun = abun / total) %>% 
    group_by(`function`, description) %>% top_n(5, wt = rel_abun) %>% ungroup()
  
  
  contrib_plot <- ggplot(contrib, aes(x = Genus, y = rel_abun, fill = str_wrap(description, 8))) + geom_bar(stat = "identity") +
    facet_wrap(~description, scales = "free", 
               labeller = labeller(description = label_wrap_gen()),
               nrow = 3) + 
    coord_flip() + scale_fill_viridis_d() + 
    theme_bw(base_size = 15) + theme(legend.position = "None") + 
    labs(x = "Genus", y = "Predicted Relative Contribution")
  return(contrib_plot)
}

contrib_6w <- generate_contrib_plot("6w")
contrib_12m <- generate_contrib_plot("12m")
contrib_6w <- contrib_6w + theme(axis.title.x = element_blank())
contrib_combine <- contrib_6w/contrib_12m + plot_annotation(tag_levels = "A")
ggsave(contrib_combine, filename = "docs/publication_figures/picrust2_contrib.png", dpi = 300, 
       width = 21, height = 13)
