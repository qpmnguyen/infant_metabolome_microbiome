library(tidyverse)
library(phyloseq)
library(rlang)
library(patchwork)
library(glue)
library(data.table)
library(pheatmap)
library(Hmisc)

dir <- "//dartfs-hpc/rc/Lab/H/HoenA/Lab/QNguyen/"


get_correlation <- function(timepoint, met_id){
  # first, obtain path data, and rename samples for simplicity
  path_data <- readRDS(file = glue("picrust2/path_results/path_abun_{timepoint}_proc.rds", timepoint = timepoint)) %>% 
    as_tibble()
  colnames(path_data) <- c("pathways", paste0("samp", 1:(ncol(path_data)-1)))
  
  # p. Methylhistidine does hot have hits in metacyc
  # obtain concentration fitted data 
  metab_data <- readRDS(file = paste0(dir, glue("ResultsFiles/data/raw_{timepoint}_tar_phyloseq_obj.rds", 
                                                timepoint = toupper(timepoint)))) %>% 
    sample_data() %>% 
    as("data.frame") %>% select(-c("p.Methylhistidine")) %>% 
    rownames_to_column("ids") %>% as_tibble()
  
  annotations <- readRDS(file = glue("picrust2/annotated_path_{timepoint}.rds", timepoint = timepoint))
  
  aug_dat <- right_join(path_data, annotations, by = "pathways")
  
  name <- colnames(metab_data)[-1][met_id]
  if (name == "Propylene.glycol"){
    filt_name <- "Propylene glycol"
  } else {
    filt_name <- name
  }
  dat <- aug_dat %>% filter(!!rlang::sym(filt_name) == 1) %>% select(c(pathways, starts_with("samp")))
  if (nrow(dat) < 1){
    result <- tibble(timepoint = timepoint, met = name, path = NA, est = NA, lower = NA, upper = NA, pval = NA)
  } else {
    result <- map_df(seq_len(nrow(dat)), ~{
      path_abun <- dat %>% dplyr::slice(.x) %>% select(-pathways) %>% unlist(., use.names = FALSE)
      if (any(is.na(path_abun))){
        tibble(timepoint = timepoint, met = name, path = NA, est = NA, lower = NA, upper = NA, pval = NA)
      } else {
        met_abun <- metab_data %>% pull(name)
        test <- stats::cor.test(path_abun, met_abun, method = "spearman", conf.level = 0.95)
        rnk_path <- rank(path_abun)
        rnk_met <- rank(met_abun)
        p_ci <- cor.test(rnk_path, rnk_met, conf.level = 0.95, method = "pearson")
        tibble(timepoint = timepoint, met = name, path = dat$pathways[.x], 
               est = test$estimate, lower = p_ci$conf.int[1], upper = p_ci$conf.int[2], pval = test$p.value)
      }
    })
  }
  return(result)
}


n_metab <- 35
grid <- cross(list(
  met_id = seq_len(n_metab), 
  timepoint = c("6w", "12m")
))


sig_test <- map_df(grid, ~{
  print(.x$timepoint)
  print(.x$met_id)
  get_correlation(met_id = .x$met_id, timepoint = .x$timepoint)
})

descript_12m <- fread(file = "picrust2/path_results/12m_descript.tsv.gz") %>% as_tibble()
descript_6w <- fread(file = "picrust2/path_results/6w_descript.tsv.gz") %>% as_tibble()

descript <- bind_rows(descript_12m, descript_6w) %>% group_by(pathway, description) %>% 
  select(pathway, description)
colnames(descript)[1] <- "path"

descript <- descript %>% distinct()

sig_test <- left_join(sig_test, descript, by = "path") 

sig_test <- sig_test %>% mutate(qval = p.adjust(pval, method = "BH"))
sig_test <- sig_test %>% select(timepoint, met, path, description, est, lower, upper, pval, qval)

saveRDS(sig_test, file = "picrust2/sig_test.rds")


# Getting entire correlation heatmap
timepoint <- "12m"

correlation_analysis <- function(timepoint){
  path <- readRDS(glue("picrust2/path_results/path_abun_{timepoint}_proc.rds", timepoint = timepoint))
  met <- readRDS(file = paste0(dir, glue("ResultsFiles/data/raw_{timepoint}_tar_phyloseq_obj.rds", 
                                         timepoint = toupper(timepoint)))) %>% sample_data()
  
  path <- path %>% as.data.frame() 
  rownames(path) <- path[,1]
  path <- path[,-1]
  path <- t(path)
  path <- as.data.frame(path)
  
  correlation <- cor(x = path, y = met, method = "spearman")
  correlation <- rcorr(x = as.matrix(path), y = as.matrix(met), "spearman")
  r <- correlation$r[1:ncol(path),-c(1:ncol(path))]
  p <- correlation$P[1:ncol(path),-c(1:ncol(path))]
  adj_p <- p.adjust(p, method = "BH") %>% matrix(nrow = ncol(path), ncol = ncol(met))
  rownames(adj_p) <- rownames(p)
  colnames(adj_p) <- colnames(p)
  return(list(r = r, p_raw = p, p_adj = adj_p))
}

cor_12m <- correlation_analysis("12m")
cor_6w <- correlation_analysis("6w")

saveRDS(cor_12m, file = "picrust2/correlation_12m.rds")
saveRDS(cor_6w, file = "picrust2/correlation_6w.rds")
