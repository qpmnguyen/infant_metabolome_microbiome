library(tidyverse)
library(data.table)
library(glue)

# Main directory 
key_dir <- "//dartfs-hpc/rc/Lab/H/HoenA/Lab/QNguyen/"

# Key file for metabolomics
key_metab <- readRDS(file = paste0(key_dir, "SourceFiles/NMR/nmr_concentrations_ST_08May2018_key.rds"))

# Metabolomics annotation files 
annotations <- read.csv(file = "picrust2/metacyc_queries/metacyc_smarttable_path_-_2021-02-10T11_23_43-08_00.csv")

# Saving metabolite names for queries 
key_metab %>% as.data.frame() %>% filter(!is.na(`KEGG Compound ID`)) %>% 
  select(c(`Metabolite`, `KEGG Compound ID`)) %>% mutate(Metabolite = tolower(Metabolite)) %>%
  unite("input", everything(), sep="$") %>% write.table(file = "picrust2/metacyc_queries/metabolite_names.tsv", col.names = FALSE,
                                                        sep = "\t", row.names = FALSE)
# Last metabolite does not have any hits
metab <- key_metab %>% as.data.frame() %>% filter(!is.na(`KEGG Compound ID`)) %>% 
  select(Metabolite) %>% slice(1:35) 

# Processing metacyc queries 
process_pathways <- function(input){
  return(strsplit(input, "//") %>% sapply(str_trim))
}

# combining metabolite keys with query names 
metab <- dplyr::bind_cols(metab, annotations)
metab <- metab %>% mutate(path_consume = map(Pathways.that.consume.compound, process_pathways)) %>% 
  mutate(path_produce = map(Pathways.that.produce.compound, process_pathways)) %>% 
  select(-c(Pathways.that.produce.compound, Pathways.that.consume.compound)) %>% 
  as_tibble()

# metacyc 
colnames(metab)[1:2] <- c("compounds", "metacyc_ids")
saveRDS(metab, file = "picrust2/annotated_metablist.rds")

# Saving pathway list to upload to metacyc
for (i in c("6w", "12m")){
  data <- fread(glue("picrust2/path_results/path_abun_{time}.tsv.gz", time = i))
  data$pathway %>% write.table(file = glue("picrust2/metacyc_queries/pathlist_{time}.tsv", time = i),
                               sep = "\t", row.names = FALSE, col.names = FALSE)
}


# Processing pathway abundance data 
process_path_abun <- function(timepoint){
  metab <- readRDS(file = "picrust2/annotated_metablist.rds")
  path <- read.csv(file = glue("picrust2/metacyc_queries/pathways_{timepoint}.csv", timepoint = timepoint))
  path <- path %>% mutate(prod = map(path_prod, process_pathways), rec = map(path_rec, process_pathways)) %>%
    select(pathways, prod, rec)
  for (i in seq_len(nrow(metab))){
    name <- metab$compounds[i]
    query <- metab$metacyc_ids[i]
    path <- path %>% mutate(!!name := map_dbl(prod, ~ifelse(query %in% .x, 1,0)))
  }
  return(path)
}


path_6w <- process_path_abun("6w")
path_12m <- process_path_abun("12m")
saveRDS(path_6w, file = "picrust2/annotated_path_6w.rds")
saveRDS(path_12m, file = "picrust2/annotated_path_12m.rds")


data <- fread("picrust2/path_results/path_abun_12m.tsv.gz")
data <- data[,lapply(.SD, function(x){
  x[x == 0] <- 1
  log(x/geometricmean(x))
}), by = pathway]


# Processing abundance tables  
# process_abun_tab <- function(tab){
#   if (!"data.table" %in% class(tab)){
#     rlang::abort("Require data.table format")
#   }
#   tab <- tab[,lapply(.SD, function(x){
#     x[x == 0] <- 1 #pseudocount
#     log(x/geometricmean(x)) # centered log-ratio
#   }),by = pathway]
#   return(tab)
# }

process_abun_tab <- function(tab){
  tab <- as.data.frame(tab)
  rownames(tab) <- tab[,1]
  tab <- tab[,-1] %>% t()
  tab[tab == 0] <- 1
  tab <- acomp(tab) %>% clr() %>% unclass() 
  tab <- t(tab) %>% as.data.frame() %>% rownames_to_column("pathways")
  return(tab)
  
}

tab_6w <- fread("picrust2/path_results/path_abun_6w.tsv.gz") %>% 
  process_abun_tab()

tab_12m <- fread("picrust2/path_results/path_abun_12m.tsv.gz") %>% 
  process_abun_tab()

saveRDS(tab_12m, file = "picrust2/path_results/path_abun_12m_proc.rds")
saveRDS(tab_6w, file = "picrust2/path_results/path_abun_6w_proc.rds")
