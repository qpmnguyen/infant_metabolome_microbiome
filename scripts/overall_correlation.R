library(PMA)
library(MCTools)
library(reshape2)
library(gplots)
library(RColorBrewer)
library(rafalib)

data <- load_main(file = "./data/data_directory.csv")
source("./scripts/utils.R")
key <- key_processing(data$tax.key, type = "NA")

generate_postCCA_correlation_plots <- function(time, method, level){
  filename = paste0("./data/",time,"_",method, "NMR_clr_tax.rds")
  data <- readRDS(file = filename)
  met <- data$met
  tax <- data$tax
  
  cv <- CCA.permute(x = tax, z = met, typex = "standard", typez = "standard", nperms = 50, niter = 5, standardize = T)
  cca <- CCA(x = tax, z = met, K = 1, penaltyx = cv$bestpenaltyx, penaltyz = cv$bestpenaltyz, niter = 50)
  tax_idx <- which(cca$u != 0)
  met_idx <- which(cca$v != 0)
  combined <- cbind(tax[,tax_idx], met[,met_idx])
  
  
  correlation <- cor(tax[,tax_idx], met[,met_idx], method = "spearman")
  rownames(correlation) <- sapply(rownames(correlation), function(x){
    idx <- which(rownames(key) == x)
    return(key[idx,level])
  })
  if (length(met_idx) == 1){
    colnames(correlation) <- colnames(met)[met_idx]
  }
  dir.1 <- ifelse(cca$u[tax_idx] > 0,"Pos","Neg")
  dir.2 <- ifelse(cca$v[met_idx] > 0, "Pos", "Neg")
  cols.1 <- palette(c('coral2','steelblue2'))[as.fumeric(dir.1)]
  cols.2 <- palette(c('coral2','steelblue2'))[as.fumeric(dir.2)]
  name <- paste0(time, "_", method, "_CCA_heatmap.png")
  png(filename = name, width = 20, height = 15, units = "in", res = 300)
  heatmap.2(correlation, col=brewer.pal(11,"RdBu"),
            RowSideColors = cols.1,
            ColSideColors = cols.2, trace = "none", key = FALSE)
  dev.off()
}

grid <- expand.grid(time <- c("12M", "6W"), method <- c("tar", "untar"))

for (i in 1:nrow(grid)){
  generate_postCCA_correlation_plots(time = grid$Var1[i], method = grid$Var2[i], level = "Genus")
}


