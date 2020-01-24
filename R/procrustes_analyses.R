library(vegan)
library(ade4)
library(optparse)

option_list <- list(
  make_option("--input", help = "Input file for data loading"),
  make_option("--output", help = "Output file for data loading")
)

opt <- parse_args(OptionParser(option_list = option_list))

# accessing dataset 
data <- readRDS(file = opt$input)
grid <- expand.grid(
  tax_dist = names(data)[1:2],
  met_dist = names(data)[3:4]
)
k <- dim(as.matrix(data[[1]]))[1] - 1
output <- list()
# Performing ordinations  
for (i in 1:nrow(grid)){
  tax_ord <- cmdscale(d = data[[as.character(grid$tax_dist[i])]], eig = T, k = 10)
  met_ord <- cmdscale(d = data[[as.character(grid$met_dist[i])]], eig = T, k = 10)
  proc_test <- protest(tax_ord, met_ord)
  results <- list(tax_ord = tax_ord, met_ord = met_ord, proc_test = proc_test)
  output[[i]] <- results
}

# saving results  
names(output) <- with(grid, paste(tax_dist, met_dist, sep = "_"))
saveRDS(output, file = opt$output)
