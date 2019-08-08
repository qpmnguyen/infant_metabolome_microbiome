library(picante)
library(MLmetrics)
library(caret)
library(compositions)
library(doParallel)
library(foreach)
library(optparse)

option_list <- list(
  make_option("--ncores", type = "integer", help = "Number of cores")
)

opt <- parse_args(OptionParser(option_list = option_list))

input <- readRDS("/dartfs-hpc/rc/home/k/f00345k/research/metadecon/datasets/real_data/raw_data.rds")

tax <- input$tax
met <- input$met

registerDoParallel(opt$ncores)
ctrl <- trainControl(method = "cv", number = 2, search = "grid", verboseIter = F)
output <- foreach (i = 1:1000) %dopar% {
  if (i %% 5 == 0 & i != 0){
    print(paste("At", (i/1000)*100, "%"))
  }
  rand_tax <- randomizeMatrix(tax, null.model = "richness", iterations = 1000)
  rand_tax <- unclass(acomp(rand_tax))
  rand_tax[rand_tax == 0] <- 1
  rand_tax <- unclass(clr(rand_tax))
  rand_tax <- scale(rand_tax) # scaling 
  met_idx <- sample(1:ncol(met), size = 1)
  rand_met <- randomizeMatrix(met[,met_idx], null.model = "richness", iterations = 1000)
  train_idx <- as.vector(caret::createDataPartition(met[,1], times = 1, p = 0.8, list = F))
  # transform and split 
  train_t <- rand_tax[train_idx,]
  train_m <- rand_met[train_idx]
  train_m <- log(train_m + 1)
  test_t <- rand_tax[-train_idx,]
  test_m <- rand_met[-train_idx]
  # run model fit procedure 
  mod <- train(x = train_t, y = train_m, trControl = ctrl, tuneLength = 10, method = 'rf')
  predictions <- predict(mod, newdata = test_t)
  predictions <- exp(predictions) - 1 # inv transform  
  # final values  
  RRSE <- MLmetrics::RRSE(y_pred = predictions, y_true = test_m)
  RMSE <-  MLmetrics::RMSE(y_pred = predictions, y_true = test_m)
  Corr <- cor(x = test_m, y = predictions, method = "spearman")
  R2 <- caret::R2(pred = predictions, obs = test_m, formula = "traditional")
  result <- data.frame(RRSE = RRSE, RMSE = RMSE, Corr = Corr, R2 = R2)
  result
}
saveRDS(output, file = "./null_permutation.rds")
