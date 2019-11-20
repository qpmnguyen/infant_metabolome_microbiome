library(PMA)
library(vegan)
library(phyloseq)
library(optparse)


data <- readRDS(file = "data/processed/16S_12M_tar_processed_prediction.rds")
tax <- data$tax
met <- data$met
mod <- CCA.permute(x = tax, z = met, typex = "standard", typez = "standard", nperms = 25, niter = 5, standardize = T, K = 36)
cca <- CCA(x = tax, z = met, K = 36, typex = "standard", typez = "standard", niter = 5, standardize = T, penaltyx = mod$bestpenaltyx, 
           penaltyz = mod$bestpenaltyz)
plot(cca$cors^2, type = "line")
cor(cca$u[,1], cca$v[,1])
length(cca$v[,1])
length(cca$u[,1])
mean(cca$v[,1]^2)*cca$cors[1]
red <- c()
for (i in 1:36){
  red[i] <- mean(cca$v[,i]^2)*cca$cors[i]
}
plot(red)

plot(cca$cors^2/sum(cca$cors^2))
