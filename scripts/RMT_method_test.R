dat <- readRDS(file = "./data/12M_untarNMR_clr_tax.rds")
tax <- dat$tax
met <- dat$met

#met <- asin(sqrt(met))
met <- log(met + 1)


sample <- prcomp(tax, center = F, scale = F)

eigen <- (sample$sdev)^2
max(eigen)

sample2 <- prcomp(met, center = F, scale = F)
eigen2 <- (sample2$sdev)
max(eigen2)

covariance <- cov(t(tax))
max(eigen(covariance)$value)

covariance2 <- cov(t(met))
max(eigen(covariance2)$value)


