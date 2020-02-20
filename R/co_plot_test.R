library(ggpubr)
library(tidyverse)
data2 <- readRDS(file = "output/analyses/ordinations/12M_tar_distance.rds")
data1 <- readRDS(file = "output/analyses/ordinations/6W_tar_distance.rds")
tax <- data1$tax_euclidean
met <- data1$met_euclidean
t1 <- cbind(tax %>% unclass() %>% as.data.frame() %>% pivot_longer(everything()) %>% mutate(tax = value) %>% select(-value,-name), 
      met %>% unclass() %>% as.data.frame() %>% pivot_longer(everything()) %>% mutate(met = value) %>% select(-value, -name))
tax <- data2$tax_euclidean
met <- data2$met_euclidean
t2 <- cbind(tax %>% unclass() %>% as.data.frame() %>% pivot_longer(everything()) %>% mutate(tax = value) %>% select(-value,-name), 
            met %>% unclass() %>% as.data.frame() %>% pivot_longer(everything()) %>% mutate(met = value) %>% select(-value, -name))

t_full <- rbind(t1 %>% mutate(time = rep("6 Weeks",nrow(t1))),
                t2 %>% mutate(time = rep("12 months",nrow(t2))))

ggplot(t_full, aes(x = tax, y = met)) + geom_point() + coord_fixed(ylim = c(0, max(tax)), xlim = c(0,max(t$tax))) + 
  theme_pubr() + geom_bin2d() + scale_fill_viridis_c() + geom_hline(yintercept = median(met)) + geom_vline(xintercept = median(tax))
ggplot(t_full, aes(met, fill = time)) + geom_density()
ggplot(t_full, aes(tax, fill = time)) + geom_density()
