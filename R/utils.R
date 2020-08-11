# Script to hold auxiliary functions 


# get directories using time, metab status, and data stage
get_dir <- function(time, mettype, data_stage, sens){
  dir <- "//dartfs-hpc/rc/lab/H/HoenA/Lab/QNguyen/ResultsFiles/data"
  if (sens == F){ # if not in the pre-defined sensitivity analysis mode 
    file_path <- glue("{dir}/{data_stage}_{time}_{mettype}_phyloseq_obj.rds",
                      dir = dir, data_stage = data_stage, time = time, mettype = mettype)  
  } else if (sens == T){
    file_path <- glue("{dir}/{data_stage}_common_{time}_{mettype}_phyloseq_obj.rds",
                      dir = dir, data_stage = data_stage, time = time, mettype = mettype)
  }
  return(file_path)
  
}

# supplemental functions 
#' @title Testing for significance of principal components
#' @description Testing for significance of eigenvalues using Tracy-Widom test statistics. Details
#'    see Frost et al. 2015
#' @param prcomp.output Principal component analysis object.  
rmt.test <- function(prcomp.output){
  library(RMTstat) # requiring RMTstat package 
  eigenvals <- (prcomp.output$sdev)^2
  n.eigenvec <- nrow(prcomp.output$rotation)
  sample.size <- nrow(prcomp.output$x)
  tws <- c()
  pvals <- c()
  for (i in 1:length(eigenvals)){
    eval <- eigenvals[i] #selecting the eigenvalue
    pdim <- n.eigenvec - i + 1 #number of dimensions of the wishart matrix
    par <- WishartMaxPar(ndf = sample.size, pdim = pdim)
    tws[i] <- (eval - par$centering)/par$scaling
    pvals[i] <- pgamma(tws[i] + 9.84801, shape=46.446, scale=0.186054, lower.tail = F)
  }
  result = list(tws = tws, pvals = pvals)
  return(result)
}