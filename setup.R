# using rstudio api to automatically initialize an R project 

# installing renv 
if (!requireNamespace("renv")){
  if (!requireNamespace("devtools")){
    install.packages('devtools')
  } else {
    devtools::install_github('rstudio/renv')
  }
}

# installing BiocManager
if (!requireNamespace('BiocManager')){
  install.packages('BiocManager')
}

# making sure that bioconductor repositories are added to options and union with previous versions
options(repos = c(getOption("repos"), BiocManager::repositories()))
renv::status() # check status of renv
renv::restore()
renv::snapshot() # snapshot after restoring