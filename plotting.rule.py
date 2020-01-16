TIME = ["6W", "12M"]
MET = ["tar", "untar"]
TAX = ["16S"]
METHODS = ["rf", "svm", "enet", "spls", "xgboost", "gbm"]
EXT = ['rds', 'svg']


rule prediction_plots:
  input:
  output:
  scripts:
    
