#!/usr/bin/env bash 
Rscript R/spearman_corr.R --input /dartfs-hpc/rc/lab/H/HoenA/Lab/QNguyen/ResultsFiles/data/processed_12M_tar_prediction_phyloseq_obj.rds \
	--output output/analyses/correlation/12M_tar_spearman.rds \
	--metric spearman --MHC BH
