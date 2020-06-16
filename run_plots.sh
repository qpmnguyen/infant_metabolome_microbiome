#!/bin/bash
for i in 6W 12M; do
	for j in tar untar; do
		Rscript R/plots_correlation.R --scca output/analyses/correlation/${i}_${j}_scca.rds \
			--corr output/analyses/correlation/${i}_${j}_spearman.rds
	done
done
