#!/usr/bin/env bash 
declare -a mod=("svm" "rf" "enet" "spls")
for val in ${mod[@]}; do 
	Rscript R/prediction_eval.R --model val --metab_type tar --time 12M --tax_type 16S
done
