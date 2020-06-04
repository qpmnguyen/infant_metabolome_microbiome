#!/usr/bin/env bash
RENV_PATHS_CACHE_HOST=/home/quangnguyen/.local/renv/cache
RENV_PATHS_CACHE_CONTAINER=/renv/cache
DATA_HOST=/mnt/HoenLab/Lab/QNguyen/ResultsFiles/data
DATA_CONTAINER=/analysis/data
SCRIPT_HOST=/home/quangnguyen/research/infant_metabolome_robust_variation_microbiome/R
SCRIPT_CONTAINER=/analysis/R
RESULTS_HOST=/home/quangnguyen/research/results/infant_metabolome
RESULTS_CONTAINER=/analysis/results
docker run --rm \
	-e PASSWORD=password \
	-e "RENV_PATHS_CACHE=${RENV_PATHS_CACHE_CONTAINER}"\
	-v "${RENV_PATHS_CACHE_HOST}:${RENV_PATHS_CACHE_CONTAINER}" \
	-v "${DATA_HOST}:${DATA_CONTAINER}" \
	-v "${SCRIPT_HOST}:${SCRIPT_CONTAINER}" \
	-v "${RESULTS_HOST}:${RESULTS_CONTAINER}" \
	-w "/analysis" \
	-p 8787:8787 rocker/rstudio:3.6.0


