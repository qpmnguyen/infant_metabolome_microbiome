#!/usr/bin/env bash 
module load python/3.7-Anaconda
source activate snakemake
snakemake -s snakefiles/correlation.rule.py
