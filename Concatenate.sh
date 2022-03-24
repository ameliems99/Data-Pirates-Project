#!/bin/bash
#SBATCH --partition=standard
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:20:00

Rscript FINAL_CODE.R
