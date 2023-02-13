#!/bin/bash
#SBATCH --mem=100G
#SBATCH -J soup_X
#SBATCH -p veryhimem
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH --output=%x.out


module load R/4.2.1

Rscript scripts/soupX_explore.R

