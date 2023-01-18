#!/bin/bash
#SBATCH --mem=5G
#SBATCH -J droplet_QC_test
#SBATCH -p all 
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH --output=%x.out


module load R/4.2.1

Rscript scripts/test.R

