#!/bin/bash
#SBATCH --mem=150G
#SBATCH -J fetal_liver
#SBATCH -p veryhimem
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH --output=%x.out


module load R/4.2.1

Rscript scripts/fetal_liver.R

