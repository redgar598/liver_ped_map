#!/bin/bash
#SBATCH --mem=50G
#SBATCH -J fetal_liver
#SBATCH -p himem
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH --output=%x.out


module load R/4.2.1

Rscript scripts/fetal_liver.R

