#!/bin/bash
#SBATCH --mem=150G
#SBATCH -t 5-00:00:00
#SBATCH -J MTM1
#SBATCH -p veryhimem
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH --output=%x.out


module load R/4.2.1




##Rscript R_functions/test.R
##Rscript R_functions/test_notarget.R

Rscript scripts/human_fetal_liver.R
