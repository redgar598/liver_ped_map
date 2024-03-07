#!/bin/bash
#SBATCH --mem=200G
#SBATCH -J test
#SBATCH -p veryhimem
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH --output=%x.out


module load R/4.2.1




Rscript R_functions/test.R
##Rscript R_functions/test_notarget.R

#Rscript scripts/02_IFALD_ped_adult_fetal_integration.R
