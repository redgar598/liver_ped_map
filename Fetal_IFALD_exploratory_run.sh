#!/bin/bash
#SBATCH --mem=150G
#SBATCH -J QC_Fetal_adult_ped_IFALD
#SBATCH -p veryhimem
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH --output=%x.out


module load R/4.2.1




##Rscript R_functions/test.R
##Rscript R_functions/test_notarget.R

Rscript scripts/IFALD_ped_adult_fetal_integration.R
