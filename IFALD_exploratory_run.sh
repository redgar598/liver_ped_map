#!/bin/bash
#SBATCH --mem=150G
#SBATCH -J QC_adult_ped_IFALD
#SBATCH -p veryhimem
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH --output=%x.out


module load R/4.2.1




##Rscript R_functions/test.R
##Rscript R_functions/test_notarget.R
##Rscript scripts/01_ped_liver_QC_with_dropletQC_IFALD.R

Rscript scripts/01_ped_liver_QC_with_dropletQC_SCINA_IFALD.R