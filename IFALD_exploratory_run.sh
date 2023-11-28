#!/bin/bash
#SBATCH --mem=150G
#SBATCH -J QC_adult_ped_IFALD
#SBATCH -p veryhimem
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH --output=%x.out


module load R/4.2.1



#Rscript scripts/test.R
##Rscript R_functions/test_notarget.R

Rscript scripts/01_ped_liver_QC_with_dropletQC_SCINA_IFALD.R

#Rscript scripts/01_ped_liver_QC_with_dropletQC_SCINA_IFALD_PBMC_only.R

#Rscript scripts/03_PBMC_correction_IFALD.R