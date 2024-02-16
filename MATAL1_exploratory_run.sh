#!/bin/bash
#SBATCH --mem=150G
#SBATCH -J MATAL1_QC
#SBATCH -p veryhimem
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH -t 5-00:00:00
#SBATCH --output=%x.out


module load R/4.2.1



#Rscript scripts/test.R
##Rscript R_functions/test_notarget.R

#Rscript scripts/01_ped_liver_QC_with_dropletQC_SCINA_IFALD.R

#Rscript scripts/01_ped_liver_QC_with_dropletQC_SCINA_IFALD_PBMC_only.R

#Rscript scripts/03_PBMC_correction_IFALD.R





Rscript scripts/MALAT1_QC_Zoe.R
