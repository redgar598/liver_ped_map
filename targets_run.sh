#!/bin/bash
#SBATCH --mem=100G
#SBATCH -J targets
#SBATCH -p veryhimem 
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH --output=%x.out


module load R/4.2.1




Rscript R_functions/test.R
##Rscript R_functions/test_notarget.R
##Rscript R_functions/ped_liver_QC.R

