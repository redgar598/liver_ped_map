#!/bin/bash
#SBATCH -t 1:00 
#SBATCH --mem=256M
#SBATCH -J test 
#SBATCH -p all 
#SBATCH -c 1 
#SBATCH -N 1 
#SBATCH --output=%x.out


module load R/4.2.1

Rscript tar_make()