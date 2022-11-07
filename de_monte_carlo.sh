#!/bin/bash
#SBATCH -t 5-00:00:00
#SBATCH --mem=20G
#SBATCH -J monte_carlo_DE
#SBATCH -p all
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH --output=%x.out


cellType_index=$1

module load R/4.2.1


Rscript R_functions/ped_liver_differential_expression_montecarlo.R $cellType_index









