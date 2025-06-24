#!/bin/bash
#SBATCH --mem=150G
#SBATCH -J BA_integration
#SBATCH -p veryhimem
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH -t 5-00:00:00
#SBATCH --output=%x.out


module load R/4.2.1

Rscript scripts/BA_pedmap_integration.R



