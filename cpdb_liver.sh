#!/bin/bash
#SBATCH --mem=150G
#SBATCH -J cpdb_ped_healthy
#SBATCH -p veryhimem
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH --output=%x.out



### generate data for cpdb
module load R/4.2.1

Rscript scripts/cellphonedb_datasetup.R


conda activate cpdb 
python -u psorasis_import.py