#!/bin/bash
#SBATCH --mem=100G
#SBATCH -J BIDCell_C85
#SBATCH -p veryhimem
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH --output=%x.out



python scripts/liver_BIDCell_C85.py