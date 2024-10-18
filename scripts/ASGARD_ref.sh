#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=3:00:00
#SBATCH -J ASGARD_ref
#SBATCH -c 8
#SBATCH --account=def-gbader
#SBATCH --output=%x.out


module load r/4.4.0

Rscript scripts/ASGARD_Ref_compute_can.R