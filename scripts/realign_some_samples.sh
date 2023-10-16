#!/bin/bash
#SBATCH --mem=150G
#SBATCH -J realign_from_bam
#SBATCH -p veryhimem
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH --output=%x.out


module load cellranger/3.1.0

cellranger bamtofastq --nthreads=8 /cluster/projects/macparland/RE/PediatricAdult/realign_samples/possorted_genome_bam.bam  /cluster/projects/macparland/RE/PediatricAdult/realign_samples


#cellranger count --id=IFALD006_realign \
#   --fastqs=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/... \
#   --sample=IFALD006_realign \
#   --transcriptome=/cluster/tools/software/centos7/cellranger/3.1.0/cellranger-tiny-ref/3.0.0 



