#!/bin/bash
#SBATCH --mem=150G
#SBATCH -J realign_from_bam
#SBATCH -p long
#SBATCH --time=10:00:00
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH --output=%x.out


module load cellranger/3.1.0

#cellranger bamtofastq --nthreads=8 /cluster/projects/macparland/RE/PediatricAdult/realign_samples/possorted_genome_bam.bam  /cluster/projects/macparland/RE/PediatricAdult/realign_samples/IFALD006


cellranger count --id=IFALD006_realign \
   --fastqs=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/IFALD006/MacParland_Sonya__HSC-FI_006_0_1_HMW73DMXX \
   --sample=bamtofastq \
   --chemistry=SC3Pv3 \
   --transcriptome=/cluster/tools/software/centos7/cellranger/3.1.0/cellranger-tiny-ref/3.0.0 


# mv output to projects
#mv /cluster/home/t117652uhn/liver_ped_map/IFALD006_realign /cluster/projects/macparland/RE/PediatricAdult/realign_samples

