#!/bin/bash
#SBATCH --mem=50G
#SBATCH -J realign_from_bam
#SBATCH -p himem
#SBATCH --time=10:00:00
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH --output=%x.out


module load cellranger/3.1.0







###########
## Build reference
###########

#wget ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

#wget ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz
#gunzip Homo_sapiens.GRCh38.93.gtf.gz

# cellranger mkgtf Homo_sapiens.GRCh38.93.gtf Homo_sapiens.GRCh38.93.filtered.gtf \
#                 --attribute=gene_biotype:protein_coding \
#                 --attribute=gene_biotype:lincRNA \
#                 --attribute=gene_biotype:antisense \
#                 --attribute=gene_biotype:IG_LV_gene \
#                 --attribute=gene_biotype:IG_V_gene \
#                 --attribute=gene_biotype:IG_V_pseudogene \
#                 --attribute=gene_biotype:IG_D_gene \
#                 --attribute=gene_biotype:IG_J_gene \
#                 --attribute=gene_biotype:IG_J_pseudogene \
#                 --attribute=gene_biotype:IG_C_gene \
#                 --attribute=gene_biotype:IG_C_pseudogene \
#                 --attribute=gene_biotype:TR_V_gene \
#                 --attribute=gene_biotype:TR_V_pseudogene \
#                 --attribute=gene_biotype:TR_D_gene \
#                 --attribute=gene_biotype:TR_J_gene \
#                 --attribute=gene_biotype:TR_J_pseudogene \
#                 --attribute=gene_biotype:TR_C_gene

# cd /cluster/projects/macparland/RE/PediatricAdult/realign_samples/reference
# cellranger mkref --genome=GRCh38 \
#                 --fasta=Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#                 --genes=Homo_sapiens.GRCh38.93.filtered.gtf \
#                 --ref-version=3.0.0





# cellranger bamtofastq --nthreads=8 /cluster/projects/macparland/RE/PediatricAdult/realign_samples/possorted_genome_bam.bam  /cluster/projects/macparland/RE/PediatricAdult/realign_samples/IFALD006


# cellranger count --id=IFALD006_realign \
#    --fastqs=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/IFALD006/MacParland_Sonya__HSC-FI_006_0_1_HMW73DMXX \
#    --sample=bamtofastq \
#    --chemistry=SC3Pv3 \
#    --transcriptome=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/reference/GRCh38





cellranger bamtofastq --nthreads=8 /cluster/projects/macparland/RE/PediatricAdult/realign_samples/possorted_genome_bam.bam  /cluster/projects/macparland/RE/PediatricAdult/realign_samples/C64


# cellranger count --id=C64_realign \
#    --fastqs=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C64/MacParland_Diana__C102_5pr_V2_somthing \
#    --sample=bamtofastq \
#    --chemistry=SC5P-R2 \
#    --transcriptome=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/reference/GRCh38

