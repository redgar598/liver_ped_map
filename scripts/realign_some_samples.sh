#!/bin/bash
#SBATCH --mem=50G
#SBATCH -J realign_compare
#SBATCH -p himem
#SBATCH --time=20:00:00
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





#cellranger bamtofastq --nthreads=8 /cluster/projects/macparland/RE/PediatricAdult/realign_samples/possorted_genome_bam.bam  /cluster/projects/macparland/RE/PediatricAdult/realign_samples/C64


# cellranger count --id=C64_realign \
#     --fastqs=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C64/McGilvray_Sonya__C64_Enriched_5pr_0_1_CDPFDANXX \
#     --sample=bamtofastq \
#     --chemistry=SC5P-R2 \
#     --transcriptome=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/reference/GRCh38



#cellranger bamtofastq --nthreads=8 /cluster/projects/macparland/RE/PediatricAdult/realign_samples/possorted_genome_bam.bam  /cluster/projects/macparland/RE/PediatricAdult/realign_samples/C39_NPC


# cellranger count --id=C39_NPC_realign \
#     --fastqs=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C39_NPC/McGilvery_Sonya__NPC_MissingLibrary_1_CB8R9ANXX \
#     --sample=bamtofastq \
#     --chemistry=SC3Pv2 \
#     --transcriptome=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/reference/GRCh38



# cellranger bamtofastq --nthreads=8 /cluster/projects/macparland/RE/PediatricAdult/realign_samples/possorted_genome_bam.bam  /cluster/projects/macparland/RE/PediatricAdult/realign_samples/C39_TLH

# cellranger count --id=C39_TLH_realign \
#     --fastqs=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C39_TLH/McGilvery_Sonya__TLH_MissingLibrary_1_CB8R9ANXX \
#     --sample=bamtofastq \
#     --chemistry=SC3Pv2 \
#     --transcriptome=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/reference/GRCh38


# cellranger bamtofastq --nthreads=8 /cluster/projects/macparland/RE/PediatricAdult/realign_samples/possorted_genome_bam.bam  /cluster/projects/macparland/RE/PediatricAdult/realign_samples/C54

# cellranger count --id=C54_realign \
#     --fastqs=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C54/McGilvery_Sonya__TLH_april_12_MissingLibrary_1_CCJ7PANXX \
#     --sample=bamtofastq \
#     --chemistry=SC3Pv2 \
#     --transcriptome=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/reference/GRCh38


# module load R/4.2.1

# Rscript scripts/compare_after_realign.R


##########
## bamtofastq for all samples
##########

##########
## remove identifying sample names
##########


# cd /cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/MacParland_Diana__C93_Frozen_Liver_220919_3pr_V3_1/outs



# MacParland_Diana__C93_Frozen_Liver_220919_3pr_V3_1  
# MacParland_Diana__C96_Frozen_Liver_220919_3pr_V3_1
# MacParland_Diana__SingleCell_C82_16Sept21    
# MacParland_Sonya__SingleCell_iFALD073_PBMC_25Jan21_3pr_v3
# MacParland_Sonya__SingleCell_iFALD073_Biopsy_25Jan21_3pr_v3
# X210026__SingleCell_C88_14Oct21_3pr_V3_1


## WILL NEED TO REALIGN THESE SO THE HTML HAS THE NEW SAMPLE NAME

#cellranger bamtofastq --nthreads=8 /cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/MacParland_Diana__C93_Frozen_Liver_220919_3pr_V3_1/outs/possorted_genome_bam.bam  /cluster/projects/macparland/RE/PediatricAdult/realign_samples/C93

cellranger count --id=C93_realign \
    --fastqs=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C93/MacParland_Diana__C93_Frozen_Liver_3pr_V3_1 \
    --sample=bamtofastq \
    --chemistry=SC3Pv3 \
    --transcriptome=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/reference/GRCh38

           
