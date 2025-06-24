## fxiing C70

salloc -c 1 -t 1:0:0 --mem 1G
cd /cluster/projects/macparland/RE/PediatricAdult/realign_samples

scp /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/MacParland_Sonya__C70_Caudate_5pr/outs/possorted_genome_bam.bam t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/realign_samples



cellranger bamtofastq --nthreads=8 /cluster/projects/macparland/RE/PediatricAdult/realign_samples/possorted_genome_bam.bam  /cluster/projects/macparland/RE/PediatricAdult/realign_samples/C70_redo


cellranger count --id=IFALD006_realign \
   --fastqs=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/IFALD006/MacParland_Sonya__HSC-FI_006_0_1_HMW73DMXX \
   --sample=bamtofastq \
   --chemistry=SC3Pv3 \
   --transcriptome=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/reference/GRCh38



40,729,361 reads, whereas the final matrices contain 607,032,905 reads. 

zgrep -c "^@" C70/MacParland_Sonya__C70_Caudate_5pr_0_1_HJHT3DRXX/bamtofastq_S1_L001_R2_005.fastq.gz

zcat C70/MacParland_Sonya__C70_Caudate_5pr_0_1_HJHT3DRXX/bamtofastq_S1_L001_R2_005.fastq.gz | wc -l


## after upload lets try this again 
sbatch realign_some_samples.sh



salloc -c 1 -t 2:0:0 --mem 10G
cd /cluster/projects/macparland/RE/PediatricAdult/realign_samples
module load cellranger/3.1.0



cellranger bamtofastq --nthreads=8 /cluster/projects/macparland/RE/PediatricAdult/realign_samples/possorted_genome_bam.bam  /cluster/projects/macparland/RE/PediatricAdult/realign_samples/C70




cellranger count --id=C70_realign \
   --fastqs="/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C70/MacParland_Sonya__C70_Caudate_5pr_0_1_HMW73DMXX,/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C70/MacParland_Sonya__C70_Caudate_5pr_0_1_HJHT3DRXX" \
   --sample=bamtofastq \
   --chemistry=SC5P-R2 \
   --transcriptome=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/reference/GRCh38




### from SRA
cd /cluster/projects/macparland/RE/PediatricAdult/realign_samples

mkdir C70_SRA_test

scp /media/redgar/Seagate\ Portable\ Drive/C70_fastq/SRR26173958.fastq.gz t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C70_SRA_test

cellranger count --id=C70_SRA_test_ranger \
   --fastqs="/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C70_SRA_test" \
   --sample=SRR26173958 \
   --chemistry=SC5P-R2 \
   --transcriptome=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/reference/GRCh38

export PATH=$PWD/sratoolkit.3.2.0-ubuntu64/bin:$PATH



## on external
export PATH=$PWD/sratoolkit.3.2.0-ubuntu64/bin:$PATH
fastq-dump --split-files --gzip SRR26173958




# if needed try veryhimem




##########################
## with extra fastqs sent
##########################

salloc -c 1 -t 1:0:0 --mem 1G
cd /cluster/projects/macparland/RE/PediatricAdult/realign_samples
mkdir C70_morefastq


salloc -c 1 -t 1:0:0 --mem 1G
cd /cluster/projects/macparland/RE/PediatricAdult/realign_samples/C70_morefastq

scp /media/redgar/Seagate\ Portable\ Drive/synology_sonya/SynologyDrive/fastq-C70/C70_Caudate_5pr_*.fastq.gz t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C70_morefastq
scp /media/redgar/Seagate\ Portable\ Drive/C70_more_fastq/fastq/C70_Caudate_5pr_10XBarcode_*.fastq.gz t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C70_morefastq

# update sample names of files
for file in C70_Caudate_5pr_10XBarcode*.fastq.gz; do
    mv "$file" "${file/C70_Caudate_5pr_10XBarcode/C70_Caudate_5pr}"
done

salloc -c 2 -t 4:0:0 --mem 10G
cd /cluster/projects/macparland/RE/PediatricAdult/realign_samples/C70_morefastq

module load cellranger/3.1.0

cellranger count --id=C70_more_fastq \
   --fastqs="/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C70_morefastq" \
   --sample=C70_Caudate_5pr \
   --chemistry=SC5P-R2 \
   --transcriptome=/cluster/projects/macparland/RE/PediatricAdult/realign_samples/reference/GRCh38


nano /cluster/projects/macparland/RE/PediatricAdult/realign_samples/C70_morefastq/realign_from_bam_C70.out

C70_Caudate_5pr_S19_L001_I1_001.fastq.gz  C70_Caudate_5pr_S19_L002_R2_001.fastq.gz  C70_Caudate_5pr_S3_L002_R1_001.fastq.gz  C70_Caudate_5pr_S5_L002_I1_001.fastq.gz
C70_Caudate_5pr_S19_L001_R1_001.fastq.gz  C70_Caudate_5pr_S3_L001_I1_001.fastq.gz   C70_Caudate_5pr_S3_L002_R2_001.fastq.gz  C70_Caudate_5pr_S5_L002_R1_001.fastq.gz
C70_Caudate_5pr_S19_L001_R2_001.fastq.gz  C70_Caudate_5pr_S3_L001_R1_001.fastq.gz   C70_Caudate_5pr_S5_L001_I1_001.fastq.gz  C70_Caudate_5pr_S5_L002_R2_001.fastq.gz
C70_Caudate_5pr_S19_L002_I1_001.fastq.gz  C70_Caudate_5pr_S3_L001_R2_001.fastq.gz   C70_Caudate_5pr_S5_L001_R1_001.fastq.gz 
C70_Caudate_5pr_S19_L002_R1_001.fastq.gz  C70_Caudate_5pr_S3_L002_I1_001.fastq.gz   C70_Caudate_5pr_S5_L001_R2_001.fastq.gz



scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C70_morefastq/C70_more_fastq/outs/*html /home/redgar/Downloads
##40,729,361 reads, whereas the final matrices contain 607,032,905 reads. 


## somehow only S5 fatsqs were counted by Erica
## with all the fastqs the reads are closer to right