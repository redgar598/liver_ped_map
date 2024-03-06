
# Instructions: https://1drv.ms/w/s!AinqViGDqMQDiFsDGm9YthNkfpvD
# Issues or requests on the cluster, please send email to hpc_support@uhnresearch.ca

# Your group directory is /cluster/projects/macparland 
# Your home directory is /cluster/home/




# when connected UHN
ssh t117652uhn@h4huhnlogin1.uhnresearch.ca

# password: long passcode also used at EBI



## running pipeline
module load R/4.2.1


## see project directory
salloc -c 1 -t 1:0:0 --mem 1G
cd /cluster/projects/macparland/RE

module load R/4.2.1

## running code
salloc -c 1 -t 1:0:0 --mem 1G
cd liver_ped_map/

# need internet for package install
#so from login
cd liver_ped_map
module load R/4.2.1


# submit job on login node
cd liver_ped_map
git pull
sbatch IFALD_exploratory_run.sh
sbatch cpdb_liver.sh
sbatch Fetal_IFALD_exploratory_run.sh



sbatch MATAL1_exploratory_run.sh



sbatch exploratory_run.sh
sbatch de_monte_carlo_notparallel.sh
sbatch metacell.sh
sbatch Fetal_IFALD_exploratory_run.sh
sbatch test_run.sh


# Python
#conda activate liver_scRNAseq
#sbatch scripts/parse_reference.sh

######
## Big data transfer
######
ssh t117652uhn@h4huhndata1.uhnresearch.ca


##########
# monte carlo cell type - age
##########
for cellType_index in {1..14}
do

sbatch de_monte_carlo.sh $cellType_index

done

##########
# monte carlo cell type - sex
##########
for cellType_index in {1..12}
do

sbatch de_monte_carlo_sex.sh $cellType_index

done





#################
## Realign some samples
#################
salloc -c 1 -t 1:0:0 --mem 1G
cd /cluster/projects/macparland/RE/PediatricAdult/realign_samples
sbatch realign_some_samples.sh


scp /media/redgar/Seagate\ Portable\ Drive/IFALD/191218_A00827_0099_AHMW73DMXX_MacParland_Sonya/MacParland_Sonya__HSC-FI_006/possorted_genome_bam.bam t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/realign_samples
scp /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/McGilvray_Sonya__C64_Enriched_5pr/outs/possorted_genome_bam.bam t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/realign_samples
scp /media/redgar/Seagate\ Portable\ Drive/Bams\ 5\ liver\ map/BAM\ for\ 5\ liver\ map/USB\ Copy_2018-05-25_101435/170712_D00355_0169_ACB8R9ANXX_McGilvery_Sonya/McGilvery_Sonya__NPC/possorted_genome_bam.bam t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/realign_samples
scp /media/redgar/Seagate\ Portable\ Drive/Bams\ 5\ liver\ map/BAM\ for\ 5\ liver\ map/USB\ Copy_2018-05-25_101435/170712_D00355_0169_ACB8R9ANXX_McGilvery_Sonya/McGilvery_Sonya__TLH/possorted_genome_bam.bam t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/realign_samples

scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/realign_samples/possorted_genome_bam.bam /media/redgar/Seagate\ Portable\ Drive/realign_samples/CD54_from_diana


cp /cluster/home/t117652uhn/liver_ped_map/scripts/realign_some_samples.sh . 



scp -r t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/realign_samples/IFALD006_realign/outs /media/redgar/Seagate\ Portable\ Drive/realign_samples/IFALD006_realign
scp -r t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C64_realign/outs /media/redgar/Seagate\ Portable\ Drive/realign_samples/C64_realign
scp -r t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C39_NPC_realign/outs /media/redgar/Seagate\ Portable\ Drive/realign_samples/C39_NPC_realign
scp -r t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C39_TLH_realign/outs /media/redgar/Seagate\ Portable\ Drive/realign_samples/C39_TLH_realign
scp -r t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/realign_samples/C54_realign/outs /media/redgar/Seagate\ Portable\ Drive/realign_samples/C54_realign

# submit job on login node
cd liver_ped_map
git pull
sbatch scripts/realign_some_samples.sh
nano realign_compare.out

nano scripts/compare_after_realign.R


scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/figures/realign_*pdf /home/redgar/Documents/liver_ped_map/figures
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/figures/jpeg/realign_*jpeg /home/redgar/Documents/liver_ped_map/figures/jpeg


scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/figures/*MALAT1.pdf /home/redgar/Documents/liver_ped_map/figures


### downlaod fetal liver data
scp /media/redgar/Seagate\ Portable\ Drive/fetal_liver/E-MTAB-7407-unix-ftp_liver.txt  t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/fetal_liver
cat ./E-MTAB-7407-unix-ftp_liver.txt | sh

cat /media/redgar/Seagate\ Portable\ Drive/fetal_liver/scp_liver_only.txt| xargs -i scp {} t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/fetal_liver

scp /media/redgar/Seagate\ Portable\ Drive/fetal_liver/E-MTAB-7407.sdrf.txt  t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/fetal_liver


######
##SCI
######
scp /media/redgar/Seagate\ Portable\ Drive/spinalcord_tutorial_data/d10x_SCI_merged.rds t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data
scp /media/redgar/Seagate\ Portable\ Drive/spinalcord_tutorial_data/d10x_78_SCI_merged.rds t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data

scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/d10x_SCI_fullintegrated.rds /media/redgar/Seagate\ Portable\ Drive/spinalcord_tutorial_data 

#################
### Moving files
#################

### big raw data move
#ssh t117652uhn@h4huhndata1.uhnresearch.ca
scp /home/redgar/Documents/liver_ped_map/data/data_transfer_updated_jan16_2023.csv t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw
scp /home/redgar/Documents/liver_ped_map/data/data_transfer_updated_mar20_2023.csv t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw
scp /home/redgar/Documents/liver_ped_map/data/data_transfer_updated_mar20_2023_IFALD.csv t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw
scp /home/redgar/Documents/liver_ped_map/data/data_transfer_updated_may15_2023_IFALD_PBMC.csv t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw
scp /home/redgar/Documents/liver_ped_map/data/data_transfer_updated_june21_2023_IFALD_PBMC.csv t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw

scp /home/redgar/Documents/liver_ped_map/data/data_transfer_updated_Nov28_2023_IFALD_PBMC.csv t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data
scp /home/redgar/Documents/liver_ped_map/data/data_transfer_updated_feb12_2024_IFALD_PBMC.csv t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data


scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/C39_NPC_june6_2017 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/
scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/C39_TLH_june6_2017 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/
scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/C54_3prV2_12Apr18 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/
scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/X210026__SingleCell_C88_14Oct21_3pr_V3_1 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/
scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/MacParland_Jawairia__C104_Bx_5pr_V2 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/
scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/MacParland_Sonya__C68_Total_liver t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/
scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/MacParlnd_Sai__C97_3pr_V3_1 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/

scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/MacParland_Sonya__SingleCell_iFALD073_Biopsy_25Jan21_3pr_v3 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/
scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/MacParland_Catia__HSC_IF030_3pr_v3 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/
scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/MacParland_Sonya__HSC-FI_006 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/
scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/MacParland_Sonya__SingleCell_iFALD073_PBMC_25Jan21_3pr_v3 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/

scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/MacParland_Diana__C105_Frozen_5pr_V2 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/
scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/MacParland_Diana__C102_5pr_V2 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/

scp -r /media/redgar/Seagate\ Portable\ Drive/MacParland_Jawairia__C115_biopsy_5pr_V2 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/
scp -r /media/redgar/Seagate\ Portable\ Drive/MacParland_Jawairia__C113_5pr_V2 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/



scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/IFALD_d10x_adult_ped_raw_PBMC.rds /media/redgar/Seagate\ Portable\ Drive/processed_data 
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/IFALD_adult_ped_integrated_PBMC.rds /media/redgar/Seagate\ Portable\ Drive/processed_data 
scp /media/redgar/Seagate\ Portable\ Drive/processed_data/IFALD_d10x_adult_ped_raw_PBMC.rds t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data





######################
### feb 2024 update
######################
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/IFALD_adult_ped_integrated.rds /media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024 
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/IFALD_d10x_adult_ped_raw.rds /media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024 
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/IFALD_adult_ped_cellRough.rds /media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024 

scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/IFALD_adult_ped_SCINA_cell_labels.RData /media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024 

scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds /media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024 
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/IFALD_adult_ped_cellRefined_withDropletQC.rds /media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024 
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/IFALD_adult_ped_PBMC_integrated.rds /media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024 

scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/statistical_analysis_*.txt /home/redgar/Documents/liver_ped_map/data/cellphonedb

scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/Fetal_IFALD_adult_ped_integrated_myeloid_only.RData /media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/Fetal_IFALD_adult_ped_integrated_KC_only.RData /media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024


scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/figures/*pdf /home/redgar/Documents/liver_ped_map/figures
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/figures/jpeg/*jpeg /home/redgar/Documents/liver_ped_map/figures/jpeg
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/figures/png/*png /home/redgar/Documents/liver_ped_map/figures/png
######################







scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/MALAT1_adult_ped_integrated.rds /media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024 
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/MALAT1_adult_ped_cellRough.rds /media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024 
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/IFALD_QC_metrics.Rdata /media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024 
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/MALAT1_QC_metrics.Rdata /media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024 


scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/IFALD_adult_ped_PBMC_integrated.rds /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/IFALD_adult_ped_PBMC_integrated.rds /media/redgar/Seagate\ Portable\ Drive/processed_data 
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/IFALD_adult_ped_PBMC_merged.rds /media/redgar/Seagate\ Portable\ Drive/processed_data 



scp /media/redgar/Seagate\ Portable\ Drive/fetal_liver/download.h5ad t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/fetal_liver
scp /media/redgar/Seagate\ Portable\ Drive/fetal_liver/fetal_liver.h5ad t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/fetal_liver
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/fetal_integrated.rds /media/redgar/Seagate\ Portable\ Drive/fetal_liver 
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/fetal_celllabels.rds /media/redgar/Seagate\ Portable\ Drive/fetal_liver 
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/d10x_fetal_raw.rds /media/redgar/Seagate\ Portable\ Drive/fetal_liver 
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/fetal_celllabels.rds /media/redgar/Seagate\ Portable\ Drive/fetal_liver 



scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/QC_Fetal_adult_ped_IFALD.out /home/redgar/Documents/liver_ped_map 
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/Fetal_IFALD_adult_ped_integrated.rds /media/redgar/Seagate\ Portable\ Drive/fetal_liver 
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/Fetal_IFALD_adult_ped_cellRough.rds /media/redgar/Seagate\ Portable\ Drive/fetal_liver 
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/Fetal_IFALD_adult_ped_pltData.RData /media/redgar/Seagate\ Portable\ Drive/fetal_liver 

scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/Fetal_ped_IFALD_adult_PCA_myeloid.RData /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/Fetal_ped_IFALD_adult_diffexpression_myeloid.RData /media/redgar/Seagate\ Portable\ Drive/fetal_liver 


scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/Fetal_IFALD_adult_ped_integrated_*_only.RData /media/redgar/Seagate\ Portable\ Drive/processed_data
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/Fetal_IFALD_adult_ped_integrated_myeloid_only.RData /media/redgar/Seagate\ Portable\ Drive/processed_data
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/Fetal_IFALD_adult_ped_raw_myeloid_only.RData /media/redgar/Seagate\ Portable\ Drive/processed_data


scp /home/redgar/Documents/liver_ped_map/data/Liver_Markers*.csv t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data

scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/IFALD_adult_ped_SCINA_markers_withcitations_cell_labels.RData /media/redgar/Seagate\ Portable\ Drive/processed_data 
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/IFALD_adult_ped_SCINA_cell_labels.RData /media/redgar/Seagate\ Portable\ Drive/processed_data 



## general file moving moving files
scp /home/redgar/Documents/tissue_MHCI/R_functions/pretty_plots.R t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn


cp -R /cluster/projects/macparland/DN/PediatricAdult /cluster/projects/macparland/RE


cp /cluster/projects/macparland/RE/PediatricAdult/input_metadata.txt /cluster/home/t117652uhn/liver_ped_map/data

#grab figs
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/figures/*pdf /home/redgar/Documents/liver_ped_map/figures
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/figures/jpeg/*jpeg /home/redgar/Documents/liver_ped_map/figures/jpeg

scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/figures/IFALD_*pdf /home/redgar/Documents/liver_ped_map/figures
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/figures/jpeg/IFALD_*jpeg /home/redgar/Documents/liver_ped_map/figures/jpeg
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/figures/jpeg/IFALD_*png /home/redgar/Documents/liver_ped_map/figures/png


#grab data
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/QC_metrics.Rdata /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/adult_ped_integrated.rds /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/adult_ped_integrated_refinedlabels_withDropletQC.rds /home/redgar/Documents/liver_ped_map/data

scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/adult_ped_cellRefined_withDropletQC.rds /home/redgar/Documents/liver_ped_map/data

scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/adult_ped_cellRough.rds /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/d10x_adult_ped_raw.rds /home/redgar/Documents/liver_ped_map/data

scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/IFALD_adult_ped_integrated.rds /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/IFALD_d10x_adult_ped_raw.rds /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/IFALD_adult_ped_cellRefined_withDropletQC.rds /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/IFALD_adult_ped_SCINA_cell_labels.RData /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/IFALD_adult_ped_cellRough.rds /home/redgar/Documents/liver_ped_map/data

scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/IFALD_adult_ped_integrated.rds /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/IFALD_adult_ped_cellRough.rds /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/IFALD_d10x_adult_ped_raw.rds /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/IFALD_adult_ped_SCINA_cell_labels.RData /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/processed_data/IFALD_adult_ped_cellRefined_withDropletQC.rds /home/redgar/Documents/liver_ped_map/data



scp /home/redgar/Documents/liver_ped_map/data/IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data
scp /home/redgar/Documents/liver_ped_map/data/IFALD_adult_ped_cellRefined_withDropletQC.rds t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data


scp /home/redgar/Documents/tissue_MHCI/R_functions/pretty_plots.R t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn

scp /home/redgar/SynologyDrive/MacParland_Sai__C98_3pr_V3_1/pretty_plots.R t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn

scp -r /home/redgar/SynologyDrive/MacParland_Diana__SingleCell_C92_10Feb22_3pr_V3_1/outs t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn
scp /home/redgar/Documents/liver_ped_map/data/input_metadata.txt t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn

scp -r /home/redgar/SynologyDrive/MacParland_Sai__C98_3pr_V3_1 t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn

cp -r /cluster/home/t117652uhn/MacParland_Sai__C98_3pr_V3_1 

scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/adult_ped_diff_motecarlo.RData /home/redgar/Documents/liver_ped_map/data

scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/*adult_ped_diff_motecarlo.RData /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/*adult_ped_diff_motecarlo_1000.RData /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/*sex_diff_motecarlo_1000.RData /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/*adult_ped_diff_motecarlo_1000_covarSex.RData /home/redgar/Documents/liver_ped_map/data


scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/sitecheck.txt /home/redgar/Documents/tools


scp -r /media/redgar/Seagate Portable Drive/synology/MacParland_Sai__C98_3pr_V3_1 t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn
scp -r /media/redgar/Seagate\ Portable\ Drive/synology/MacParland_Sai__C98_3pr_V3_1 t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn

scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/cell_rough_maxmean.RData.RData /home/redgar/Documents/liver_ped_map/data



wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166504/suppl/GSE166504_cell_raw_counts.20220204.txt.gz'


scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/spatial_scRNAseq/data/reference_datasets/mouse_liver_GSE166504/GSE166504_cell_raw_counts_Hepatocyte_Chow_Animal2.csv /home/redgar/Documents/spatial_scRNAseq/data/reference_datasets/mouse_liver_GSE166504

scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/fetal_scores.RData /home/redgar/Documents/liver_ped_map/data


scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/d10x_adult_ped_raw_noSoupX.rds /media/redgar/Seagate\ Portable\ Drive/processed_data


scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/figures/fetal/*pdf /home/redgar/Documents/liver_ped_map/figures/fetal
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/figures/fetal/jpeg/*jpeg /home/redgar/Documents/liver_ped_map/figures/fetal/jpeg




mv /cluster/home/t117652uhn/liver_ped_map/data/IFALD_adult_ped_integrated.rds /cluster/projects/macparland/RE/PediatricAdult/processed_data
mv /cluster/home/t117652uhn/liver_ped_map/data/IFALD_adult_ped_cellRough.rds /cluster/projects/macparland/RE/PediatricAdult/processed_data
mv /cluster/home/t117652uhn/liver_ped_map/data/IFALD_d10x_adult_ped_raw.rds /cluster/projects/macparland/RE/PediatricAdult/processed_data

mv /cluster/home/t117652uhn/liver_ped_map/data/d10x_fetal_raw.rds /cluster/projects/macparland/RE/PediatricAdult/processed_data


scp /home/redgar/Documents/liver_ped_map/data/IFALD_B_cell_labels.rds t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data



scp -r t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/TA/Caudate_vs_Flush/Flush /media/redgar/Seagate\ Portable\ Drive/Caudate_Flush 
scp -r t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/TA/Caudate_vs_Flush/Caudate /media/redgar/Seagate\ Portable\ Drive/Caudate_Flush 



scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/QC_adult_ped_IFALD.out /home/redgar/Documents/liver_ped_map
