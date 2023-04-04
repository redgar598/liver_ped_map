
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


# submit job on login node?
cd liver_ped_map
git pull
sbatch exploratory_run.sh
sbatch de_monte_carlo_notparallel.sh
sbatch metacell.sh


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
for cellType_index in {1..12}
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
### Moving files
#################

### big raw data move
#ssh t117652uhn@h4huhndata1.uhnresearch.ca
scp /home/redgar/Documents/liver_ped_map/data/data_transfer_updated_jan16_2023.csv t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw
scp /home/redgar/Documents/liver_ped_map/data/data_transfer_updated_mar20_2023.csv t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw

scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/C39_NPC_june6_2017 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/
scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/C39_TLH_june6_2017 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/
scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/C54_3prV2_12Apr18 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/
scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/X210026__SingleCell_C88_14Oct21_3pr_V3_1 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/
scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/MacParland_Jawairia__C104_Bx_5pr_V2 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/
scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/MacParland_Sonya__C68_Total_liver t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/
scp -r /media/redgar/Seagate\ Portable\ Drive/ped_liver_map_raw/MacParlnd_Sai__C97_3pr_V3_1 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/ped_liver_map_raw/

scp /media/redgar/Seagate\ Portable\ Drive/fetal_liver/download.h5ad t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/fetal_liver




## general file moving moving files
scp /home/redgar/Documents/tissue_MHCI/R_functions/pretty_plots.R t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/input_metadata.txt /home/redgar/Documents/liver_ped_map/data


cp -R /cluster/projects/macparland/DN/PediatricAdult /cluster/projects/macparland/RE


cp /cluster/projects/macparland/RE/PediatricAdult/input_metadata.txt /cluster/home/t117652uhn/liver_ped_map/data

#grab figs
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/figures/*pdf /home/redgar/Documents/liver_ped_map/figures
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/figures/jpeg/*jpeg /home/redgar/Documents/liver_ped_map/figures/jpeg

#grab data
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/QC_metrics.Rdata /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/adult_ped_integrated_refinedlabels_withDropletQC.rds /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/adult_ped_cellRefined_withDropletQC.rds /home/redgar/Documents/liver_ped_map/data

scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/adult_ped_cellRough.rds /home/redgar/Documents/liver_ped_map/data
scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/d10x_adult_ped_raw.rds /home/redgar/Documents/liver_ped_map/data

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

scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/DGE_compare.Rdata /home/redgar/Documents/liver_ped_map/data


scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data/d10x_adult_ped_raw_noSoupX.rds /media/redgar/Seagate\ Portable\ Drive/processed_data