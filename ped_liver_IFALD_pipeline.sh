# submit job on login node
cd liver_ped_map
git pull
sbatch IFALD_exploratory_run.sh
sbatch cpdb_liver.sh
sbatch Fetal_IFALD_exploratory_run.sh


###################################3


01_ped_liver_QC_with_dropletQC_SCINA_IFALD.R
01_fetal_liver.R

02_ped_liver_QC_with_dropletQC_SCINA_IFALD_PBMC_only.R 

03_IFALD_myeloid.R
03_IFALD_Bcells.R 
IFALD_presentation_plots.R
04_PBMC_correction_IFALD 
03_IFALD_HSC.R

cellphonedb_datasetup.R
cpdb_network_plt.R
cpdb_network_plt_IFALD.R


## running




## to do 
02_IFALD_ped_adult_fetal_integration.R
04_Fetal_PCA_roughwork
cellxgene_ped_adult_ifald.R
06_adolescent.R

cpdb_age

	cpdb_IFALD

	#03_IFALD_Tcells.R
	#03 colangiocyted
	#03 LSEC
	#03 hepatocytes
	#

	#04_differential_expression_montecarlo.R
	#04_differential_expression_montecarlo_withFetal_IFALD.R
	#05_differential_expression_age_interpretation.R

	06_cell_composistion.R


