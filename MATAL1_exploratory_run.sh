#!/bin/bash
#SBATCH --mem=150G
#SBATCH -J MATAL1_QC
#SBATCH -p veryhimem
#SBATCH -c 32 
#SBATCH -N 1 
#SBATCH -t 5-00:00:00
#SBATCH --output=%x.out


module load R/4.2.1



#Rscript scripts/test.R
##Rscript R_functions/test_notarget.R

#Rscript scripts/01_ped_liver_QC_with_dropletQC_SCINA_IFALD.R

#Rscript scripts/01_ped_liver_QC_with_dropletQC_SCINA_IFALD_PBMC_only.R

#Rscript scripts/03_PBMC_correction_IFALD.R





Rscript scripts/MALAT1_QC_Zoe.R




# [1] "6550 cell originally 4568 after MALAT1 filter"
# [1] "11071 cell originally 614 after MALAT1 filter"
# [1] "12332 cell originally 6196 after MALAT1 filter"
# [1] "10620 cell originally 6076 after MALAT1 filter"
# [1] "1684 cell originally 753 after MALAT1 filter"
# [1] "7936 cell originally 5871 after MALAT1 filter"
# [1] "12821 cell originally 470 after MALAT1 filter"
# [1] "19340 cell originally 1508 after MALAT1 filter"
# [1] "10285 cell originally 684 after MALAT1 filter"
# [1] "15993 cell originally 7236 after MALAT1 filter"
# [1] "65540 cell originally 2213 after MALAT1 filter"
# [1] "10000 cell originally 4 after MALAT1 filter"
# [1] "10217 cell originally 5696 after MALAT1 filter"
# [1] "21057 cell originally 6798 after MALAT1 filter"
# [1] "11111 cell originally 1475 after MALAT1 filter"
# [1] "7604 cell originally 583 after MALAT1 filter"
# [1] "6515 cell originally 66 after MALAT1 filter"
# [1] "19349 cell originally 3351 after MALAT1 filter"
# [1] "56622 cell originally 2 after MALAT1 filter"
# [1] "8404 cell originally 316 after MALAT1 filter


# scp t117652uhn@h4huhnlogin1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/figures/*MALAT1.pdf /home/redgar/Documents/liver_ped_map/figures



#