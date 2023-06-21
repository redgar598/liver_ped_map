### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(reshape2)
library(gtools)
library(colorspace)
library(cowplot)
library(SeuratDisk)


source("scripts/00_pretty_plots.R")
source("scripts/00_fanciest_UMAP.R")
source("scripts/00_plot_gene_exp.R")


###############
## Ped Healthy Output
###############
means<-read.table(here("data/cellphonedb/statistical_analysis_significant_means_06_20_2023_17:10:20.txt"), sep="\t", header=T)
means[1:5,1:10]


CCR_sig<-means[grep("CCR",means$interacting_pair), ]
lapply(1:nrow(CCR_sig), function(x) cbind(CCR_sig[x,c(1:12)],CCR_sig[x,13:496][which(!(is.na(CCR_sig[x,13:496])))]) )

CCL_sig<-means[grep("CCL3|CCL4",means$interacting_pair), ]
CCL_sig[,c(1:12)]
lapply(1:nrow(CCL_sig), function(x) cbind(CCL_sig[x,c(1:12)],CCL_sig[x,13:496][which(!(is.na(CCL_sig[x,13:496])))]) )

unlist(means[grep("CCL3_CCR5",means$interacting_pair), ])

unlist(means[grep("CCL15_CCR1",means$interacting_pair), ])



load(here("data","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))



