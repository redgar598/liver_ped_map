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
library(RColorBrewer)


source("scripts/00_pretty_plots.R")
source("scripts/00_plot_gene_exp.R")
source("scripts/00_fanciest_UMAP.R")

load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

fanciest_UMAP(d10x.combined,"HSC",F)

unique(d10x.combined$Age)
d10x.combined$Age[which(d10x.combined$Age%in%c(0.33,0.58))]<-"0-11 months"
d10x.combined$Age[which(d10x.combined$Age%in%c(2))]<-"1-4"
d10x.combined$Age[which(d10x.combined$Age%in%c(9,11,12))]<-"5-14"
d10x.combined$Age[which(d10x.combined$Age%in%c(16,17))]<-"15-19"

d10x.combined$file<-NULL
d10x.combined$orig.ident<-NULL
d10x.combined$BMI<-NULL
d10x.combined$Approx_bam_GB<-NULL
d10x.combined$old.ident<-NULL
d10x.combined$integrated_snn_res.0.5<-NULL
d10x.combined$CellType_rough<-NULL
d10x.combined$second_best_cell<-NULL
d10x.combined$SCINA_broad<-NULL
d10x.combined$cluster_consensus<-NULL
d10x.combined$age_id<-NULL
d10x.combined$SCINA_refined<-NULL

save(d10x.combined, file=here("/media/redgar/Seagate Portable Drive/HCA_Ped_data/","IFALD_adult_ped_integrated.rds"))


load("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/IFALD_adult_ped_integrated_healthy_only.rds")

fanciest_UMAP(d10x.combined,"HSC",F)
