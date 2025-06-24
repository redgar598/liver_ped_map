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


source("scRNA_seq_scripts/00_pretty_plots.R")
source("scRNA_seq_scripts/00_plot_gene_exp.R")
source("scRNA_seq_scripts/00_fanciest_UMAP.R")

load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))


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
d10x.combined$FreshorFrozen<-NULL


colnames(d10x.combined@meta.data)[which(colnames(d10x.combined@meta.data)=="individual")]<-"sample"
d10x.combined$sample<-gsub("_caud3pr|_caud5pr|_bx5pr","", d10x.combined$sample)
d10x.combined$sample<-gsub("TLH","_TLH",d10x.combined$sample)
d10x.combined$sample<-gsub("NPC","_NPC",d10x.combined$sample)

d10x.combined$Tissue<-"Caudate"
d10x.combined$Tissue[which(d10x.combined$sample%in%c("C105","C102","C113"))]<-"Right Lobe"
d10x.combined$Tissue[which(d10x.combined$sample%in%c("IFALD006","IFALD030"))]<-"Liver"
d10x.combined$Tissue[which(d10x.combined$sample%in%c("IFALD073"))]<-"Left Lobe"

d10x.combined$Dissociation<-as.factor(d10x.combined$Perfused)
levels(d10x.combined$Dissociation)<-c("Manual","Perfusion")
d10x.combined$Perfused<-NULL

DefaultAssay(d10x.combined) <- "RNA"
d10x.combined <- NormalizeData(d10x.combined,scale.factor = 10000, normalization.method = "LogNormalize")

save(d10x.combined, file=here("/media/redgar/Seagate Portable Drive/HCA_Ped_data/","IFALD_adult_ped_integrated.rds"))


##############
### healthy only map
##############
load("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/IFALD_adult_ped_integrated_healthy_only.rds")


unique(d10x.combined_healthy$Age)
d10x.combined_healthy$Age[which(d10x.combined_healthy$Age%in%c(0.33,0.58))]<-"0-11 months"
d10x.combined_healthy$Age[which(d10x.combined_healthy$Age%in%c(2))]<-"1-4"
d10x.combined_healthy$Age[which(d10x.combined_healthy$Age%in%c(9,11,12))]<-"5-14"
d10x.combined_healthy$Age[which(d10x.combined_healthy$Age%in%c(16,17))]<-"15-19"

d10x.combined_healthy$file<-NULL
d10x.combined_healthy$orig.ident<-NULL
d10x.combined_healthy$BMI<-NULL
d10x.combined_healthy$Approx_bam_GB<-NULL
d10x.combined_healthy$old.ident<-NULL
d10x.combined_healthy$integrated_snn_res.0.5<-NULL
d10x.combined_healthy$CellType_rough<-NULL
d10x.combined_healthy$second_best_cell<-NULL
d10x.combined_healthy$SCINA_broad<-NULL
d10x.combined_healthy$cluster_consensus<-NULL
d10x.combined_healthy$age_id<-NULL
d10x.combined_healthy$SCINA_refined<-NULL
d10x.combined_healthy$FreshorFrozen<-NULL

##

colnames(d10x.combined_healthy@meta.data)[which(colnames(d10x.combined_healthy@meta.data)=="individual")]<-"sample"
d10x.combined_healthy$sample<-gsub("_caud3pr|_caud5pr|_bx5pr","", d10x.combined_healthy$sample)
d10x.combined_healthy$sample<-gsub("TLH","_TLH",d10x.combined_healthy$sample)
d10x.combined_healthy$sample<-gsub("NPC","_NPC",d10x.combined_healthy$sample)

d10x.combined_healthy$Tissue<-"Caudate"
d10x.combined_healthy$Tissue[which(d10x.combined_healthy$sample%in%c("C105","C102","C113"))]<-"Right Lobe"
d10x.combined_healthy$Tissue[which(d10x.combined_healthy$sample%in%c("IFALD006","IFALD030"))]<-"Liver"
d10x.combined_healthy$Tissue[which(d10x.combined_healthy$sample%in%c("IFALD073"))]<-"Left Lobe"

  
d10x.combined_healthy$Dissociation<-as.factor(d10x.combined_healthy$Perfused)
levels(d10x.combined_healthy$Dissociation)<-c("Manual","Perfusion")
d10x.combined_healthy$Perfused<-NULL

DefaultAssay(d10x.combined_healthy) <- "RNA"
d10x.combined_healthy <- NormalizeData(d10x.combined_healthy,scale.factor = 10000, normalization.method = "LogNormalize")

save(d10x.combined_healthy, file=here("/media/redgar/Seagate Portable Drive/HCA_Ped_data/","Adult_ped_integrated.rds"))

