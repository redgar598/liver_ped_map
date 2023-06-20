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

load(here("data","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))




## Counts
DefaultAssay(d10x.combined)<-"RNA"
d10x.combined.subset<-subset(d10x.combined, subset = age_condition == "Ped Healthy")
rm(d10x.combined)
gc()


              #levels(d10x.combined.subset$CellType_refined)[c(6,12)]<-c("CLNK_T_cell","KC_like")
              
              ##just two cell types
              # d10x.combined.subset<-subset(d10x.combined.subset, subset = CellType_refined %in% c("KC_like","CLNK_T_cell"))
              # d10x.combined.subset <- d10x.combined.subset[, sample(colnames(d10x.combined.subset), size = 200, replace=F)]


d10x.combined.subset <- NormalizeData(d10x.combined.subset, scale.factor = 10000, normalization.method = "LogNormalize")

SaveH5Seurat(d10x.combined.subset, filename = here("../../../projects/macparland/RE/PediatricAdult/processed_data","ped_healthy.h5Seurat"), overwrite = T)
Convert(here("../../../projects/macparland/RE/PediatricAdult/processed_data","ped_healthy.h5Seurat"), dest = "h5ad", overwrite = T)



## meta data
meta_data<-d10x.combined.subset@meta.data[,c("age_condition","CellType_refined")]
meta_data$age_condition<-as.character(meta_data$age_condition)
meta_data$Cell<-rownames(meta_data)

meta_data$CellType_refined<-gsub("\n"," ", meta_data$CellType_refined)

write.csv(meta_data[,c("Cell","CellType_refined")], 
          file=paste(here("../../../projects/macparland/RE/PediatricAdult/processed_data/"),"Ped_healthy_meta.csv", sep=""),  
          row.names = F,  quote = F)


## microenvironment
micro_env<-data.frame(cell_type = unique(meta_data$CellType_refined), microenvironment = c("Env1"))

write.csv(micro_env, 
          file=paste(here("../../../projects/macparland/RE/PediatricAdult/processed_data/"),"Ped_healthy_microenv.csv", sep=""),  
          row.names = F,  quote = F)
