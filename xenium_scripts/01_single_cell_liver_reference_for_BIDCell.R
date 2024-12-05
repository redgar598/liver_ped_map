## Load Libraries
library(here)
library(Seurat)

library(SCINA)
library(reshape2)
library(dplyr)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

library(sp)
#library(rgeos)
library(scales)
library(gridExtra)
#library(rhdf5)

source("scripts/00_pretty_plots.R")
options(future.globals.maxSize= 8000 * 1024^2)


#only healthy
d10x<-readRDS(file = here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_d10x_adult_ped_raw.rds"))

load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

d10x <- subset(d10x, subset = Treatment == "Healthy")

levels(d10x$CellType_refined)[which(levels(d10x$CellType_refined)=="Myeloid Erythrocytes\n(phagocytosis)" )]<-"Myeloid Erythrocytes" 
levels(d10x$CellType_refined)[which(levels(d10x$CellType_refined)=="Macrophage\n(MHCII high)" )]<-"Macrophage MHCII high" 

d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

#############
## load xenium data
#############
samples<-c("C94_2","C94_3","C94_4")

d10x.list <- sapply(samples, function(smple){
  #load(file=paste(here("data/"),smple, "_object_raw.RData",sep=""))
  load(file=paste("/media/redgar/Seagate Portable Drive/xenium_liver/",smple,"_object_raw.RData",sep=""))
  xenium.obj$sample<-smple
  xenium.obj
})

d10x.list

xenium.obj <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "xenium_liver")
rm(d10x.list)
gc()


metadata_add<-xenium.obj@meta.data

## 477 genes in both
reference_subsample <- subset(d10x, features = intersect(rownames(d10x), rownames(xenium.obj)))
reference_subsample




###############
## reference csv
###############

## The reference csv file contains average expressions for all of the genes in the spatial transcriptomic dataset of different cell types. 
## You may choose an appropriate list of cell types to include for your data.

head(AverageExpression(object = reference_subsample, group.by = c('CellType_refined'))$RNA)

cell_typeMeans<-t(as.data.frame(AverageExpression(object = reference_subsample, group.by = c('CellType_refined'))$RNA))
cell_typeMeans<-as.data.frame(cell_typeMeans)

# No spaces in cell type names?
rownames(cell_typeMeans)<-gsub(" |-|\\+", "_", rownames(cell_typeMeans))
unique(rownames(cell_typeMeans))

cell_typeMeans_allindex<-cell_typeMeans
cell_typeMeans_allindex$ct_idx<-seq(1:nrow(cell_typeMeans_allindex))
cell_typeMeans_allindex$ct_idx<-cell_typeMeans_allindex$ct_idx-1
cell_typeMeans_allindex$cell_type<-rownames(cell_typeMeans_allindex)
cell_typeMeans_allindex$atlas<-"Ped_map"
rownames(cell_typeMeans_allindex)<-cell_typeMeans_allindex$ct_idx

write.csv(cell_typeMeans_allindex, file="/media/redgar/Seagate Portable Drive/xenium_liver/ped_map_sc_means.csv", quote = F)




###############
## Positive Negative Markers
###############

## The positive and negative markers were those with expressions in the highest and lowest 10 percentile 
## for each cell type of a tissue sample. We found that removing positive markers that were common to at 
## least a third of cell types in each dataset was appropriate across various datasets. 


                  ## remove any genes 0 in all cell types (useless)
                  #cell_typeMeans<-cell_typeMeans[,which(colSums(cell_typeMeans)!=0)]


positive_genes <- data.frame(matrix(0, nrow = nrow(cell_typeMeans), ncol = ncol(cell_typeMeans)))
colnames(positive_genes) <- colnames(cell_typeMeans)
rownames(positive_genes) <- rownames(cell_typeMeans)

lapply(1:nrow(cell_typeMeans), function(x){
  positive<-quantile(cell_typeMeans[x,], 0.9)
  positive<-positive$`90%`
  positive_genes[x,which(cell_typeMeans[x,]>=positive)]<<-1})
rowSums(positive_genes)
write.csv(positive_genes, file="/media/redgar/Seagate Portable Drive/xenium_liver/ped_map_positive_genes.csv", quote = F)



negative_genes <- data.frame(matrix(0, nrow = nrow(cell_typeMeans), ncol = ncol(cell_typeMeans)))
colnames(negative_genes) <- colnames(cell_typeMeans)
rownames(negative_genes) <- rownames(cell_typeMeans)


lapply(1:nrow(cell_typeMeans), function(x){
  negative<-quantile(cell_typeMeans[x,], 0.1)
  negative<-negative$`10%`
  negative_genes[x,which(cell_typeMeans[x,]<=negative)]<<-1})
rowSums(negative_genes)
write.csv(negative_genes, file="/media/redgar/Seagate Portable Drive/xenium_liver/ped_map_negative_genes.csv", quote = F)





