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

source("scripts/00_pretty_plots.R")

load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
cell_label_IFALD<-cell_label

# load(here("data","adult_ped_cellRefined_withDropletQC.rds"))
# cell_label_og<-cell_label

load(here("data","adult_ped_integrated_refinedlabels_withDropletQC.rds"))
meta_original<-d10x.combined@meta.data


cell_label_IFALD$index<-rownames(cell_label_IFALD)
cell_label_IFALD$index<-sapply(1:length(cell_label_IFALD$index), function(x){
  paste(strsplit(cell_label_IFALD$index[x],"-")[[1]][1], "_",cell_label_IFALD$individual[x], sep="")
})

meta_original$cell<-sapply(1:length(meta_original$cell), function(x){
  paste(strsplit(meta_original$cell[x],"-")[[1]][1], "_",meta_original$individual[x], sep="")
})

dim(cell_label_IFALD)
length(which(cell_label_IFALD$index %in% meta_original$cell))
dim(meta_original)
length(which(meta_original$cell %in% cell_label_IFALD$index))

meta_original<-meta_original[which(meta_original$cell %in% cell_label_IFALD$index),]
cell_label_IFALD<-cell_label_IFALD[which(cell_label_IFALD$index %in% meta_original$cell),]
cell_label_IFALD<-cell_label_IFALD[match(meta_original$cell, cell_label_IFALD$index),]
identical(meta_original$cell, cell_label_IFALD$index)

meta_original<-meta_original[,c("cell","CellType_refined")]
colnames(meta_original)<-c("cell","CellType_original")
meta_original$index<-rownames(meta_original)

cell_label_IFALD<-merge(meta_original,cell_label_IFALD[,c("index","CellType_refined")], by.x="cell", by.y="index")
rownames(cell_label_IFALD)<-cell_label_IFALD$index

d10x.combined <- AddMetaData(d10x.combined, metadata = cell_label_IFALD)



tile_plt<-melt(table(d10x.combined$CellType_refined, d10x.combined$CellType_original))

tile_plt$value[which(tile_plt$value==0)]<-NA
ggplot(tile_plt, aes(Var1, Var2, fill=value))+geom_tile()+geom_text(aes(label=value))+
  scale_fill_continuous_sequential(palette = "Blues", rev=T, na.value = "white")+theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


DimPlot(d10x.combined, group.by = "CellType_original")+colscale_cellType+ggtitle("")
DimPlot(d10x.combined, group.by = "CellType_refined")+colscale_cellType+ggtitle("")
