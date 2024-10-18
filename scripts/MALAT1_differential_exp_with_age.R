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



##############
## MALAT1
##############
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_d10x_adult_ped_raw.rds"))

load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))


cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

de_celltype<-lapply(1:length(as.character(unique(d10x$CellType_refined))), function(x){
  celltype<-as.character(unique(d10x$CellType_refined))[x]
  print(celltype)
  d10x_raw_KC<-subset(d10x, subset = CellType_refined %in% celltype)
  
  Idents(d10x_raw_KC)<-d10x_raw_KC$age_condition
  table(d10x_raw_KC$age_condition)
  
  ## age differential
  de<-FindMarkers(d10x_raw_KC, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F, features = c("MALAT1","B2M"), logfc.threshold = 0,min.pct = 0,min.cells.feature = 3 )
  de$celltype<-celltype
  de[which(rownames(de)=="MALAT1"),]
})

de_celltype<-do.call(rbind, de_celltype)

de_celltype[which(de_celltype$p_val_adj<0.005 & abs(de_celltype$avg_log2FC)>0.5),]

vln_malat1<-VlnPlot(d10x, features = "MALAT1", split.by = "age_condition",group.by = "CellType_refined")
save_plts(vln_malat1, "MALAT1_age", w=15, h=4)
