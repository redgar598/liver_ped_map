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


options(stringsAsFactors = FALSE)

source("R_functions/pretty_plots.R")


### load integrate for UMAP etc
load(here("data","adult_ped_integrated.rds"))

d10x.combined_CD3<-subset(d10x.combined, subset = CellType_rough %in% c("CD3_Tcell"))
d10x.combined_CD3 <- RunPCA(d10x.combined_CD3, npcs = 30, verbose = FALSE)
d10x.combined_CD3 <- RunUMAP(d10x.combined_CD3, reduction = "pca", dims = 1:30)
d10x.combined_CD3 <- FindNeighbors(d10x.combined_CD3, reduction = "pca", dims = 1:30)
d10x.combined_CD3 <- FindClusters(d10x.combined_CD3, resolution = 0.3)
CD3_umap<-DimPlot(d10x.combined_CD3, label=T)
CD3_umap

##############
### entropy in clusters
##############
source(here("R_functions/entropy_d10x.R"))

plt_entropy_individual<-entropy_d10(d10x.combined_CD3, "individual")
plt_entropy_age<-entropy_d10(d10x.combined_CD3, "AgeGroup")
plt_entropy_chem<-entropy_d10(d10x.combined_CD3, "Chemistry")

entropy_plt(plt_entropy_individual, "individual", d10x.combined_CD3)
entropy_plt(plt_entropy_age, "AgeGroup", d10x.combined_CD3)
entropy_plt(plt_entropy_chem, "Chemistry", d10x.combined_CD3)



plt_entropy_individual<-entropy_d10(d10x.combined, "individual")
plt_entropy_age<-entropy_d10(d10x.combined, "AgeGroup")
plt_entropy_chem<-entropy_d10(d10x.combined, "Chemistry")

entropy_plt(plt_entropy_individual, "individual", d10x.combined)
entropy_plt(plt_entropy_age, "AgeGroup", d10x.combined)
entropy_plt(plt_entropy_chem, "Chemistry", d10x.combined)
