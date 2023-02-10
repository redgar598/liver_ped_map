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


load(here("data","adult_ped_integrated_refinedlabels_withDropletQC.rds"))



d10x.combined_adult<-subset(d10x.combined, subset = AgeGroup == "Adult")
rm(d10x.combined)
gc()
d10x.combined_adult <- RunPCA(d10x.combined_adult, npcs = 30, verbose = FALSE)
d10x.combined_adult <- RunUMAP(d10x.combined_adult, reduction = "pca", dims = 1:30)

all_refined_cluster_umap_nolab<-DimPlot(d10x.combined_adult, reduction = "umap", pt.size=0.25, group.by = "CellType_refined")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-11, y=-12, label=paste0("n = ",comma(ncol(d10x.combined_adult))))
all_refined_cluster_umap_nolab
save_plts(all_refined_cluster_umap_nolab, "refined_cellType_map_nolabel_adultsOnly", w=7,h=5)