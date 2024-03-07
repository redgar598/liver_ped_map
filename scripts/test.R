### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(reshape2)
library(gtools)
library(SoupX)
library(colorspace)
library(cowplot)
library(DropletQC)
library(anndata)
library(RColorBrewer)

source("scripts/00_pretty_plots.R")
source("scripts/00_entropy_d10x.R")



load(here("../../../projects/macparland/RE/PediatricAdult/processed_dataFetal_IFALD_adult_ped_integrated.rds"))


d10x.fetal_ped_IFALD$CellType_harmonized<-d10x.fetal_ped_IFALD$CellType_refined
levels(d10x.fetal_ped_IFALD$CellType_harmonized)[which(levels(d10x.fetal_ped_IFALD$CellType_harmonized)=="Kupffer Cell")]<-"KC Like"
levels(d10x.fetal_ped_IFALD$CellType_harmonized)[which(levels(d10x.fetal_ped_IFALD$CellType_harmonized)%in%c("Mono-Mac","Monocyte-DC precursor","Monocyte" ))]<-"Mono-Mac"
levels(d10x.fetal_ped_IFALD$CellType_harmonized)[which(levels(d10x.fetal_ped_IFALD$CellType_harmonized)%in%c("DC2","DC1" ))]<-"Macrophage\n(MHCII high)"
levels(d10x.fetal_ped_IFALD$CellType_harmonized)[which(levels(d10x.fetal_ped_IFALD$CellType_harmonized)%in%c("NK-like cells","NK" ))]<-"NK cell"

d10x.combined_myeloid<-subset(d10x.fetal_ped_IFALD, subset = CellType_refined %in% c("RR Myeloid","KC Like","Macrophage\n(MHCII high)","Cycling Myeloid",
                                                                                     "CDC1","VCAM1+ Erythroblastic Island Macrophage",
                                                                                     "pDC precursor","Neutrophil-myeloid progenitor","Mono-NK",
                                                                                     "Myeloid Erythrocytes\n(phagocytosis)","HSC/MPP","Monocyte-DC precursor","Mono-Mac",
                                                                                     "Kupffer Cell","Neutrophil-myeloid progenitor","Mono-NK",
                                                                                     "VCAM1+ Erythroblastic Island Macrophage","Monocyte","Erythroblastic Island Macrophage"))
d10x.combined_bcell<-subset(d10x.fetal_ped_IFALD, subset = CellType_harmonized %in% c("Mature B-cells","pro B cell","pre pro B cell","B cell",
                                                                                      "pre B cell","Plasma cells"))
d10x.combined_tcell<-subset(d10x.fetal_ped_IFALD, subset = CellType_harmonized %in% c("NK cell","CD3+ T-cells","CLNK T-cells","Cycling T-cells",
                                                                                      "ILC precursor","Early lymphoid/T lymphocyte","Mono-NK",
                                                                                      "gd T-cells",""))
d10x.combined_HSC<-subset(d10x.fetal_ped_IFALD, subset = CellType_harmonized %in% c("HSC"))

save(d10x.combined_bcell, file=here("../../../projects/macparland/RE/PediatricAdult/processed_data/Fetal_IFALD_adult_ped_integrated_bcell_only.RData"))
save(d10x.combined_tcell, file=here("../../../projects/macparland/RE/PediatricAdult/processed_data/Fetal_IFALD_adult_ped_integrated_tcell_only.RData"))
save(d10x.combined_HSC, file=here("../../../projects/macparland/RE/PediatricAdult/processed_data/Fetal_IFALD_adult_ped_integrated_HSC_only.RData"))


########
## Myeloid overlapping
########

d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)
d10x.combined_myeloid <- RunUMAP(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindNeighbors(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindClusters(d10x.combined_myeloid, resolution = 0.2)


myeloid_cluster_umap<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=T)
myeloid_cluster_umap
save_plts(myeloid_cluster_umap, "IFALD_fetal_myeloid_cluster_umap", w=5,h=4)

myeloid_cluster_umap_individual<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=T, split.by = "age_condition", ncol=2)
myeloid_cluster_umap_individual
save_plts(myeloid_cluster_umap_individual, "IFALD_fetal_myeloid_cluster_umap_individual", w=10,h=6)

DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=T,
        group.by = "CellType_refined",split.by = "age_condition", ncol=2)+colscale_cellType_fetal

DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=T,
        group.by = "CellType_refined")+colscale_cellType_fetal


d10x.combined_myeloid$CellType_harmonized<-as.factor(d10x.combined_myeloid$CellType_refined)
levels(d10x.combined_myeloid$CellType_harmonized)[which(levels(d10x.combined_myeloid$CellType_harmonized)=="Kupffer Cell")]<-"KC Like"
# levels(d10x.combined_myeloid$CellType_harmonized)[which(levels(d10x.combined_myeloid$CellType_harmonized)%in%c("Mono-Mac"))]<-"Mono-Mac"
# levels(d10x.combined_myeloid$CellType_harmonized)[which(levels(d10x.combined_myeloid$CellType_harmonized)%in%c("Monocyte" ))]<-"Monocyte"


myeloid_cluster_umap<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=F,
                              group.by = "CellType_harmonized",split.by = "age_condition", ncol=2)+colscale_cellType_fetal_combo
myeloid_cluster_umap
save_plts(myeloid_cluster_umap, "IFALD_fetal_myeloid_cluster_umap_groups", w=10,h=6)

######## Subset to just the overlapping cell types
d10x.combined_myeloid<-subset(d10x.combined_myeloid, subset = CellType_harmonized %in% c("Mono-Mac","KC Like","Macrophage\n(MHCII high)","CDC1","Monocyte"))
d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)
d10x.combined_myeloid <- RunUMAP(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindNeighbors(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindClusters(d10x.combined_myeloid, resolution = 0.2)

save(d10x.combined_myeloid, file=here("../../../projects/macparland/RE/PediatricAdult/processed_data/Fetal_IFALD_adult_ped_integrated_myeloid_only.RData"))

