### Load libraries
library(here)
library(Seurat)
library(SeuratDisk)

load(here("data","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

DefaultAssay(d10x.combined)<-"RNA"

str(d10x.combined@meta.data)

i <- sapply(d10x.combined@meta.data, is.factor)
d10x.combined@meta.data[i] <- lapply(d10x.combined@meta.data[i], as.character)

SaveH5Seurat(d10x.combined, filename = here("data","adult_ped_ifald.h5Seurat"), overwrite=T)
Convert(here("data","adult_ped_ifald.h5Seurat"), dest = "h5ad",overwrite=T)



##############
## myeloid only
##############

load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

d10x.combined_myeloid<-subset(d10x.combined, subset = CellType_refined %in% c("Mono-Mac","Macrophage\n(MHCII high)","KC Like","CDC1","Cycling Myeloid","Myeloid Erythrocytes\n(phagocytosis)"))
rm(d10x.combined)
gc()
d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)
d10x.combined_myeloid <- RunUMAP(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindNeighbors(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindClusters(d10x.combined_myeloid, resolution = 0.3)
d10x.combined_myeloid <- FindClusters(d10x.combined_myeloid, resolution = 0.6)


DefaultAssay(d10x.combined_myeloid)<-"RNA"
str(d10x.combined_myeloid@meta.data)

i <- sapply(d10x.combined_myeloid@meta.data, is.factor)
d10x.combined_myeloid@meta.data[i] <- lapply(d10x.combined_myeloid@meta.data[i], as.character)

SaveH5Seurat(d10x.combined_myeloid, filename = here("data","adult_ped_ifald_myeloid.h5Seurat"), overwrite=T)
Convert(here("data","adult_ped_ifald_myeloid.h5Seurat"), dest = "h5ad",overwrite=T)
