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
