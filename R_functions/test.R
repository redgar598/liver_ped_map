library(here)
library(Seurat)

load(here("data/adult_ped_integrated_refinedlabels.rds"))

cell_label<-d10x.combined@meta.data
save(cell_label, file=paste(here("data/"),"adult_ped_cellRefined.rds", sep=""))
