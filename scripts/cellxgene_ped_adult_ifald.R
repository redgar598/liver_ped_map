load(here("data","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

DefaultAssay(d10x.combined)<-"RNA"

SaveH5Seurat(d10x.combined, filename = here("data","adult_ped_ifald.h5Seurat"), overwrite=T)
Convert(here("data","adult_ped_ifald.h5Seurat"), dest = "h5ad",overwrite=T)
