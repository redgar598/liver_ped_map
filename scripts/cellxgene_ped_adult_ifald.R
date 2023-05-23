load(here("data","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

SaveH5Seurat(d10x.combined, filename = here("data","adult_ped_ifald.h5Seurat"))
Convert(here("data","adult_ped_ifald.h5Seurat"), dest = "h5ad")
