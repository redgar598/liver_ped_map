library(SeuratDisk)
library(Seurat)
library(here)

#https://developmental.cellatlas.io/fetal-liver
fetal_liver_h5ad <- Convert(here("/media/redgar/Seagate Portable Drive/fetal_liver/download.h5ad"), ".h5seurat", overwrite= TRUE)
fetal_liver <- LoadH5Seurat(fetal_liver_h5ad, assays = "RNA")

print(dim(fetal_liver))
print(fetal_liver)

#saveRDS(here("data/adult_intestine","local.rds"))
