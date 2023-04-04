library(SeuratDisk)
library(Seurat)
library(here)
# 
#                                   #https://developmental.cellatlas.io/fetal-liver
#                                   fetal_liver_h5ad <- Convert(here("/media/redgar/Seagate Portable Drive/fetal_liver/download.h5ad"), ".h5seurat", overwrite= TRUE)
#                                   fetal_liver_h5ad <- Convert(here("../../../projects/macparland/RE/PediatricAdult/fetal_liver/download.h5ad"), ".h5seurat", overwrite= TRUE)
#                                   fetal_liver_h5ad <- Convert(here("../../Downloads/fetal_liver_alladata.h5ad"), ".h5seurat", overwrite= TRUE)
#                                   
#                                   fetal_liver <- LoadH5Seurat(fetal_liver_h5ad, assays = "RNA")
#                                   
#                                   print(dim(fetal_liver))
#                                   print(fetal_liver)
#                                   
#                                   saveRDS(fetal_liver, here("../../../projects/macparland/RE/PediatricAdult/fetal_liver","fetal_liver.rds"))
#                                   
# 


# library(anndata)
# #https://developmental.cellatlas.io/fetal-liver
# 
# fetal_liver_h5ad <- read_h5ad(here("/media/redgar/Seagate Portable Drive/fetal_liver/download.h5ad"))
# write_h5ad(fetal_liver_h5ad, here("/media/redgar/Seagate Portable Drive/fetal_liver/fetal_liver.h5ad"))
# 

fetal_liver_h5ad <- Convert(here("../../../projects/macparland/RE/PediatricAdult/fetal_liver/fetal_liver.h5ad"), ".h5seurat", overwrite= TRUE)
fetal_liver <- LoadH5Seurat(fetal_liver_h5ad, assays = "RNA")

print(dim(fetal_liver))
print(fetal_liver)

saveRDS(fetal_liver, here("../../../projects/macparland/RE/PediatricAdult/fetal_liver","fetal_liver.rds"))
