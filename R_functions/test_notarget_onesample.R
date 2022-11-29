library(here)
library(Seurat)


dataset_loc <- here("../../../projects/macparland/RE/PediatricAdult")

samples<-list.files(dataset_loc)
print(samples)

y=1


print(file.path(dataset_loc,paste(samples[y], sep=""),"filtered_feature_bc_matrix"))
d10x_read <- Read10X(file.path(dataset_loc,paste(samples[y], sep=""),"filtered_feature_bc_matrix"))
colnames(d10x_read) <- paste(sapply(strsplit(colnames(d10x_read),split="-"),'[[',1L),samples[y],sep="-")

#' Initialize the Seurat object with the raw (non-normalized data).
d10x <- CreateSeuratObject(counts = d10x_read, project = "ped_adult_map", min.cells = 3, min.features = 0)
d10x