library(here)
library(Seurat)


dataset_loc <- here("../../../projects/macparland/RE/PediatricAdult")

samples<-list.files(dataset_loc)
print(samples)

d10x.data <- sapply(1:length(samples), function(y){
  print(file.path(dataset_loc,paste(samples[y], sep=""),"filtered_feature_bc_matrix"))
  d10x <- Read10X(file.path(dataset_loc,paste(samples[y], sep=""),"filtered_feature_bc_matrix"))
  #colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})
d10x.data<-do.call("cbind",d10x.data)

#' Initialize the Seurat object with the raw (non-normalized data).
d10x <- CreateSeuratObject(counts = d10x.data, project = "ped_adult_map", min.cells = 3, min.features = 0)
d10x