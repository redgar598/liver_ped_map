library(here)
library(Seurat)


dataset_loc <- here("../../../projects/macparland/RE/PediatricAdult")

samples<-list.files(dataset_loc)
samples<-samples[-grep("meta",samples)]
print(samples)

meta<-read.table(here(dataset_loc,"input_metadata.txt"), header=T)
#meta<-read.table(here("data/input_metadata.txt"), header=T)


d10x.list <- sapply(1:length(samples), function(y){
  print(file.path(dataset_loc,paste(samples[y], sep=""),"filtered_feature_bc_matrix"))
  d10x <- Read10X(file.path(dataset_loc,paste(samples[y], sep=""),"filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),samples[y],sep="-")
  print(dim(d10x))
  #' Initialize the Seurat object with the raw (non-normalized data).
  d10x<-CreateSeuratObject(counts = d10x, project = "ped_adult_map", min.cells = 3, min.features = 0)
  
  print(head(colnames(d10x)))
  #add meta data to each seurat object
  meta_cell<-data.frame(cell=colnames(d10x), individual=sapply(colnames(d10x), function(x) strsplit(x,"-")[[1]][2]))
  print(head(meta_cell))
  meta_cell_add<-merge(meta_cell, meta, by.x="individual", by.y="Sample_ID")
  meta_cell_add<-meta_cell_add[match(colnames(d10x), meta_cell_add$cell),]
  print(identical(meta_cell_add$cell, colnames(d10x)))
  print(head(meta_cell_add))
  
  d10x<- AddMetaData(d10x, meta_cell_add$donor, col.name = "individual")
  d10x
  })

d10x.list
