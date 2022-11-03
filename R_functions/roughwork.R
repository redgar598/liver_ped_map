### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(reshape2)
library(gtools)
library(SoupX)




source("R_functions/pretty_plots.R")


dataset_loc <- here("../../../projects/macparland/RE/PediatricAdult")

samples<-list.files(dataset_loc)
samples<-samples[-grep("meta",samples)]
print(samples)

#meta<-read.table(here("data/input_metadata.txt"), header=T)
meta<-read.table(here(dataset_loc,"input_metadata.txt"), header=T)
meta$Sample_ID[which(meta$Sample_ID=="C85_caud3pr")]<-"C85_caud5pr"

y=1
  print(file.path(dataset_loc,paste(samples[y], sep=""),"filtered_feature_bc_matrix"))
  d10x <- load10X(file.path(dataset_loc,paste(samples[y], sep=""),"filtered_feature_bc_matrix"))
  print(class(d10x))
  #colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),samples[y],sep="-")
  print(dim(d10x))
  
  # Load data and estimate soup profile
  # Estimate rho
  d10x_soup = autoEstCont(d10x)
  print(dim(d10x_soup))
  print(class(d10x_soup))
  # Clean the data
  out = adjustCounts(d10x_soup)
  print(dim(out))
  print(clas(out))
  
  #' Initialize the Seurat object with the raw (non-normalized data).
  d10x<-CreateSeuratObject(counts = out, project = "ped_adult_map", min.cells = 3, min.features = 0)

  #add meta data to each seurat object
  meta_cell<-data.frame(cell=colnames(d10x), individual=sapply(colnames(d10x), function(x) strsplit(x,"-")[[1]][2]))
  # print(head(meta_cell))
  meta_cell_add<-merge(meta_cell, meta, by.x="individual", by.y="Sample_ID")
  meta_cell_add<-meta_cell_add[match(colnames(d10x), meta_cell_add$cell),]
  # print(identical(meta_cell_add$cell, colnames(d10x)))
  rownames(meta_cell_add)<-meta_cell_add$cell
  d10x<- AddMetaData(d10x, meta_cell_add)
  d10x




