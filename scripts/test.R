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
library(colorspace)
library(cowplot)
library(DropletQC)




source("scripts/00_pretty_plots.R")
source("scripts/00_entropy_d10x.R")

#dataset_loc <- here("/media/redgar/Seagate Portable Drive/ped_liver_map_raw")
dataset_loc <- here("../../../projects/macparland/RE/PediatricAdult/ped_liver_map_raw")

samples<-list.files(dataset_loc)
samples<-samples[-grep("data_transfer",samples)]
print(samples)

#meta<-read.table(here("data/data_transfer_updated_jan16_2023.csv"), header=T, sep=",")
meta<-read.table(here(dataset_loc,"data_transfer_updated_jan16_2023.csv"), header=T, sep=",")


y=1

  caud<-meta$Sample_ID[which(meta$file == samples[y])]
  print(caud)
  print(file.path(dataset_loc,paste(samples[y],"/outs", sep=""),"filtered_feature_bc_matrix"))
  d10x <- Read10X(file.path(dataset_loc,paste(samples[y],"/outs", sep=""),"filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),caud,sep="-")
  # print(dim(d10x))
  #' Initialize the Seurat object with the raw (non-normalized data).
  d10x<-CreateSeuratObject(counts = d10x, project = "ped_adult_map", min.cells = 0, min.features = 0)
  
  ## SoupX needs clusters so quickly make clusters for each sample
  d10x    <- SCTransform(d10x, verbose = F)
  d10x    <- RunPCA(d10x, verbose = F)
  d10x    <- RunUMAP(d10x, dims = 1:30, verbose = F)
  d10x    <- FindNeighbors(d10x, dims = 1:30, verbose = F)
  d10x    <- FindClusters(d10x, verbose = T)
  meta_clusters    <- d10x@meta.data
  
  sc = load10X(file.path(dataset_loc,paste(samples[y],"/outs", sep="")))
  sc = setClusters(sc, setNames(meta_clusters$seurat_clusters, rownames(meta_clusters)))
  
  ######
  ## Load data and estimate soup profile
  ######
  # Estimate rho
  sc = autoEstCont(sc)
  #Genes with highest expression in background. These are often enriched for ribosomal proteins.
  print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20))
  # Clean the data
  sc = adjustCounts(sc)
  
  d10x = CreateSeuratObject(sc)
  
  
  #add meta data to each seurat object
  meta_cell<-data.frame(cell=colnames(d10x), individual=caud)
  meta_cell_add<-merge(meta_cell, meta, by.x="individual", by.y="Sample_ID")
  meta_cell_add<-meta_cell_add[match(colnames(d10x), meta_cell_add$cell),]
  print(identical(meta_cell_add$cell, colnames(d10x)))
  rownames(meta_cell_add)<-meta_cell_add$cell
  d10x<- AddMetaData(d10x, meta_cell_add)
  
  
  
  ######
  ## dropletQC
  ######
  nf1 <- nuclear_fraction_tags(
    outs = file.path(dataset_loc,paste(samples[y],"/outs", sep="")),
    tiles = 1, cores = 1, verbose = FALSE)
  head(nf1)
  
  print(identical(rownames(nf1), colnames(d10x)))
  d10x<- AddMetaData(d10x, nf1)
  d10x
  
  nf.umi <- data.frame(nf=d10x$nuclear_fraction,
                       umi=d10x$nCount_RNA)
  
  # Run identify_empty_drops
  empty_drop <- identify_empty_drops(nf_umi=nf.umi)
  empty_drop$individual<-d10x$individual
  empty_drop_damagedcell <- identify_damaged_cells(empty_drop, verbose = FALSE, output_plots = F)
  
  head(empty_drop_damagedcell[[1]])
  table(empty_drop_damagedcell[[1]]$cell_status)
  
  print(identical(rownames(empty_drop_damagedcell[[1]]), colnames(d10x)))
  d10x<- AddMetaData(d10x, empty_drop_damagedcell[[1]])
  d10x$nf<-NULL
  d10x$umi<-NULL
  d10x
print(head(d10x))