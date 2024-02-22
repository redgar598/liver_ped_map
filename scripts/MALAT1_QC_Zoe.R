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
library(SCINA)




source("scripts/00_pretty_plots.R")
source("scripts/00_entropy_d10x.R")
source("scripts/00_fanciest_UMAP.R")


#cellranger_out_loc <- here("/media/redgar/Seagate Portable Drive/ped_liver_map_raw")
#cellranger_out_realigned <- here("/media/redgar/Seagate Portable Drive/realign_samples")
cellranger_out_loc <- here("../../../projects/macparland/RE/PediatricAdult/ped_liver_map_raw")
cellranger_out_realigned <- here("../../../projects/macparland/RE/PediatricAdult/realign_samples")

samples<-c(list.files(cellranger_out_loc),list.files(cellranger_out_realigned))
print(samples)

meta<-read.table(here("data/data_transfer_updated_feb12_2024_IFALD_PBMC.csv"), header=T, sep=",")

samples<-samples[which(samples%in%meta$file)]
samples<-samples[-grep("PBMC",samples)]


d10x.list <- sapply(1:length(samples), function(y){
  if(length(grep("realign",samples[y]))==1){
    dataset_loc<-cellranger_out_realigned
  }else{
    dataset_loc<-cellranger_out_loc
  }

  caud<-meta$Sample_ID[which(meta$file == samples[y])]
  print(caud)
  print(file.path(dataset_loc,paste(samples[y],"/outs", sep=""),"filtered_feature_bc_matrix"))
  d10x <- Read10X(file.path(dataset_loc,paste(samples[y],"/outs", sep=""),"filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),caud,sep="-")
  # print(dim(d10x))
  #' Initialize the Seurat object with the raw (non-normalized data).
  d10x<-CreateSeuratObject(counts = d10x, project = "ped_adult_map", min.cells = 0, min.features = 0)

  #add meta data to each seurat object
  meta_cell<-data.frame(cell=colnames(d10x), individual=caud)
  meta_cell_add<-merge(meta_cell, meta, by.x="individual", by.y="Sample_ID")
  meta_cell_add$relALBChange<-NA
  meta_cell_add<-meta_cell_add[match(colnames(d10x), meta_cell_add$cell),]
  print(identical(meta_cell_add$cell, colnames(d10x)))
  rownames(meta_cell_add)<-meta_cell_add$cell
  d10x<- AddMetaData(d10x, meta_cell_add)
  
  d10x <- NormalizeData(d10x)

})

d10x.list



#######################################################################################################################################

d10x.list<-lapply(1:length(d10x.list), function(x){
  sobj<-d10x.list[[x]]
  
  # This script filters cells based on MALAT1 expression
  # Code is written assuming the use of a Seurat object, but
  # should be able to be applied to any single-cell object
  
  # Look at histogram of MALAT1 counts
  hist(sobj@assays$RNA@data["MALAT1",], freq = FALSE, breaks=100)
  
  # If there is a peak (even a small peak) at zero, followed by a dip, then a more
  # normal distribution, run the following code.
  
  # Choose rough x value that surrounds the minimum between 0 and where slope starts to rise again; this is max_counts
  # Counts is vector of MALAT1 counts, min and max are range of count values to look between
  # Lower bw if minimum doesn't look totally accurate
  # Change lwd to change thickness of red line
  # Change breaks to change bins of histogram
  define_malat1_threshold <- function(counts, max_counts, min_counts = 0, bw = 0.1, lwd = 2, breaks = 100) {
    pdf(here("figures",paste(unique(sobj$individual),"MALAT1.pdf", sep="")))
    # Visualise MALAT1 histogram
    hist(counts, breaks = breaks, freq = FALSE, main=unique(sobj$individual))
    # Calculate the density values
    density_data <- density(counts, bw = bw)
    # Isolate x value at minimum
    
    ######## Updated this line as in one sample the first values in counts was negative so threw off the indexing
    threshold <- density_data$x[density_data$x >= min_counts & density_data$x <= max_counts][which(density_data$y[density_data$x >= min_counts & density_data$x <= max_counts] == min(density_data$y[density_data$x >= min_counts & density_data$x <= max_counts]))]
   
    # Visualise on histogram
    abline(v = threshold, col = "red", lwd = lwd)
    dev.off() 
    return(threshold)
  }
  
  # Run this function on MALAT1 reads, eg:
  threshold <- define_malat1_threshold(sobj@assays$RNA@data["MALAT1",], max_counts = 2)
  
  # Use this value to subset your cells
  cells_subset <- WhichCells(sobj, expression = MALAT1 > threshold, slot = "counts")
  sobj_subset <- subset(sobj, cells = cells_subset)

  print(paste(unique(sobj$individual),":",ncol(sobj), "cell originally", ncol(sobj_subset), "after MALAT1 filter"))

  sobj_subset

})


#######################################################################################################################################


## cell counts
plt_count_raw<-lapply(1:length(d10x.list), function(x) {
  df<-data.frame(raw_cell_count=nrow(d10x.list[[x]]@meta.data),individual=unique(d10x.list[[x]]@meta.data$individual))
  df})
plt_count_raw<-do.call(rbind, plt_count_raw)
print(plt_count_raw)



#'## QC
#'The percentage of reads that map to the mitochondrial genome
#'Low-quality / dying cells often exhibit extensive mitochondrial contamination
#'We calculate mitochondrial QC metrics with the PercentageFeatureSet function, which calculates the percentage of counts originating from a set of features
#'We use the set of all genes starting with MT- as a set of mitochondrial genes

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
invisible(lapply(1:length(d10x.list), function(x){
  d10x.list[[x]][["percent.mt"]] <<- PercentageFeatureSet(d10x.list[[x]], pattern = "^MT-")}))

# Show QC metrics for the first 5 cells
print(head(d10x.list[[2]]@meta.data, 5))

#'Low-quality cells or empty droplets will often have very few genes
#'Cell doublets or multiplets may exhibit an aberrantly high gene count


# Visualize QC metrics
#nFeature number of unique genes
#nCount number of total molecules
plt_QC_data<-do.call(rbind, lapply(1:length(d10x.list), function(x) d10x.list[[x]]@meta.data))
save(plt_QC_data, file=here("data","MALAT1_QC_metrics.Rdata"))
