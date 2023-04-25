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
library(anndata)

source("scripts/00_pretty_plots.R")
source("scripts/00_entropy_d10x.R")

#################
## Load raw QC'ed data
#################
d10x_fetal<-readRDS(here("data","d10x_fetal_raw.rds"))
#d10x_fetal<-readRDS("/media/redgar/Seagate Portable Drive/fetal_liver/d10x_fetal_raw.rds")

d10x_ped_IFALD<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))
## add cell type labels
load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x_ped_IFALD), cell_label$index),]
identical(colnames(d10x_ped_IFALD), cell_label$index)
d10x_ped_IFALD <- AddMetaData(d10x_ped_IFALD, metadata = cell_label)


#################
## Match meta data
#################
d10x_fetal@meta.data$Barcodes<-NULL
colnames(d10x_fetal@meta.data)[which(colnames(d10x_fetal@meta.data)=="Extract.Name")]<-"Sample"
colnames(d10x_fetal@meta.data)[which(colnames(d10x_fetal@meta.data)=="Characteristics.individual.")]<-"individual"
d10x_fetal@meta.data$Treatment<-"Healthy"
d10x_fetal@meta.data$Tissue<-"TLH"
d10x_fetal@meta.data$chemistry<-"3pr"
colnames(d10x_fetal@meta.data)[which(colnames(d10x_fetal@meta.data)=="Characteristics.sex.")]<-"Sex"
colnames(d10x_fetal@meta.data)[which(colnames(d10x_fetal@meta.data)=="Characteristics.age.")]<-"Age"
d10x_fetal@meta.data$AgeGroup<-"Fetal"
d10x_fetal@meta.data$FreshorFrozen<-"fresh"
d10x_fetal@meta.data$BMI<-NA
d10x_fetal@meta.data$relALBChange<-NA
d10x_fetal@meta.data$nuclear_fraction<-NA
d10x_fetal@meta.data$cell_status<-NA
colnames(d10x_fetal@meta.data)[which(colnames(d10x_fetal@meta.data)=="Cell.Labels")]<-"CellType_refined"
d10x_fetal@meta.data$age_condition<-paste(d10x_fetal$AgeGroup, d10x.combined$Treatment, sep=" ")


d10x_ped_IFALD@meta.data$file<-NULL        
d10x_ped_IFALD@meta.data$Approx_bam_GB<-NULL   
d10x_ped_IFALD@meta.data$Characteristics.facs.sorting.<-NA
d10x_ped_IFALD@meta.data$Sample<-d10x_ped_IFALD@meta.data$individual
d10x_ped_IFALD@meta.data$individual<-sapply(1:nrow(d10x_ped_IFALD@meta.data), function(x) strsplit(d10x_ped_IFALD@meta.data$individual[x],"_")[[1]][1])
d10x_ped_IFALD@meta.data$CellType_rough <-NULL    
d10x_ped_IFALD@meta.data$second_best_cell  <-NULL    
d10x_ped_IFALD@meta.data$S.Score   <-NULL 
d10x_ped_IFALD@meta.data$G2M.Score <-NULL 
d10x_ped_IFALD@meta.data$Phase     <-NULL 
d10x_ped_IFALD@meta.data$old.ident<-NULL 
d10x_ped_IFALD@meta.data$integrated_snn_res.0.5 <-NULL 
d10x_ped_IFALD@meta.data$seurat_clusters<-NULL 
d10x_ped_IFALD@meta.data$age_id<-NULL 
d10x_ped_IFALD@meta.data$index<-NULL 


head(d10x_ped_IFALD@meta.data)
head(d10x_fetal@meta.data)

d10x_fetal@meta.data<-d10x_fetal@meta.data[,colnames(d10x_ped_IFALD@meta.data)]

head(d10x_fetal@meta.data)
head(d10x_ped_IFALD@meta.data)






#################
## Merge
#################
d10x <- merge(d10x_ped_IFALD,d10x_fetal, merge.data=TRUE, project = "IFALD_fetal_adult_ped_map")

###############
## Integrate by donor
###############
#https://satijalab.org/seurat/articles/integration_rpca.html
print("RUNNING INTEGRATION")

## run integration across donor and hopefully that will also smooth out differences with chemistry and batch?
d10x.list<- SplitObject(d10x, split.by = "sample")

# normalize, identify variable features and score cell cycle for each dataset independently
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

d10x.list <- lapply(X = d10x.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = d10x.list)
d10x.list <- lapply(X = d10x.list, FUN = function(x) {
  #x <- ScaleData(x, features = features, verbose = FALSE)
  x <- ScaleData(x, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

## Identify anchors
chem.anchors <- FindIntegrationAnchors(object.list = d10x.list, anchor.features = features, reduction = "rpca")
d10x.fetal_ped_IFALD <- IntegrateData(anchorset = chem.anchors)

DefaultAssay(d10x.fetal_ped_IFALD) <- "integrated"

print("INTEGRATED")


# Run the standard workflow for visualization and clustering
d10x.fetal_ped_IFALD <- ScaleData(d10x.fetal_ped_IFALD, verbose = FALSE)
d10x.fetal_ped_IFALD <- RunPCA(d10x.fetal_ped_IFALD, npcs = 30, verbose = FALSE)
d10x.fetal_ped_IFALD <- RunUMAP(d10x.fetal_ped_IFALD, reduction = "pca", dims = 1:30)
d10x.fetal_ped_IFALD <- RunTSNE(d10x.fetal_ped_IFALD, dims = 1:30)

d10x.fetal_ped_IFALD <- FindNeighbors(d10x.fetal_ped_IFALD, reduction = "pca", dims = 1:30)
d10x.fetal_ped_IFALD <- FindClusters(d10x.fetal_ped_IFALD, resolution = 0.5)

d10x.fetal_ped_IFALD


##############
## Save integrated to look local
##############
save(d10x.fetal_ped_IFALD, file=paste(here("data/"),"Fetal_IFALD_adult_ped_integrated.rds", sep=""))
cell_label<-d10x.fetal_ped_IFALD@meta.data
save(cell_label, file=paste(here("data/"),"Fetal_IFALD_adult_ped_cellRough.rds", sep=""))

