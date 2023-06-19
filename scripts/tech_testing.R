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
source("scripts/00_plot_gene_exp.R")


d10x_liver<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))

tech1<-subset(d10x_liver, subset = individual == "C39_caud3prTLH")
tech1<-tech1[, sample(colnames(tech1), size = 500, replace=F)]

tech2<-subset(d10x_liver, subset = individual == "C39_caud3prNPC")
tech2<-tech2[, sample(colnames(tech2), size = 500, replace=F)]

query1<-subset(d10x_liver, subset = individual == "C93_caud3pr")
query1<-query1[, sample(colnames(query1), size = 500, replace=F)]

query2<-subset(d10x_liver, subset = individual == "C97_caud3pr")
query2<-query2[, sample(colnames(query2), size = 500, replace=F)]


###############
## Tech integration
###############
d10x.list<- list(tech1 = tech1, tech2 = tech2)

## run integration across donor and hopefully that will also smooth out differences with chemistry?
## data is already split by donor

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
  x <- ScaleData(x, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


## Identify anchors
tech.anchors <- FindIntegrationAnchors(object.list = d10x.list, anchor.features = features, reduction = "rpca")
d10x.combined <- IntegrateData(anchorset = tech.anchors)

tech.anchors.callnames<-tech.anchors@anchors
tech.anchors.callnames$cellname1<-tech1$cell[tech.anchors.callnames$cell1]
tech.anchors.callnames$cellname2<-tech2$cell[tech.anchors.callnames$cell2]


head(chem.anchors@anchors)
tech1$highlight_anchors<-NA
tech1$highlight_anchors[chem.anchors@anchors$cell1]<-"anchor1"
tech2$highlight_anchors<-NA
tech2$highlight_anchors[chem.anchors@anchors$cell2]<-"anchor2"

d10x.merged<-merge(tech1, tech2)

## Add cell type
load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x.merged), cell_label$index),]
identical(colnames(d10x.merged), cell_label$index)

d10x.merged <- AddMetaData(d10x.merged, metadata = cell_label)

d10x.merged <- NormalizeData(d10x.merged)
d10x.merged <- FindVariableFeatures(d10x.merged, selection.method = "vst", nfeatures = 2000)
d10x.merged <- ScaleData(d10x.merged) 
d10x.merged <- RunPCA(d10x.merged, npcs = 30, verbose = FALSE)
d10x.merged <- RunUMAP(d10x.merged, reduction = "pca", dims = 1:30)



DimPlot(d10x.merged, group.by = "CellType_refined")+colscale_cellType
DimPlot(d10x.merged, group.by = "highlight_anchors")



##############
##### integrate all samples with only tech anchors
##############
rm(d10x.combined)

d10x.merge1<-merge(tech1, query1)
d10x.merge2<-merge(tech2, query2)

d10x.list<- list(tech1 = d10x.merge1, tech2 = d10x.merge2)

## run integration across donor and hopefully that will also smooth out differences with chemistry?
## data is already split by donor

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
  x <- ScaleData(x, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


head(tech.anchors@object.list)
head(tech.anchors@anchors)


## Identify anchors
query.anchors <- FindIntegrationAnchors(object.list = d10x.list, anchor.features = features, reduction = "rpca")

query_techanchors<-query.anchors
query_techanchors@anchors<-tech.anchors@anchors

# replace anchors with tech anchors, update cell index to full data set though
cell_match1<-data.frame(cellname1 = as.character(d10x.list[[1]]$cell[which(d10x.list[[1]]$cell%in%tech.anchors.callnames$cellname1)]),
           index1 = which(d10x.list[[1]]$cell%in%tech.anchors.callnames$cellname1))
cell_match<-merge(cell_match1, tech.anchors.callnames, by.x="cellname1", by.y="cellname1")

cell_match2<-data.frame(cellname2 = as.character(d10x.list[[2]]$cell[which(d10x.list[[2]]$cell%in%tech.anchors.callnames$cellname2)]),
                       index2 = which(d10x.list[[2]]$cell%in%tech.anchors.callnames$cellname2))
cell_match<-merge(cell_match2, cell_match, by.x="cellname2", by.y="cellname2")

cell_match<-cell_match[,c("cell1", "cell2","score", "dataset1", "dataset2")]
cell_match<-cell_match[order(cell_match$dataset1, cell_match$cell1),]

query_techanchors@anchors<-cell_match

d10x.combined <- IntegrateData(anchorset = query_techanchors)

d10x.combined <- FindVariableFeatures(d10x.combined, selection.method = "vst", nfeatures = 2000)
d10x.combined <- ScaleData(d10x.combined) 
d10x.combined <- RunPCA(d10x.combined, npcs = 30, verbose = FALSE)
d10x.combined <- RunUMAP(d10x.combined, reduction = "pca", dims = 1:30)

DimPlot(d10x.combined)


## check other anchors produces a different result
rm(d10x.combined)
query.anchors <- FindIntegrationAnchors(object.list = d10x.list, anchor.features = features, reduction = "rpca")
d10x.combined <- IntegrateData(anchorset = query.anchors)

d10x.combined <- FindVariableFeatures(d10x.combined, selection.method = "vst", nfeatures = 2000)
d10x.combined <- ScaleData(d10x.combined) 
d10x.combined <- RunPCA(d10x.combined, npcs = 30, verbose = FALSE)
d10x.combined <- RunUMAP(d10x.combined, reduction = "pca", dims = 1:30)

DimPlot(d10x.combined)



###########
## highlight anchors
###########
head(chem.anchors@anchors)
d10x.merge1$highlight_anchors<-NA
d10x.merge1$highlight_anchors[chem.anchors@anchors$cell1]<-"anchor1"
d10x.merge2$highlight_anchors<-NA
d10x.merge2$highlight_anchors[chem.anchors@anchors$cell2]<-"anchor2"

d10x.merged<-merge(d10x.merge1, d10x.merge2)

## Add cell type
load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x.merged), cell_label$index),]
identical(colnames(d10x.merged), cell_label$index)

d10x.merged <- AddMetaData(d10x.merged, metadata = cell_label)

d10x.merged <- NormalizeData(d10x.merged)
d10x.merged <- FindVariableFeatures(d10x.merged, selection.method = "vst", nfeatures = 2000)
d10x.merged <- ScaleData(d10x.merged) 
d10x.merged <- RunPCA(d10x.merged, npcs = 30, verbose = FALSE)
d10x.merged <- RunUMAP(d10x.merged, reduction = "pca", dims = 1:30)



DimPlot(d10x.merged, group.by = "CellType_refined")+colscale_cellType
DimPlot(d10x.merged, group.by = "highlight_anchors")


