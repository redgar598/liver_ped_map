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

#d10x<-readRDS(here("/media/redgar/Seagate Portable Drive/spinalcord_tutorial_data/d10x_SCI_merged.rds"))
d10x<-readRDS(here("data/d10x_SCI_merged.rds"))
d10x78<-readRDS(here("data/d10x_78_SCI_merged.rds"))

d10x<- merge(d10x, y= d10x78, merge.data=TRUE, project = "SCI")
d10x

d10x.list.sample <- SplitObject(d10x, split.by = "sample_ID")

# normalize, identify variable features and score cell cycle for each dataset independently
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

d10x.list.sample <- lapply(X = d10x.list.sample, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = d10x.list.sample)
d10x.list.sample <- lapply(X = d10x.list.sample, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

## Identify anchors
sample.anchors <- FindIntegrationAnchors(object.list = d10x.list.sample, anchor.features = features, reduction = "rpca")
d10x.combined <- IntegrateData(anchorset = sample.anchors)

DefaultAssay(d10x.combined) <- "integrated"

print("INTEGRATED")

# Run the standard workflow for visualization and clustering
d10x.combined <- ScaleData(d10x.combined, verbose = FALSE)
d10x.combined <- RunPCA(d10x.combined, npcs = 30, verbose = FALSE)
d10x.combined <- RunUMAP(d10x.combined, reduction = "pca", dims = 1:30)
d10x.combined <- RunTSNE(d10x.combined, dims = 1:30)

d10x.combined <- FindNeighbors(d10x.combined, reduction = "pca", dims = 1:30)
d10x.combined <- FindClusters(d10x.combined, resolution = 0.5)

d10x.combined

save(d10x.combined, file=here("data/d10x_SCI_fullintegrated.rds"))