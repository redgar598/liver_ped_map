#'### Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(here)
library(ggplot2)
library(reshape2)
library(gridExtra)
#library(limma)
library(cowplot)
library(gtools)
library(ggsignif)
library(stats)
library(scales)
library(colorspace)



options(stringsAsFactors = FALSE)

source("scripts/00_pretty_plots.R")



## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","d10x_adult_ped_raw.rds"))


######
## add cell type labels
######
load(here("data","adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)




##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")
Idents(d10x) <- "AgeGroup"
cell_types<-unique(d10x$CellType_rough)


#########
## Signature Genes
#########
myeloid_immune_supressive<-c("CTSB","CD163","MS4A7","FOLR2","GPNMB","VSIG4","HMOX1","MSR1")
inflammatory_macs<-c("CD74","HLA-DRA","TYROBP","C1QC","HLA-DPA1","HLA-DPB1","LYZ","S100A6")
exhausted_tcells<-c("TOX","PDCD1","LAG3","TNFRSF9","CXCL13","ENTPD1","HAVCR2","CD38")

recent_recruit_myeloid<-c("S100A8","S100A9","CD68","LYZ")
kuffer_signature<-c("VSIG4","MARCO","CD5L","HMOX1")
#kuffer_signature<-c("CD74","HLA-DRA","TYROBP","C1QC","HLA-DPA1","HLA-DPB1","S100A6","MARCO","CD5L")



######
## Score Signatures
######
d10x <- AddModuleScore(
  object = d10x,
  features = list(myeloid_immune_supressive),
  ctrl = 5,
  name = 'myeloid_immune_supressive_score'
)

d10x <- AddModuleScore(
  object = d10x,
  features = list(inflammatory_macs),
  ctrl = 5,
  name = 'inflammatory_macs_score'
)

d10x <- AddModuleScore(
  object = d10x,
  features = list(exhausted_tcells),
  ctrl = 5,
  name = 'exhausted_tcells_score'
)

d10x <- AddModuleScore(
  object = d10x,
  features = list(recent_recruit_myeloid),
  ctrl = 5,
  name = 'recently_recruited_myeloid'
)

d10x <- AddModuleScore(
  object = d10x,
  features = list(kuffer_signature),
  ctrl = 5,
  name = 'kuffer_like_score'
)


score_data<-d10x@meta.data[,c("myeloid_immune_supressive_score1","inflammatory_macs_score1","exhausted_tcells_score1","recently_recruited_myeloid1","kuffer_like_score1")]
rm(d10x)
gc()

### load integrate for UMAP etc
load(here("data","adult_ped_integrated_refinedlabels_withDropletQC.rds"))
score_data<-score_data[match(rownames(d10x.combined@meta.data),rownames(score_data)),]
identical(rownames(d10x.combined@meta.data),rownames(score_data))
d10x.combined <- AddMetaData(d10x.combined, metadata = score_data)




#######
## just myeloid cells
#######
d10x.combined_KC<-subset(d10x.combined, subset = CellType_refined %in% c("KC Like","RR Myeloid"))
d10x.combined_KC<-subset(d10x.combined_KC, subset = individual %in% c("C70_caud5pr","C64_caud5pr","C85_caud3pr","C93_caud3pr", "C96_caud3pr", "C82_caud3pr"))

rm(d10x.combined)
gc()

d10x.combined_KC <- RunPCA(d10x.combined_KC, npcs = 30, verbose = FALSE)
d10x.combined_KC <- RunUMAP(d10x.combined_KC, reduction = "pca", dims = 1:30)
d10x.combined_KC <- FindNeighbors(d10x.combined_KC, reduction = "pca", dims = 1:30)
d10x.combined_KC <- FindClusters(d10x.combined_KC, resolution = 0.1)
DimPlot(d10x.combined_KC, label=T)
DimPlot(d10x.combined_KC, label=T, reduction="pca", group.by = "AgeGroup")

FeaturePlot(d10x.combined_KC, features = c("CD86","MRC1","ALB"),reduction="pca", min.cutoff = "q9", pt.size=1)
FeaturePlot(d10x.combined_KC, features = c("recently_recruited_myeloid1","kuffer_like_score1","ALB"),reduction="pca", min.cutoff = "q9", pt.size=1)

ggplot(d10x.combined_KC@meta.data, aes(AgeGroup, kuffer_like_score1))+geom_violin()+geom_boxplot()+geom_point()
ggplot(d10x.combined_KC@meta.data, aes(recently_recruited_myeloid1, kuffer_like_score1))+geom_point()
ggplot(d10x.combined_KC@meta.data, aes(recently_recruited_myeloid1, kuffer_like_score1))+geom_point()+facet_wrap(~AgeGroup)
