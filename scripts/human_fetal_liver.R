### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(reshape2)
library(gtools)
library(colorspace)
library(cowplot)
library(RColorBrewer)


#################
## Load raw QC'ed data
#################
d10x_fetal<-readRDS(here("../../../projects/macparland/RE/PediatricAdult/processed_data","d10x_fetal_raw.rds"))
#d10x_fetal<-readRDS("/media/redgar/Seagate Portable Drive/fetal_liver/d10x_fetal_raw.rds")

Idents(d10x_fetal) <- "Cell.Labels"
DefaultAssay(d10x_fetal)<-"RNA"
all.markers<-FindAllMarkers(d10x_fetal)

all.markers[grep("MTM1", all.markers$gene),]

save(all.markers,  file="../../../projects/macparland/RE/random_side/Fetal_differential_MTM1.RData")

Mtm1_sig<-all.markers[grep("MTM1", all.markers$gene),][which(all.markers[grep("Mtm1", all.markers$gene),]$p_val_adj<0.005),]
Mtm1_sig

d10x_fetal <- NormalizeData(d10x_fetal)
d10x_fetal <- FindVariableFeatures(d10x_fetal, selection.method = "vst", nfeatures = 2000)
d10x_fetal <- ScaleData(d10x_fetal) 
d10x_fetal <- RunPCA(d10x_fetal, ndims.print = 1:10, nfeatures.print = 10)
d10x_fetal <- RunUMAP(d10x_fetal, dims = 1:30)

umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x_fetal, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x_fetal@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")



DefaultAssay(d10x_fetal)<-"RNA"
gene="MTM1"
gene_exp<-FetchData(d10x_fetal, vars=gene)
gene_exp$cell<-rownames(gene_exp)
plt_myeloid<-merge(plt_myeloid, gene_exp, by='cell')

save(plt_myeloid, file="../../../projects/macparland/RE/random_side/Fetal_UMAP_mtm1.RData")
