#'---
#'title: scRNAseq Differential Expression
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---



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
#library(ggsignif)


options(stringsAsFactors = FALSE)

source("R_functions/pretty_plots.R")



## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","d10x_adult_ped_raw.rds"))


######
## add cell type labels from split analysis
######
#immune
load(here("output","immune_iterative_label.Rdata"))
#stromal
load(here("output","strom_celltype_label.Rdata"))
#epithelial
load(here("output","epi_celltype_label.Rdata"))

cell_label<-rbind(epi_cell_labels, immune_cell_labels, stromal_cell_labels)
cell_label$cluster_ID<-as.character(cell_label$cluster_ID)
cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

## remove doublets
d10x<-subset(d10x, subset = cluster_ID != "Doublet")
d10x<-subset(d10x, subset = cluster_ID != "Paneth (UC only)")
d10x.cd.ctrl<-subset(d10x, subset = orig.ident %in% c("CD","Ctrl"))

## testing factor
d10x.cd.ctrl$cell_stim<-paste(d10x.cd.ctrl$cluster_ID, d10x.cd.ctrl$orig.ident, sep = "_")
Idents(d10x.cd.ctrl) <- "cell_stim"

table(d10x.cd.ctrl$cluster_ID, d10x.cd.ctrl$orig.ident)



#MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene,
#consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
#and linear regression for the continuous process (i.e., the expression level). 

cell_types<-unique(d10x.cd.ctrl$cluster_ID)

contrasts_celltype_stim<-do.call(rbind,lapply(1:length(cell_types), function(x){
  combinations(n = 2, r = 2, v = d10x.cd.ctrl$cell_stim[grep(cell_types[x],d10x.cd.ctrl$cell_stim)], repeats.allowed = FALSE)}))

contrasts_celltype_stim

nrow(contrasts_celltype_stim)

contrasts_celltype_stim[37,]<-c("enterocyte_CD","enterocyte_Ctrl")
contrasts_celltype_stim

d10x.cd.ctrl
# this is 922,965 tests across all comparisons (24,945 genes, 37 comparisons)

diff_exp_all<-lapply(1:nrow(contrasts_celltype_stim), function(x){
  de<-FindMarkers(d10x.cd.ctrl, ident.1 = contrasts_celltype_stim[x,1], ident.2 = contrasts_celltype_stim[x,2], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
  print(paste(contrasts_celltype_stim[x,1],"vs", contrasts_celltype_stim[x,2],":", nrow(de), sep=" "))
  de$gene<-rownames(de)
  rownames(de)<-NULL
  de<-de[,c(6,1:5)]
  de$cell.1<-contrasts_celltype_stim[x,1]
  de$cell.2<-contrasts_celltype_stim[x,2]
  de})


diff_exp_all<-do.call(rbind, diff_exp_all)

save(diff_exp_all, file=here("data","primary_diff_genes.RData"))
#load(file=here("../../../codon/scRNAseq_codon/data","primary_diff_genes.RData"))
