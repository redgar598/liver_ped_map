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

source("scripts/00_pretty_plots.R")



## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","d10x_adult_ped_raw.rds"))


######
## add cell type labels
######
load(here("data","adult_ped_cellRefined.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)




##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")


## testing factor
levels(d10x$CellType_refined)[which(levels(d10x$CellType_refined)=="LSEC\n(Hepatocyte Like)")]<-"LSEC_hep"
levels(d10x$CellType_refined)[which(levels(d10x$CellType_refined)=="LSEC")]<-"LSEC_nothep"
levels(d10x$CellType_refined)[which(levels(d10x$CellType_refined)=="Neutrophil\n(DEFA+)")]<-"Neutrophil_DEFA"
levels(d10x$CellType_refined)[which(levels(d10x$CellType_refined)=="Neutrophil")]<-"Neutrophil_notDEFA"

d10x$cell_sex<-paste(as.character(d10x$CellType_refined), d10x$Sex, sep = "_")
Idents(d10x) <- "cell_sex"

table(d10x$CellType_refined, d10x$Sex)


#MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene,
#consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
#and linear regression for the continuous process (i.e., the expression level). 

cell_types<-unique(as.character(d10x$CellType_refined))
cell_types<-cell_types[-grep("Hepatocyte Like",cell_types)]
cell_types<-cell_types[-grep("Neutrophil_DEFA",cell_types)]# on 11 male defa neutrophils so not meaningful to compare
cell_types[grep("CD3",cell_types)]<-"CD3"

contrasts_celltype_sex<-do.call(rbind,lapply(1:length(cell_types), function(x){
  combinations(n = 2, r = 2, v = d10x$cell_sex[grep(cell_types[x],d10x$cell_sex)], repeats.allowed = FALSE)
}))

contrasts_celltype_sex

nrow(contrasts_celltype_sex)


###########
## Monte carlo the DE
###########

d10x_F<-subset(d10x, subset = Sex == "F")
ncol(d10x_F)
d10x_M<-subset(d10x, subset = Sex == "M")
ncol(d10x_M)


### paralize
command_args <- commandArgs(trailingOnly = TRUE)
cell_type_indx <- as.numeric(command_args[1])
cell_type<-cell_types[cell_type_indx]

samp_num=1000


#DE_monte_carlo<-lapply(cell_types, function(cell_type){

contrasts_celltype<-contrasts_celltype_sex[grep(cell_type, contrasts_celltype_sex)]
cell_type<-as.character(unique(d10x$CellType_refined)[grep(cell_type,unique(d10x$CellType_refined))])

d10x_F_celltype<-subset(d10x_F, subset = CellType_refined == cell_type)
ncol(d10x_F_celltype)
d10x_M_celltype<-subset(d10x_M, subset = CellType_refined == cell_type)
ncol(d10x_M_celltype)


de_lists<-sapply(1:samp_num, function(x){
  set.seed(x)
  
  ## make downsampled
  if(ncol(d10x_F_celltype)<ncol(d10x_M_celltype)){
    M_cells_random <- d10x_M_celltype[, sample(colnames(d10x_M_celltype), size = ncol(d10x_F_celltype), replace=F)]
    d10_DE<-merge(d10x_F_celltype, M_cells_random)
  }else{
    F_cells_random <- d10x_F_celltype[, sample(colnames(d10x_F_celltype), size = ncol(d10x_M_celltype), replace=F)]
    d10_DE<-merge(d10x_M_celltype, F_cells_random)}
  
  ## run DE 
  de<-FindMarkers(d10_DE, ident.1 = contrasts_celltype[1], ident.2 = contrasts_celltype[2], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
  print(paste(contrasts_celltype[1],"vs", contrasts_celltype[2],":", nrow(de), sep=" "))
  de$gene<-rownames(de)
  rownames(de)<-NULL
  de<-de[,c(6,1:5)]
  de$cell.1<-contrasts_celltype[1]
  de$cell.2<-contrasts_celltype[2]
  
  de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]$gene
})

sig_gene_count<-unlist(de_lists)
if(length(sig_gene_count)==0){NA}else{
  sig_gene_count<-as.data.frame(table(sig_gene_count))
  colnames(sig_gene_count)<-c("gene","sig_count")
  
  sig_gene_count$monte_carlo_sig<-sapply(1:nrow(sig_gene_count),function(x){
    1-(sig_gene_count$sig_count[x]+1)/(samp_num+1)
  })
  
  sig_gene_count$cell<-cell_type
  sig_gene_count}#})

#DE_monte_carlo<-do.call(rbind, DE_monte_carlo)
DE_monte_carlo<-sig_gene_count
head(DE_monte_carlo)
DE_monte_carlo<-DE_monte_carlo[which(!(is.na(DE_monte_carlo$gene))),]

save(DE_monte_carlo, file=here("data",paste(cell_type,"sex_diff_motecarlo_1000.RData",sep="_")))

