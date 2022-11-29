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
## add cell type labels
######
load(here("data","adult_ped_cellRough.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)




##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")


## testing factor
d10x$cell_age<-paste(d10x$CellType_rough, d10x$AgeGroup, sep = "_")
Idents(d10x) <- "cell_age"

table(d10x$CellType_rough, d10x$AgeGroup)


#MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene,
#consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
#and linear regression for the continuous process (i.e., the expression level). 

cell_types<-unique(d10x$CellType_rough)

contrasts_celltype_age<-do.call(rbind,lapply(1:length(cell_types), function(x){
  combinations(n = 2, r = 2, v = d10x$cell_age[grep(cell_types[x],d10x$cell_age)], repeats.allowed = FALSE)}))

contrasts_celltype_age

nrow(contrasts_celltype_age)

###########
## Monte carlo the DE
###########


d10x_adult<-subset(d10x, subset = AgeGroup == "Adult")
ncol(d10x_adult)
d10x_ped<-subset(d10x, subset = AgeGroup == "Ped")
ncol(d10x_ped)


### paralize
command_args <- commandArgs(trailingOnly = TRUE)
cell_type_indx <- as.numeric(command_args[1])
cell_type<-cell_types[cell_type_indx]

samp_num=10


#DE_monte_carlo<-lapply(cell_types, function(cell_type){
  
  contrasts_celltype<-contrasts_celltype_age[grep(cell_type, contrasts_celltype_age)]
  
  d10x_adult_celltype<-subset(d10x_adult, subset = CellType_rough == cell_type)
  ncol(d10x_adult_celltype)
  d10x_ped_celltype<-subset(d10x_ped, subset = CellType_rough == cell_type)
  ncol(d10x_ped_celltype)

  
  de_lists<-sapply(1:samp_num, function(x){
    set.seed(x)
    
    ## make downsampled
    if(ncol(d10x_adult_celltype)<ncol(d10x_ped_celltype)){
      ped_cells_random <- d10x_ped_celltype[, sample(colnames(d10x_ped_celltype), size = ncol(d10x_adult_celltype), replace=F)]
      d10_DE<-merge(d10x_adult_celltype, ped_cells_random)
    }else{
      adult_cells_random <- d10x_adult_celltype[, sample(colnames(d10x_adult_celltype), size = ncol(d10x_ped_celltype), replace=F)]
      d10_DE<-merge(d10x_ped_celltype, adult_cells_random)}
    
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
DE_monte_carlo<-DE_monte_carlo[which(!(is.na(DE_monte_carlo$gene))),]

save(DE_monte_carlo, file=here("data",paste(cell_type,"adult_ped_diff_motecarlo_10.RData",sep="_")))


