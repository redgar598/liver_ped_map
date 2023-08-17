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
library(RColorBrewer)

source("scripts/00_pretty_plots.R")
source("scripts/00_entropy_d10x.R")



################
## Differential Expression
################
load(here("/media/redgar/Seagate Portable Drive/processed_data/Fetal_IFALD_adult_ped_raw_myeloid_only.RData"))

levels(d10x_myeloid$CellType_harmonized)[which(levels(d10x_myeloid$CellType_harmonized)=="Macrophage\n(MHCII high)")]<-"Marcophage_MHCII"

d10x_myeloid$CellType_harmonized<-as.character(d10x_myeloid$CellType_harmonized)
d10x_myeloid$cell_age<-paste(d10x_myeloid$CellType_harmonized, d10x_myeloid$age_condition, sep = "_")
Idents(d10x_myeloid) <- "cell_age"

table(d10x_myeloid$CellType_harmonized, d10x_myeloid$age_condition)

#MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene,
#consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
#and linear regression for the continuous process (i.e., the expression level). 

cell_types<-unique(as.character(d10x_myeloid$CellType_harmonized))

# too few Macrophage (CLEC9A high) 
cell_types<-cell_types[-grep("CLEC9A ",cell_types)]


contrasts_celltype_age<-do.call(rbind,lapply(1:length(cell_types), function(x){
  combinations(n = 3, r = 2, v = d10x_myeloid$cell_age[grep(cell_types[x],d10x_myeloid$cell_age)], repeats.allowed = FALSE)}))

contrasts_celltype_age

nrow(contrasts_celltype_age)


###########
## Monte carlo the DE
###########
d10x_adult<-subset(d10x_myeloid, subset = age_condition == "Adult Healthy")
ncol(d10x_adult)
d10x_ped<-subset(d10x_myeloid, subset = age_condition == "Ped Healthy")
ncol(d10x_ped)
d10x_fetal<-subset(d10x_myeloid, subset = age_condition == "Fetal Healthy")
ncol(d10x_fetal)


### paralize
command_args <- commandArgs(trailingOnly = TRUE)
cell_type_indx <- as.numeric(command_args[1])
cell_type<-cell_types[cell_type_indx]

samp_num=1000


#DE_monte_carlo<-lapply(cell_types, function(cell_type){

contrasts_celltype<-as.data.frame(contrasts_celltype_age)[grep(cell_type, as.data.frame(contrasts_celltype_age)$V1),]

d10x_adult_celltype<-subset(d10x_adult, subset = CellType_harmonized == cell_type)
ncol(d10x_adult_celltype)
d10x_ped_celltype<-subset(d10x_ped, subset = CellType_harmonized == cell_type)
ncol(d10x_ped_celltype)
d10x_fetal_celltype<-subset(d10x_fetal, subset = CellType_harmonized == cell_type)
ncol(d10x_fetal_celltype)


de_lists<-sapply(1:samp_num, function(x){
  set.seed(x)
  
  smallest<-c("adult","ped","fetal")[which(c(ncol(d10x_adult_celltype),ncol(d10x_ped_celltype),ncol(d10x_fetal_celltype))==min(c(ncol(d10x_adult_celltype),ncol(d10x_ped_celltype),ncol(d10x_fetal_celltype))))]

  ## make downsampled
  if(smallest=="adult"){
    ped_cells_random <- d10x_ped_celltype[, sample(colnames(d10x_ped_celltype), size = ncol(d10x_adult_celltype), replace=F)]
    fetal_cells_random <- d10x_fetal_celltype[, sample(colnames(d10x_fetal_celltype), size = ncol(d10x_adult_celltype), replace=F)]
    d10_DE<-merge(fetal_cells_random, ped_cells_random)
    d10_DE<-merge(d10x_adult_celltype, d10_DE)
  }else{
    if(smallest=="ped"){
      fetal_cells_random <- d10x_fetal_celltype[, sample(colnames(d10x_fetal_celltype), size = ncol(d10x_ped_celltype), replace=F)]
      adult_cells_random <- d10x_adult_celltype[, sample(colnames(d10x_adult_celltype), size = ncol(d10x_ped_celltype), replace=F)]
      d10_DE<-merge(fetal_cells_random, adult_cells_random)
      d10_DE<-merge(d10x_ped_celltype, d10_DE)
    }else{
      if(smallest=="fetal"){
        ped_cells_random <- d10x_ped_celltype[, sample(colnames(d10x_ped_celltype), size = ncol(d10x_fetal_celltype), replace=F)]
        adult_cells_random <- d10x_adult_celltype[, sample(colnames(d10x_adult_celltype), size = ncol(d10x_fetal_celltype), replace=F)]
        d10_DE<-merge(adult_cells_random, ped_cells_random)
        d10_DE<-merge(d10x_fetal_celltype, d10_DE)
      }}}
  
  ## run DE 
  de<-FindMarkers(d10_DE, ident.1 = contrasts_celltype$V1[3], ident.2 = contrasts_celltype$V2[3], test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
  print(paste(contrasts_celltype[1],"vs", contrasts_celltype[2],":", nrow(de), sep=" "))
  de$gene<-rownames(de)
  rownames(de)<-NULL
  de<-de[,c(6,1:5)]
  de$cell.1<-contrasts_celltype$V1[3]
  de$cell.2<-contrasts_celltype$V2[3]
  
  # more relaxed FDR and no FC 
  #de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]$gene
  de[which(de$p_val_adj < 0.05),]$gene
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

save(DE_monte_carlo, file=here("data",paste(cell_type,"adult_ped_diff_motecarlo_1000_covarSex.RData",sep="_")))




## age differential KC
d10x_raw_KC<-subset(d10x_myeloid, subset = CellType_harmonized %in% c("KC Like"))

DefaultAssay(d10x_raw_KC) <- "RNA"
Idents(d10x_raw_KC)<-d10x_raw_KC$age_condition
table(d10x_raw_KC$age_condition)

de_fetal_ped_KC<-FindMarkers(d10x_raw_KC, ident.1 = "Fetal Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
de_ped_adult_KC<-FindMarkers(d10x_raw_KC, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
de_fetal_adult_KC<-FindMarkers(d10x_raw_KC, ident.1 = "Fetal Healthy", ident.2 = "Adult Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)


## age differential RR
d10x_raw_RR<-subset(d10x_myeloid, subset = CellType_harmonized %in% c("RR Myeloid"))

DefaultAssay(d10x_raw_RR) <- "RNA"
Idents(d10x_raw_RR)<-d10x_raw_RR$age_condition
table(d10x_raw_RR$age_condition)

de_fetal_ped_RR<-FindMarkers(d10x_raw_RR, ident.1 = "Fetal Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
de_ped_adult_RR<-FindMarkers(d10x_raw_RR, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
de_fetal_adult_RR<-FindMarkers(d10x_raw_RR, ident.1 = "Fetal Healthy", ident.2 = "Adult Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)

## age differential MHCII
d10x_raw_MHC<-subset(d10x_myeloid, subset = CellType_harmonized %in% c("Macrophage\n(MHCII high)"))

DefaultAssay(d10x_raw_MHC) <- "RNA"
Idents(d10x_raw_MHC)<-d10x_raw_MHC$age_condition
table(d10x_raw_MHC$age_condition)

de_fetal_ped_MHC<-FindMarkers(d10x_raw_MHC, ident.1 = "Fetal Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
de_ped_adult_MHC<-FindMarkers(d10x_raw_MHC, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
de_fetal_adult_MHC<-FindMarkers(d10x_raw_MHC, ident.1 = "Fetal Healthy", ident.2 = "Adult Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
