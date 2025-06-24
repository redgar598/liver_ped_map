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


##############
## Differential expression with age in KC
##############
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_d10x_adult_ped_raw.rds"))

load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
d10x$Sex[which(d10x$individual%in%c("C113","C115"))]<-"M"

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")


d10x_raw_KC<-subset(d10x, subset = CellType_refined %in% c("KC Like"))
rm(d10x)
gc()

ggplot(d10x_raw_KC@meta.data, aes(Age))+geom_histogram(bins = 50)+facet_wrap(~Treatment)


Idents(d10x_raw_KC)<-d10x_raw_KC$age_condition
table(d10x_raw_KC$age_condition)

## age differential
de_ped_adult<-FindMarkers(d10x_raw_KC, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de_ped_adult<-de_ped_adult[which(de_ped_adult$p_val_adj < 0.005 & abs(de_ped_adult$avg_log2FC) > 1),]


## early late peds
d10x_raw_KC@meta.data$age_cats<-as.factor(paste(d10x_raw_KC@meta.data$Age,d10x_raw_KC@meta.data$Treatment ))
levels(d10x_raw_KC@meta.data$age_cats)<-c("Early IFALD", "Early IFALD", 
                                          "Early Healthy", "Early Healthy",
                                          "Late Healthy", "Late Healthy",
                                          "Early Healthy" , "Adult Healthy" ,
                                          "Adult Healthy", "Adult Healthy" ,
                                          "Adult Healthy", "Adult Healthy" ,
                                          "Adult Healthy" ,"Early Healthy" , 
                                          "Early IFALD")

Idents(d10x_raw_KC)<-d10x_raw_KC$age_cats
table(d10x_raw_KC$age_cats)

## age differential early late peds
de_late_early<-FindMarkers(d10x_raw_KC, ident.1 = "Late Healthy", ident.2 = "Early Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de_late_early<-de_late_early[which(de_late_early$p_val_adj < 0.005 & abs(de_late_early$avg_log2FC) > 1),]


## age differential early peds to adult
de_adult_early<-FindMarkers(d10x_raw_KC, ident.1 = "Adult Healthy", ident.2 = "Early Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de_adult_early<-de_adult_early[which(de_adult_early$p_val_adj < 0.005 & abs(de_adult_early$avg_log2FC) > 1),]


## overlap
intersect(rownames(sig_de_late_early), rownames(sig_de_ped_adult))
VlnPlot(d10x_raw_KC, group.by = "age_cats", features =intersect(rownames(sig_de_late_early), rownames(sig_de_ped_adult)))


intersect(rownames(sig_de_adult_early), rownames(sig_de_ped_adult))
VlnPlot(d10x_raw_KC, group.by = "age_cats", features =intersect(rownames(sig_de_adult_early), rownames(sig_de_ped_adult))[1:10])



####################
## specific analysis IFLAD versus peds in HSCs
####################
d10x<-readRDS(file = here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_d10x_adult_ped_raw.rds"))
load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
d10x$Sex[which(d10x$individual%in%c("C113","C115"))]<-"M"

d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

d10x_raw_hsc<-subset(d10x, subset = CellType_rough %in% c("HSC"))
rm(d10x)
gc()

d10x_raw_hsc@meta.data$age_cats_IFALD<-as.factor(paste(d10x_raw_hsc@meta.data$Age,d10x_raw_hsc@meta.data$Treatment ))
levels(d10x_raw_hsc@meta.data$age_cats_IFALD)<-c("0-11 months IFALD", "0-11 months IFALD", 
                                          "5-14 years Healthy", "5-14 years Healthy",
                                          "15-19 years Healthy", "15-19 years Healthy",
                                          "1-4 years Healthy" , "26 years Healthy" ,
                                          "57 years Healthy", "61 years Healthy" ,
                                          "65 years Healthy", "67 years Healthy" ,
                                          "69 years Healthy" ,"5-14 Healthy" , 
                                          "5-14 IFALD")

Idents(d10x_raw_hsc)<-d10x_raw_hsc$age_cats_IFALD
table(d10x_raw_hsc$age_cats_IFALD, d10x_raw_hsc$Age)



## age differential early late peds
de_early<-FindMarkers(d10x_raw_hsc, ident.1 = "0-11 months IFALD", ident.2 = "1-4 years Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de_early<-de_early[which(de_early$p_val_adj < 0.005 & abs(de_early$avg_log2FC) > 1),]

write.csv(sig_de_early, file=here("data","differential_IFALDHSC_only_youngest.csv"))


# ## age differential early peds to adult
# de_mid<-FindMarkers(d10x_raw_hsc, ident.1 = "Mid IFALD", ident.2 = "Mid Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
# sig_mid<-de_mid[which(de_mid$p_val_adj < 0.005 & abs(de_mid$avg_log2FC) > 1),]
# intersect(full_HSC_ped_IFALD$X, rownames(sig_mid))
# 
# sig_mid[grep("COL1A1|PDGFRA|CXCL12|IGFBP3", rownames(sig_mid)),]
# VlnPlot(d10x_raw_hsc, features=c("COL1A1","PDGFRA","CXCL12","IGFBP3"))
# 


## both comparisons confounded with sex differences


full_HSC_ped_IFALD<-read.csv(here("data/differential_IFALDHSC.csv"))

intersect(full_HSC_ped_IFALD$X, rownames(sig_de_early))

sig_de_early[grep("COL1A1|PDGFRA|CXCL12|IGFBP3", rownames(sig_de_early)),]
VlnPlot(d10x_raw_hsc, features=c("COL1A1","PDGFRA","CXCL12","IGFBP3"))


DefaultAssay(d10x_raw_hsc) <- "RNA"
meta_myeloid<-d10x_raw_hsc@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)

gene_exp<-FetchData(d10x_raw_hsc, vars=c("COL1A1","PDGFRA","CXCL12","IGFBP3"))
gene_exp$cell<-rownames(gene_exp)
plt_myeloid<-merge(meta_myeloid, gene_exp, by='cell')

plt_myeloid<-plt_myeloid[which(plt_myeloid$age_cats_IFALD%in%c("0-11 months IFALD", "1-4 years Healthy")),]

plt_myeloid<-melt(plt_myeloid[,c("age_cats_IFALD"  ,  "COL1A1", "PDGFRA"  ,  "CXCL12"  ,  "IGFBP3")])

early_only<-ggplot(plt_myeloid, aes(age_cats_IFALD,value))+
  geom_violin(fill="grey90",color="white")+geom_boxplot(aes(fill=age_cats_IFALD),width=0.1)+#fillscale_agecondition+
  theme_bw()+th_present+xlab("Age Group")+ylab("Expression")+facet_wrap(~variable)+
  theme(legend.position = "none")+scale_fill_manual(values=c("#B4EB65","cornflowerblue"))+ylim(0,7)
save_plts(early_only, "early_HSC_IFALD", w=5, h=5)
