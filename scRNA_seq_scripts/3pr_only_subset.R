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



Idents(d10x_raw_KC)<-d10x_raw_KC$age_condition
table(d10x_raw_KC$age_condition)

## age differential
de<-FindMarkers(d10x_raw_KC, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]




####################
## just 3' n=3 peds and n=6 adults
####################
d10x_raw_KC_3pr<-subset(d10x_raw_KC, subset = chemistry %in% c("3pr"))

## age differential
de_3pr<-FindMarkers(d10x_raw_KC_3pr, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de_3pr<-de_3pr[which(de_3pr$p_val_adj < 0.005 & abs(de_3pr$avg_log2FC) > 1),]



### overlap
length(rownames(sig_de_3pr))
length(rownames(sig_de))
length(intersect(rownames(sig_de), rownames(sig_de_3pr)))

49/72

de_3pr$gene<-rownames(de_3pr)
colnames(de_3pr)<-paste(colnames(de_3pr),"_3pr", sep = "")
de$gene<-rownames(de)


de_combine_plt<-merge(de_3pr,de, by.x="gene_3pr", by.y="gene")
cor_r<-cor.test(de_combine_plt$avg_log2FC, de_combine_plt$avg_log2FC_3pr, method = "spearman")
cor_r$estimate

de_combine_plt$Sig<-"Neither"
de_combine_plt$Sig[which(de_combine_plt$gene_3pr%in%rownames(sig_de))]<-"Only in Full (5' and 3') Cohort"
de_combine_plt$Sig[which(de_combine_plt$gene_3pr%in%rownames(sig_de_3pr))]<-"Only in 3' Samples"
de_combine_plt$Sig[which(de_combine_plt$gene_3pr%in%intersect(rownames(sig_de), rownames(sig_de_3pr)))]<-"Both"


de_combine_plt_cor<-ggplot(de_combine_plt, aes(avg_log2FC,avg_log2FC_3pr))+
  geom_vline(xintercept = c(1,-1), color="grey")+stat_smooth(method="lm", se=F, color="black")+
  geom_point(shape=21, aes(fill=Sig))+theme_bw()+
  scale_fill_manual(values=c("#f09c16","grey","#7fc97f","#386cb0"),name="Significant in which cohort")+
  theme(legend.position = "bottom")+guides(fill=guide_legend(ncol=2))+
  xlab("Fold Change (log2)\nFull (5' and 3') Cohort")+  ylab("Fold Change (log2)\n3' Samples Only")+
  geom_hline(yintercept = c(1,-1), color="grey")+xlim(-4,4)+ylim(-4,4)+
  annotate("text", x=3, y=-3, label= paste("italic(r[s]) ==  ", round(cor_r$estimate,2)),
           parse = TRUE)
de_combine_plt_cor
save_plts(de_combine_plt_cor, "FC_cor_3pr", w=6, h=6.5)  


