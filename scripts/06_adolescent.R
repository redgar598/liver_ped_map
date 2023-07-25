####################
# Adolescent
####################


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


source("scripts/00_pretty_plots.R")
source("scripts/00_plot_gene_exp.R")
source("scripts/00_fanciest_UMAP.R")
source("scripts/00_entropy_d10x.R")


load(here("data","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))


d10x.combined_myeloid<-subset(d10x.combined, subset = CellType_refined %in% c("RR Myeloid","Macrophage\n(MHCII high)","KC Like","Macrophage\n(CLEC9A high)","Cycling Myeloid","Myeloid Erythrocytes\n(phagocytosis)"))
rm(d10x.combined)
gc()
d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)
d10x.combined_myeloid <- RunUMAP(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindNeighbors(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindClusters(d10x.combined_myeloid, resolution = 0.3)


d10x.combined_myeloid$age_group_adolescent<-d10x.combined_myeloid$Age
d10x.combined_myeloid$age_group_adolescent<-as.factor(d10x.combined_myeloid$age_group_adolescent)

relevel_age<-as.numeric(as.character(levels(d10x.combined_myeloid$age_group_adolescent)))
relevel_age<-sapply(1:length(relevel_age), function(x) if(relevel_age[x]>=16 & relevel_age[x]<18){"16-17"}else{if(relevel_age[x]>=18){">17"}else{"<16"}})
levels(d10x.combined_myeloid$age_group_adolescent)<-relevel_age
  
d10x.combined_myeloid$age_condition_adolescent<-paste(d10x.combined_myeloid$age_group_adolescent, d10x.combined_myeloid$Treatment)
d10x.combined_myeloid$age_condition_adolescent<-factor(d10x.combined_myeloid$age_condition_adolescent, levels=c("<16 IFALD","<16 Healthy","16-17 Healthy",">17 Healthy"))


table(d10x.combined_myeloid$age_condition_adolescent, d10x.combined_myeloid$CellType_refined)
table(d10x.combined_myeloid$age_condition_adolescent, d10x.combined_myeloid$individual)
table(d10x.combined_myeloid$Age, d10x.combined_myeloid$individual)


cell_num_myeloid<-as.data.frame(table(d10x.combined_myeloid$age_condition_adolescent))
colnames(cell_num_myeloid)<-c("age_condition_adolescent","CellCount")
myeloid_cluster_umap<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=F,split.by = "age_condition_adolescent", group.by = "CellType_refined", ncol=2)+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  geom_text(aes(y=-7, x=-9,label=paste0("n = ",comma(CellCount))),cell_num_myeloid, hjust=-0.1, size=3)
myeloid_cluster_umap
save_plts(myeloid_cluster_umap, "IFALD_myeloid_map_adolescent", w=7,h=6)



plt_entropy_age_condition<-entropy_d10(d10x.combined_myeloid, "age_condition_adolescent")
entropy_myeloidage<-entropy_plt(plt_entropy_age_condition, "age_condition_adolescent", d10x.combined_myeloid)
entropy_myeloidage
save_plts(entropy_myeloidage, "entropy_age_myeloid_adolescent", w=15,h=10)



##############
## Differential expression with age in KC
##############
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))

d10x$age_group_adolescent<-d10x$Age
d10x$age_group_adolescent<-as.factor(d10x$age_group_adolescent)
relevel_age<-as.numeric(as.character(levels(d10x$age_group_adolescent)))
relevel_age<-sapply(1:length(relevel_age), function(x) if(relevel_age[x]>=16 & relevel_age[x]<18){"16-17"}else{if(relevel_age[x]>=18){">17"}else{"<16"}})
levels(d10x$age_group_adolescent)<-relevel_age
d10x$age_condition_adolescent<-paste(d10x$age_group_adolescent, d10x$Treatment)
d10x$age_condition_adolescent<-factor(d10x$age_condition_adolescent, levels=c("<16 IFALD","<16 Healthy","16-17 Healthy",">17 Healthy"))

load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
d10x_raw_KC<-subset(d10x, subset = CellType_refined %in% c("KC Like"))

Idents(d10x_raw_KC)<-d10x_raw_KC$age_condition_adolescent
table(d10x_raw_KC$age_condition_adolescent)

## age differential
de<-FindMarkers(d10x_raw_KC, ident.1 = ">17 Healthy", ident.2 = "<16 Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]
sig_de_pos<-sig_de[which(sig_de$avg_log2FC>0),]
sig_de_neg<-sig_de[which(sig_de$avg_log2FC<0),]

KC_sig_de_full<-read.csv(here("data/differential_age_KC.csv"), row.names = 1)
KC_sig_de_full_pos<-KC_sig_de_full[which(KC_sig_de_full$avg_log2FC>0),]
KC_sig_de_full_neg<-KC_sig_de_full[which(KC_sig_de_full$avg_log2FC<0),]

# overall same direction 23/52 44%
length(intersect(rownames(KC_sig_de_full_pos),rownames(sig_de_pos))) #5/6
length(intersect(rownames(KC_sig_de_full_neg),rownames(sig_de_neg))) #18/46


### Plot key genes
gene_exp<-FetchData(d10x_raw_KC, vars=c("CCL4","CCL3","IL1B"))
gene_exp$cell<-rownames(gene_exp)  
gene_exp<-melt(gene_exp)
meta_myeloid<-d10x_raw_KC@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, gene_exp, by='cell')

ggplot(plt_myeloid, aes(age_condition_adolescent,log(value)))+
  geom_violin(fill="grey80",color="white")+geom_boxplot(aes(fill=age_condition_adolescent),width=0.1)+
  theme_bw()+th_present+xlab("Age Group")+ylab("Expression (log)")+
  theme(legend.position = "none")+facet_wrap(~variable)

ggplot(plt_myeloid, aes(Age,log(value)))+
  geom_point()+
  theme_bw()+th_present+xlab("Age Group")+ylab("Expression (log)")+
  theme(legend.position = "none")+facet_wrap(Treatment~variable)+stat_smooth(method="lm")




##############
## Differential expression with age and IFALD in MHCII
##############
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))

d10x$age_group_adolescent<-d10x$Age
d10x$age_group_adolescent<-as.factor(d10x$age_group_adolescent)
relevel_age<-as.numeric(as.character(levels(d10x$age_group_adolescent)))
relevel_age<-sapply(1:length(relevel_age), function(x) if(relevel_age[x]>=16 & relevel_age[x]<18){"16-17"}else{if(relevel_age[x]>=18){">17"}else{"<16"}})
levels(d10x$age_group_adolescent)<-relevel_age
d10x$age_condition_adolescent<-paste(d10x$age_group_adolescent, d10x$Treatment)
d10x$age_condition_adolescent<-factor(d10x$age_condition_adolescent, levels=c("<16 IFALD","<16 Healthy","16-17 Healthy",">17 Healthy"))

load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
d10x_raw_mhcII<-subset(d10x, subset = CellType_refined %in% c("Macrophage\n(MHCII high)"))

Idents(d10x_raw_mhcII)<-d10x_raw_mhcII$age_condition_adolescent
table(d10x_raw_mhcII$age_condition_adolescent)

## age differential
de<-FindMarkers(d10x_raw_mhcII, ident.1 = ">17 Healthy", ident.2 = "<16 Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]
sig_de_pos<-sig_de[which(sig_de$avg_log2FC>0),]
sig_de_neg<-sig_de[which(sig_de$avg_log2FC<0),]

MHCII_sig_de_full<-read.csv(here("data/differential_age_MHCII.csv"), row.names = 1)
MHCII_sig_de_full_pos<-MHCII_sig_de_full[which(MHCII_sig_de_full$avg_log2FC>0),]
MHCII_sig_de_full_neg<-MHCII_sig_de_full[which(MHCII_sig_de_full$avg_log2FC<0),]

# overall same direction 41/50 82%
length(intersect(rownames(MHCII_sig_de_full_pos),rownames(sig_de_pos))) #5/5
length(intersect(rownames(MHCII_sig_de_full_neg),rownames(sig_de_neg))) #36/45



##############
## Differential expression with age and IFALD in RR
##############
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))

d10x$age_group_adolescent<-d10x$Age
d10x$age_group_adolescent<-as.factor(d10x$age_group_adolescent)
relevel_age<-as.numeric(as.character(levels(d10x$age_group_adolescent)))
relevel_age<-sapply(1:length(relevel_age), function(x) if(relevel_age[x]>=16 & relevel_age[x]<18){"16-17"}else{if(relevel_age[x]>=18){">17"}else{"<16"}})
levels(d10x$age_group_adolescent)<-relevel_age
d10x$age_condition_adolescent<-paste(d10x$age_group_adolescent, d10x$Treatment)
d10x$age_condition_adolescent<-factor(d10x$age_condition_adolescent, levels=c("<16 IFALD","<16 Healthy","16-17 Healthy",">17 Healthy"))

load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
d10x_raw_RR<-subset(d10x, subset = CellType_refined %in% c("RR Myeloid"))

Idents(d10x_raw_RR)<-d10x_raw_RR$age_condition_adolescent
table(d10x_raw_RR$age_condition_adolescent)

## age differential
de<-FindMarkers(d10x_raw_RR, ident.1 = ">17 Healthy", ident.2 = "<16 Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]
sig_de_pos<-sig_de[which(sig_de$avg_log2FC>0),]
sig_de_neg<-sig_de[which(sig_de$avg_log2FC<0),]

RR_sig_de_full<-read.csv(here("data/differential_age_RR.csv"), row.names = 1)
RR_sig_de_full_pos<-RR_sig_de_full[which(RR_sig_de_full$avg_log2FC>0),]
RR_sig_de_full_neg<-RR_sig_de_full[which(RR_sig_de_full$avg_log2FC<0),]

# overall same direction 28/45 62%
length(intersect(rownames(RR_sig_de_full_pos),rownames(sig_de_pos))) #6/7
length(intersect(rownames(RR_sig_de_full_neg),rownames(sig_de_neg))) #22/38

