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


source("scripts/00_pretty_plots.R")
source("scripts/00_fanciest_UMAP.R")
source("scripts/00_plot_gene_exp.R")
source("scripts/00_entropy_d10x.R")



load(here("data","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

fancyUMAP_all<-fanciest_UMAP(d10x.combined,"KC Like",F)
save_plts(fancyUMAP_all, "IFALD_KC_highlight_umap_fancy", w=6,h=4)

d10x.combined_myeloid<-subset(d10x.combined, subset = CellType_refined %in% c("RR Myeloid","Macrophage\n(MHCII high)","KC Like","Macrophage\n(CLEC9A high)","Cycling Myeloid","Myeloid Erythrocytes\n(phagocytosis)"))
rm(d10x.combined)
gc()
d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)
d10x.combined_myeloid <- RunUMAP(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindNeighbors(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindClusters(d10x.combined_myeloid, resolution = 0.3)

myeloid_subtype<-DimPlot(d10x.combined_myeloid, label=T)
myeloid_subtype
save_plts(myeloid_subtype, "IFALD_myeloid_map_clusters", w=7,h=6)

myeloid_cluster_umap<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=T, group.by = "CellType_refined")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-7, y=-9, label=paste0("n = ",comma(ncol(d10x.combined_myeloid))))
myeloid_cluster_umap
save_plts(myeloid_cluster_umap, "IFALD_myeloid_map_celltype", w=7,h=6)

cell_num_myeloid<-as.data.frame(table(d10x.combined_myeloid$age_condition))
colnames(cell_num_myeloid)<-c("age_condition","CellCount")
myeloid_cluster_umap<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=F,split.by = "age_condition", group.by = "CellType_refined", ncol=2)+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  geom_text(aes(y=-7, x=-9,label=paste0("n = ",comma(CellCount))),cell_num_myeloid, hjust=-0.1, size=3)
myeloid_cluster_umap
save_plts(myeloid_cluster_umap, "IFALD_myeloid_map", w=7,h=6)

myeloid_cluster_umap<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=F,split.by = "age_id", group.by = "CellType_refined", ncol=4)+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-14, y=-12, label=paste0("n = ",comma(ncol(d10x.combined_myeloid))))
myeloid_cluster_umap
save_plts(myeloid_cluster_umap, "IFALD_myeloid_map_individual", w=12,h=10)

plt_entropy_age_condition<-entropy_d10(d10x.combined_myeloid, "age_condition")
entropy_myeloidage<-entropy_plt(plt_entropy_age_condition, "age_condition", d10x.combined_myeloid)
entropy_myeloidage
save_plts(entropy_myeloidage, "entropy_age_myeloid", w=15,h=10)

#################
### add PC
#################
embed_PC12<-as.data.frame(Embeddings(d10x.combined_myeloid, reduction = "pca"))[,1:2]
identical(rownames(embed_PC12), colnames(d10x.combined_myeloid))
d10x.combined_myeloid <- AddMetaData(d10x.combined_myeloid, metadata = embed_PC12)

FeaturePlot(d10x.combined_myeloid, reduction="umap", features = c("PC_2","ALB"))
FeaturePlot(d10x.combined_myeloid, reduction="umap", features = c("PC_1","ALB"))

DimPlot(d10x.combined_myeloid, reduction = "pca", pt.size=0.25, label=F, group.by = "CellType_refined", split.by = "age_condition")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-7, y=-9, label=paste0("n = ",comma(ncol(d10x.combined_myeloid))))

FeaturePlot(d10x.combined_myeloid, reduction="pca", features = c("ALB"))

ggplot(d10x.combined_myeloid@meta.data, aes(age_condition, PC_2))+
  geom_violin(fill="lightgrey",color="white")+geom_boxplot(aes(fill=age_condition),width=0.25)+facet_wrap(~CellType_refined)+
  fillscale_agecondition+theme_bw()


Loadings<-as.data.frame(Loadings(d10x.combined_myeloid, reduction = "pca"))
embed<-as.data.frame(Embeddings(d10x.combined_myeloid, reduction = "pca"))
vars <- Stdev(d10x.combined_myeloid, reduction = "pca")^2
Importance<-vars/sum(vars)
print(Importance[1:10])

rownames(Loadings)[order(Loadings$PC_2)][1:20]
rownames(Loadings)[rev(order(Loadings$PC_2))][1:20]

FeaturePlot(d10x.combined_myeloid, reduction="pca", features = c("CXCL12"))


ggplot(d10x.combined_myeloid@meta.data, aes(x=PC_1, y=PC_2)) +
  geom_point(aes(color=CellType_refined),alpha=0.1, shape=16, size=2) + 
  geom_density_2d(bins=10, color="grey30")+
  theme_bw()+colscale_cellType+facet_wrap(~age_condition)

FeaturePlot(d10x.combined_myeloid, reduction="pca", features = c("HLA-DRA"), split.by = "age_condition")


plot_gene_PCA(d10x.combined_myeloid, "HLA-DRA",0.7, T, "Macrophage\n(MHCII high)")
plot_gene_PCA(d10x.combined_myeloid, "MARCO",0.7, T, "Macrophage\n(MHCII high)")
plot_gene_PCA(d10x.combined_myeloid, "VSIG4",0.6, T, "Macrophage\n(MHCII high)")
plot_gene_PCA(d10x.combined_myeloid, "S100A8",0.6, T, "Macrophage\n(MHCII high)")
plot_gene_PCA(d10x.combined_myeloid, "S100A9",0.6, T, "Macrophage\n(MHCII high)")
plot_gene_PCA(d10x.combined_myeloid, "CD163",0.6, T, "Macrophage\n(MHCII high)")
plot_gene_PCA(d10x.combined_myeloid, "MRC1",0.6, T, "Macrophage\n(MHCII high)")
plot_gene_PCA(d10x.combined_myeloid, "LYVE1",0.6, T, "Macrophage\n(MHCII high)")


##############
## Differential expression with age in KC
##############
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))

load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")


d10x_raw_KC<-subset(d10x, subset = CellType_refined %in% c("KC Like"))

Idents(d10x_raw_KC)<-d10x_raw_KC$age_condition
table(d10x_raw_KC$age_condition)

## age differential
de<-FindMarkers(d10x_raw_KC, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]
sig_de[which(sig_de$avg_log2FC>0),]
sig_de[which(sig_de$avg_log2FC<0),]

### IFALD differential
de_IFALD<-FindMarkers(d10x_raw_KC, ident.1 = "Ped IFALD", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de_IFALD<-de_IFALD[which(de_IFALD$p_val_adj < 0.005 & abs(de_IFALD$avg_log2FC) > 1),]
sig_de_IFALD[which(sig_de_IFALD$avg_log2FC>0),]
sig_de_IFALD[which(sig_de_IFALD$avg_log2FC<0),]


###
## pathway adult/IFALD versus healthy ped
###
source("scripts/00_GSEA_function.R")
#http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/
GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

### Age
de$gene<-rownames(de)
gene_list = de$avg_log2FC
names(gene_list) = de$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

res = GSEA(gene_list, GO_file, pval = 0.05)

plt_path<-res$Results
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Enrichment_Cell<-"Up-regulated in \nAdult"
plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \nPed"

plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))

plt_path$direction_label<-as.factor(plt_path$Enrichment)
levels(plt_path$direction_label)<-c(0.1,-0.1)
plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))

# top and bottom 15
plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])

myeloid_GSEA<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_hline(yintercept=16.5, color="grey")+scale_fill_manual(values=c("#D64A56","cornflowerblue"))
myeloid_GSEA
save_plts(myeloid_GSEA, "GSEA_Age_KC_myeloid", w=15,h=7)

## IFALD
de_IFALD$gene<-rownames(de_IFALD)
gene_list = de_IFALD$avg_log2FC
names(gene_list) = de_IFALD$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

res = GSEA(gene_list, GO_file, pval = 0.05)

plt_path<-res$Results
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Enrichment_Cell<-"Up-regulated in \nIFALD"
plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \nHealthy"

plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))

plt_path$direction_label<-as.factor(plt_path$Enrichment)
levels(plt_path$direction_label)<-c(0.1,-0.1)
plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))

# top and bottom 15
plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])

myeloid_GSEA<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_hline(yintercept=16.5, color="grey")+scale_fill_manual(values=c("cornflowerblue","#374eb8"))
myeloid_GSEA
save_plts(myeloid_GSEA, "GSEA_IFALD_KC_myeloid", w=15,h=7)

###############
## This works with enrichment map as the generic file
##############
plt_path<-res$Results
plt_path$Description<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Phenotype<-"+1"
plt_path$Phenotype[which(plt_path$Enrichment=="Down-regulated")]<-"-1"
plt_path$leadingEdge<-NULL

plt_path$GO.ID<-plt_path$pathway
plt_path$FDR<-plt_path$padj
plt_path$p.Val<-plt_path$pval

plt_path<-plt_path[,c("GO.ID","Description","p.Val","FDR")]

write.table(plt_path, file=here("data/KC_age.txt"), sep="\t", row.names = F, quote = F)


##############
## Differential expression with age and IFALD in RR
##############
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))

load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

d10x_raw_RR<-subset(d10x, subset = CellType_refined %in% c("RR Myeloid"))

Idents(d10x_raw_RR)<-d10x_raw_RR$age_condition
table(d10x_raw_RR$age_condition)

## age differential
de<-FindMarkers(d10x_raw_RR, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]
sig_de[which(sig_de$avg_log2FC>0),]
sig_de[which(sig_de$avg_log2FC<0),]

### IFALD differential
de_IFALD<-FindMarkers(d10x_raw_RR, ident.1 = "Ped IFALD", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de_IFALD<-de_IFALD[which(de_IFALD$p_val_adj < 0.005 & abs(de_IFALD$avg_log2FC) > 1),]
sig_de_IFALD[which(sig_de_IFALD$avg_log2FC>0),]
sig_de_IFALD[which(sig_de_IFALD$avg_log2FC<0),]


###
## pathway adult/IFALD versus healthy ped
###
source("scripts/00_GSEA_function.R")
GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

### Age
de$gene<-rownames(de)
gene_list = de$avg_log2FC
names(gene_list) = de$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

res = GSEA(gene_list, GO_file, pval = 0.05)

plt_path<-res$Results
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Enrichment_Cell<-"Up-regulated in \nAdult"
plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \nPed"

plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))

plt_path$direction_label<-as.factor(plt_path$Enrichment)
levels(plt_path$direction_label)<-c(0.1,-0.1)
plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))

# top and bottom 15
plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])

myeloid_GSEA<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_hline(yintercept=16.5, color="grey")+scale_fill_manual(values=c("#D64A56","cornflowerblue"))
myeloid_GSEA
save_plts(myeloid_GSEA, "GSEA_Age_RR_myeloid", w=15,h=7)

## IFALD:: ran nothing sig
de_IFALD$gene<-rownames(de_IFALD)
gene_list = de_IFALD$avg_log2FC
names(gene_list) = de_IFALD$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

res = GSEA(gene_list, GO_file, pval = 0.05)
res

plt_path<-res$Results
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Enrichment_Cell<-"Up-regulated in \nIFALD"
plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \nHealthy"

plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))

plt_path$direction_label<-as.factor(plt_path$Enrichment)
levels(plt_path$direction_label)<-c(0.1,-0.1)
plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))

# top and bottom 15
plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])

myeloid_GSEA<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_hline(yintercept=16.5, color="grey")+scale_fill_manual(values=c("cornflowerblue","#374eb8"))
myeloid_GSEA
save_plts(myeloid_GSEA, "GSEA_IFALD_RR_myeloid", w=15,h=7)



##############
## Differential expression with age and IFALD in MHCII
##############
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))

load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

d10x_raw_mhcII<-subset(d10x, subset = CellType_refined %in% c("Macrophage\n(MHCII high)"))

Idents(d10x_raw_mhcII)<-d10x_raw_mhcII$age_condition
table(d10x_raw_mhcII$age_condition)

## age differential
de<-FindMarkers(d10x_raw_mhcII, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]
sig_de[which(sig_de$avg_log2FC>0),]
sig_de[which(sig_de$avg_log2FC<0),]

### IFALD differential
de_IFALD<-FindMarkers(d10x_raw_mhcII, ident.1 = "Ped IFALD", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de_IFALD<-de_IFALD[which(de_IFALD$p_val_adj < 0.005 & abs(de_IFALD$avg_log2FC) > 1),]
sig_de_IFALD[which(sig_de_IFALD$avg_log2FC>0),]
sig_de_IFALD[which(sig_de_IFALD$avg_log2FC<0),]



###
## pathway adult/IFALD versus healthy ped
###
source("scripts/00_GSEA_function.R")
GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

### Age
de$gene<-rownames(de)
gene_list = de$avg_log2FC
names(gene_list) = de$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

res = GSEA(gene_list, GO_file, pval = 0.05)

plt_path<-res$Results
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Enrichment_Cell<-"Up-regulated in \nAdult"
plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \nPed"

plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))

plt_path$direction_label<-as.factor(plt_path$Enrichment)
levels(plt_path$direction_label)<-c(0.1,-0.1)
plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))

# top and bottom 15
plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])

myeloid_GSEA<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_hline(yintercept=16.5, color="grey")+scale_fill_manual(values=c("#D64A56","cornflowerblue"))
myeloid_GSEA
save_plts(myeloid_GSEA, "GSEA_Age_MHCII_myeloid", w=15,h=7)


## IFALD
de_IFALD$gene<-rownames(de_IFALD)
gene_list = de_IFALD$avg_log2FC
names(gene_list) = de_IFALD$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

res = GSEA(gene_list, GO_file, pval = 0.05)

plt_path<-res$Results
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Enrichment_Cell<-"Up-regulated in \nIFALD"
plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \nHealthy"

plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))

plt_path$direction_label<-as.factor(plt_path$Enrichment)
levels(plt_path$direction_label)<-c(0.1,-0.1)
plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))

# top and bottom 15
#plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])

myeloid_GSEA<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_hline(yintercept=16.5, color="grey")+scale_fill_manual(values=c("cornflowerblue","#374eb8"))
myeloid_GSEA
save_plts(myeloid_GSEA, "GSEA_IFALD_MHCII_myeloid", w=15,h=7)




##############
## Pathway heatmap with age and IFALD in all cell types
##############
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))

load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")


d10x_raw_KC<-subset(d10x, subset = CellType_refined %in% c("KC Like"))
d10x_raw_RR<-subset(d10x, subset = CellType_refined %in% c("RR Myeloid"))
d10x_raw_mhcII<-subset(d10x, subset = CellType_refined %in% c("Macrophage\n(MHCII high)"))

Idents(d10x_raw_KC)<-d10x_raw_KC$age_condition
Idents(d10x_raw_RR)<-d10x_raw_RR$age_condition
Idents(d10x_raw_mhcII)<-d10x_raw_mhcII$age_condition

## age differential
de_RR<-FindMarkers(d10x_raw_RR, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
de_KC<-FindMarkers(d10x_raw_KC, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
de_MHCII<-FindMarkers(d10x_raw_mhcII, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)

sig_de_age_RR<-de_RR[which(de_RR$p_val_adj < 0.005 & abs(de_RR$avg_log2FC) > 1),]
sig_de_age_KC<-de_KC[which(de_KC$p_val_adj < 0.005 & abs(de_KC$avg_log2FC) > 1),]
sig_de_age_MHCII<-de_MHCII[which(de_MHCII$p_val_adj < 0.005 & abs(de_MHCII$avg_log2FC) > 1),]

write.csv(sig_de_age_RR, file=here("data","differential_age_RR.csv"))
write.csv(sig_de_age_KC, file=here("data","differential_age_KC.csv"))
write.csv(sig_de_age_MHCII, file=here("data","differential_age_MHCII.csv"))

write.csv(de_RR, file=here("data","all_gene_stats_age_RR.csv"))
write.csv(de_KC, file=here("data","all_gene_stats_age_KC.csv"))
write.csv(de_MHCII, file=here("data","all_gene_stats_age_MHCII.csv"))

## Overlap table
sig_de_age_KC$cell<-"KC Like"
sig_de_age_RR$cell<-"RR Myeloid"
sig_de_age_MHCII$cell<-"Macrophage (MHCII high)"

sig_de_age_KC$gene<-rownames(sig_de_age_KC)
sig_de_age_RR$gene<-rownames(sig_de_age_RR)
sig_de_age_MHCII$gene<-rownames(sig_de_age_MHCII)

sig_de_age<-rbind(sig_de_age_KC,sig_de_age_RR, sig_de_age_MHCII)

sig_de_age_pos<-sig_de_age[which(sig_de_age$avg_log2FC>0),]
sig_de_age_neg<-sig_de_age[which(sig_de_age$avg_log2FC<0),]

summary_tbl_pos<-as.data.frame(sig_de_age_pos %>%
                             select(gene, cell) %>% 
                             group_by(gene) %>%
                             mutate(all_cells = paste(cell, collapse = " | "))%>%
                             select(gene, all_cells))
summary_tbl_pos<-summary_tbl_pos[!duplicated(summary_tbl_pos),]
summary_tbl_pos$all_cells<-gsub("\n"," ", summary_tbl_pos$all_cells)
summary_tbl_pos$direction<-"Up in adults"

summary_tbl_neg<-as.data.frame(sig_de_age_neg %>%
                             select(gene, cell) %>% 
                             group_by(gene) %>%
                             mutate(all_cells = paste(cell, collapse = " | "))%>%
                             select(gene, all_cells))
summary_tbl_neg<-summary_tbl_neg[!duplicated(summary_tbl_neg),]
summary_tbl_neg$all_cells<-gsub("\n"," ", summary_tbl_neg$all_cells)
summary_tbl_neg$direction<-"Up in peds"

summary_tbl<-rbind(summary_tbl_neg, summary_tbl_pos)

summary_tbl$all_cells<-factor(summary_tbl$all_cells, levels=c("KC Like | RR Myeloid | Macrophage (MHCII high)",
                                                              "KC Like | RR Myeloid", "KC Like | Macrophage (MHCII high)", "RR Myeloid | Macrophage (MHCII high)",
                                                              "KC Like" , "RR Myeloid","Macrophage (MHCII high)"))
summary_tbl[order(summary_tbl$direction, summary_tbl$all_cells),]

write.csv(file=here("data","Significant_genes_adult_ped_myeloid.csv"),summary_tbl[order(summary_tbl$direction, summary_tbl$all_cells),])

## IFALD differential
de_IFALD_RR<-FindMarkers(d10x_raw_RR, ident.1 = "Ped IFALD", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
de_IFALD_KC<-FindMarkers(d10x_raw_KC, ident.1 = "Ped IFALD", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
de_IFALD_MHCII<-FindMarkers(d10x_raw_mhcII, ident.1 = "Ped IFALD", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)

sig_de_IFALD_RR<-de_IFALD_RR[which(de_IFALD_RR$p_val_adj < 0.005 & abs(de_IFALD_RR$avg_log2FC) > 1),]
sig_de_IFALD_KC<-de_IFALD_KC[which(de_IFALD_KC$p_val_adj < 0.005 & abs(de_IFALD_KC$avg_log2FC) > 1),]
sig_de_IFALD_MHCII<-de_IFALD_MHCII[which(de_IFALD_MHCII$p_val_adj < 0.005 & abs(de_IFALD_MHCII$avg_log2FC) > 1),]

write.csv(sig_de_IFALD_RR, file=here("data","differential_IFALD_RR.csv"))
write.csv(sig_de_IFALD_KC, file=here("data","differential_IFALD_KC.csv"))
write.csv(sig_de_IFALD_MHCII, file=here("data","differential_IFLAD_MHCII.csv"))


## pathway adult/IFALD versus healthy ped
source("scripts/00_GSEA_function.R")
GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

### Age
de_RR$gene<-rownames(de_RR)
gene_list = de_RR$avg_log2FC
names(gene_list) = de_RR$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
res_RR = GSEA(gene_list, GO_file, pval = 0.05)
res_RR<-res_RR$Results
res_RR$test<-"RR"

de_KC$gene<-rownames(de_KC)
gene_list = de_KC$avg_log2FC
names(gene_list) = de_KC$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
res_KC = GSEA(gene_list, GO_file, pval = 0.05)
res_KC<-res_KC$Results
res_KC$test<-"KC"

de_MHCII$gene<-rownames(de_MHCII)
gene_list = de_MHCII$avg_log2FC
names(gene_list) = de_MHCII$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
res_MHCII = GSEA(gene_list, GO_file, pval = 0.05)
res_MHCII<-res_MHCII$Results
res_MHCII$test<-"MHCII"

myeloid_res_overlap<-rbind(res_MHCII, res_KC, res_RR)

myeloid_res_overlap$pathway<-sapply(1:nrow(myeloid_res_overlap), function(x) strsplit(myeloid_res_overlap$pathway[x], "%")[[1]][1])
myeloid_res_overlap$Enrichment_Cell<-"Up-regulated in \nAdult"
myeloid_res_overlap$Enrichment_Cell[which(myeloid_res_overlap$Enrichment=="Down-regulated")]<-"Up-regulated in \nPed"

myeloid_res_overlap$label<-lapply(1:nrow(myeloid_res_overlap), function(x) paste0(myeloid_res_overlap$leadingEdge[x][[1]][1:4], collapse = ", "))

myeloid_res_overlap$direction_label<-as.factor(myeloid_res_overlap$Enrichment)
levels(myeloid_res_overlap$direction_label)<-c(0.1,-0.1)
myeloid_res_overlap$direction_label<-as.numeric(as.character(myeloid_res_overlap$direction_label))

# common<-names(table(myeloid_res_overlap$pathway)[which(table(myeloid_res_overlap$pathway)>2)])
# myeloid_res_overlap<-myeloid_res_overlap[which(myeloid_res_overlap$pathway%in%common),]
# 
# myeloid_GSEA_common<-ggplot(myeloid_res_overlap, aes(test, NES))+
#   geom_bar(aes(fill=Enrichment_Cell), stat="identity")+
#   theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
#   scale_fill_manual(values=c("#D64A56","cornflowerblue"))+
#   facet_wrap(~pathway, ncol=1)
# myeloid_GSEA_common

RR<-myeloid_res_overlap[which(myeloid_res_overlap$test=="RR"),]
KC<-myeloid_res_overlap[which(myeloid_res_overlap$test=="KC"),]
MHCII<-myeloid_res_overlap[which(myeloid_res_overlap$test=="MHCII"),]

three<-names(table(myeloid_res_overlap$pathway)[which(table(myeloid_res_overlap$pathway)>2)])
two<-intersect(KC$pathway, MHCII$pathway)[which(!(intersect(KC$pathway, MHCII$pathway)%in%three))]
two2<-intersect(RR$pathway, MHCII$pathway)[which(!(intersect(RR$pathway, MHCII$pathway)%in%c(three,two)))]
two3<-intersect(RR$pathway, KC$pathway)[which(!(intersect(RR$pathway, KC$pathway)%in%c(three,two, two2)))]
one<-myeloid_res_overlap[which(!(myeloid_res_overlap$pathway%in%c(three, two, two2, two3))),]
one<-as.data.frame(one %>% group_by(test) %>% slice_max(abs(NES), n = 5))
one<-one$pathway

myeloid_res_overlap_plt<-myeloid_res_overlap[which(myeloid_res_overlap$pathway%in%c(three, two, two2, two3, one)),]
myeloid_res_overlap_plt$pathway<-factor(myeloid_res_overlap_plt$pathway, levels=rev(c(three, two, two2, two3, one)))

myeloid_res_overlap_plt_unique<-myeloid_res_overlap_plt[,c("pathway","label")]
myeloid_res_overlap_plt_unique<-myeloid_res_overlap_plt_unique[!duplicated(myeloid_res_overlap_plt_unique$pathway),]

GSEA_all_myeloid<-plot_grid(ggplot(myeloid_res_overlap_plt, aes(test,pathway, fill=NES))+
                              geom_tile()+
                              th+theme_classic()+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8, 'RdBu')), space='Lab')(100), 
                                                                      name="",breaks = c(min(myeloid_res_overlap_plt$NES),
                                                                                         max(myeloid_res_overlap_plt$NES)),
                                                                      labels = c("Enriched\nPed",
                                                                                 "Enriched\nAdult"))+
                              ylab("")+xlab("")+ggtitle("gsea NES")+
                              geom_hline(yintercept = 5.5, color="grey")+
                              geom_hline(yintercept = 10.5, color="grey")+
                              geom_hline(yintercept = 15.5, color="grey")+
                              geom_hline(yintercept = 26.5, color="grey")+
                              geom_hline(yintercept = 31.5, color="grey")+
                              geom_hline(yintercept = 34.5, color="grey")+
                              theme(legend.position="top",legend.justification = "left", legend.title.align = 7,
                                    title = element_text(size=9)),
                            ggplot(myeloid_res_overlap_plt_unique, aes(0,pathway))+
                              geom_text(aes(label=label, hjust = 0), color="grey40", size=3)+
                              th+theme_classic()+theme(axis.text = element_blank(),
                                                       axis.line = element_line(colour="white"),
                                                       axis.ticks = element_line(colour="white"))+
                              ylab("")+xlab("")+xlim(0,1)+
                              geom_hline(yintercept = 5.5, color="grey")+
                              geom_hline(yintercept = 10.5, color="grey")+
                              geom_hline(yintercept = 15.5, color="grey")+
                              geom_hline(yintercept = 26.5, color="grey")+
                              geom_hline(yintercept = 31.5, color="grey")+
                              geom_hline(yintercept = 34.5, color="grey"),
                            align="h",axis="tb", rel_widths = c(2,1))
save_plts(GSEA_all_myeloid, "Myeloid_age_GSEA_heat", w=10.5,h=10)




###########
## Composition myeloid "types"
###########
myeloid_cluster_umap<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=T, group.by = "CellType_refined")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-7, y=-9, label=paste0("n = ",comma(ncol(d10x.combined_myeloid))))
myeloid_cluster_umap
save_plts(myeloid_cluster_umap, "IFALD_myeloid_map_celltype_subcluster_for_differential", w=7,h=6)

fancy_myeloid<-fanciest_UMAP(d10x.combined_myeloid,"KC Like",F)
save_plts(fancy_myeloid, "IFALD_myeloid_diff_KC_highlight", w=4,h=3)

fancy_myeloid<-fanciest_UMAP(d10x.combined_myeloid, NA,F)
save_plts(fancy_myeloid, "IFALD_myeloid_UMAP", w=4,h=3)

fancy_myeloid<-fanciest_UMAP(d10x.combined_myeloid, NA,T)
save_plts(fancy_myeloid, "IFALD_myeloid_UMAP_split", w=6,h=4)



table(d10x.combined_myeloid$CellType_refined, d10x.combined_myeloid$age_condition)

cell_counts<-d10x.combined_myeloid@meta.data %>% 
  group_by(individual, CellType_refined, age_condition,Age) %>% 
  summarise(count=length(unique(cell))) %>% 
  group_by(individual) %>%
  mutate(countT= sum(count)) %>%
  group_by(CellType_refined, add=TRUE) %>%
  mutate(per=100*count/countT)


cell_counts$label<-sapply(1:nrow(cell_counts), function(x){
  if(length(grep("NPC", cell_counts$individual[x]))==1){
    paste(cell_counts$Age[x], "\n(", strsplit(cell_counts$individual[x],"_")[[1]][1]," NPC)", sep="")
  }else{if(length(grep("TLH", cell_counts$individual[x])==1)){
    paste(cell_counts$Age[x], "\n(", strsplit(cell_counts$individual[x],"_")[[1]][1]," TLH)", sep="")
  }else{
    paste(cell_counts$Age[x], "\n(", strsplit(cell_counts$individual[x],"_")[[1]][1],")", sep="")}}
})

cell_counts$label<-factor(cell_counts$label, c("2\n(C104)","9\n(C105)","11\n(C85)","12\n(C93)", "16\n(C102)","17\n(C64)", "17\n(C96)",
                                               "0.33\n(IFALD030)","0.58\n(IFALD073)", "9\n(IFALD006)", 
                                               "26\n(C82)","48\n(C70)", "57\n(C97)","61\n(C68)",
                                               "65\n(C39 NPC)", "65\n(C39 TLH)", "67\n(C54)","69\n(C88)"))

cell_counts_min<-cell_counts[,c("individual","age_condition","Age","countT","label")][!duplicated(cell_counts[,c("individual","age_condition","Age","countT","label")]),]


cell_counts$CellType_refined<-factor(cell_counts$CellType_refined, levels=c("Macrophage\n(CLEC9A high)","Macrophage\n(MHCII high)",
                                                                            "KC Like","RR Myeloid","Myeloid Erythrocytes\n(phagocytosis)","Cycling Myeloid"))
myeloid_composistion<-ggplot(cell_counts, aes(label, per))+geom_bar(aes(fill=CellType_refined),stat = "identity", color="black")+
  theme_bw()+th+fillscale_cellType+xlab("Age\n(Sample ID)")+ylab("Percent of Cells in Sample")+
  facet_grid(.~age_condition, scale="free_x", space = "free")+
  geom_text(aes(label=countT, y=102), data=cell_counts_min)
myeloid_composistion
save_plts(myeloid_composistion, "myeloid_composistion", w=16,h=10)


myeloid_cluster_umap<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=F, group.by = "CellType_refined")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=10, y=-15, label=paste0("n = ",comma(ncol(d10x.combined_myeloid))))
myeloid_cluster_umap
save_plts(myeloid_cluster_umap, "IFALD_myeloid_map", w=7,h=5)

cell_num_myeloid<-as.data.frame(table(d10x.combined_myeloid$age_condition))
colnames(cell_num_myeloid)<-c("age_condition","CellCount")
myeloid_cluster_umap<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=F,split.by = "age_condition", group.by = "CellType_refined", ncol=2)+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  geom_text(aes(x=7, y=-15,label=paste0("n = ",comma(CellCount))),cell_num_myeloid, hjust=-0.1, size=3)
myeloid_cluster_umap
save_plts(myeloid_cluster_umap, "IFALD_myeloid_map_split", w=8,h=6)



###############
## plot individual genes
###############
myeloid_differential<-plot_grid(plot_gene_UMAP(d10x.combined_myeloid,"CCL3", 0),
                                plot_gene_UMAP(d10x.combined_myeloid,"CCL4", 0),
                                plot_gene_UMAP(d10x.combined_myeloid,"IL1B", 0),
                                plot_gene_UMAP(d10x.combined_myeloid,"HLA-DRB1", 0))
myeloid_differential
save_plts(myeloid_differential, "IFALD_myeloid_diff_genes_myeloidonly_tidy", w=7,h=5)


myeloid_differential_violin<-plot_grid(plot_gene_violin(d10x.combined_myeloid,"CCL3"),
                                       plot_gene_violin(d10x.combined_myeloid,"CCL4"),
                                       plot_gene_violin(d10x.combined_myeloid,"IL1B"),
                                       plot_gene_violin(d10x.combined_myeloid,"HLA-DRB1"))
myeloid_differential_violin
save_plts(myeloid_differential_violin, "IFALD_myeloid_diff_genes_myeloidonly_violin", w=8,h=9)



myeloid_differential_violin<-plot_grid(plot_gene_violin(d10x.combined_myeloid,"LY96"),
                                       plot_gene_violin(d10x.combined_myeloid,"CETP"),
                                       plot_gene_violin(d10x.combined_myeloid,"CCL4"),
                                       plot_gene_violin(d10x.combined_myeloid,"CTSB"))
myeloid_differential_violin
save_plts(myeloid_differential_violin, "IFALD_myeloid_diff_genes_withIFALD_myeloidonly_violin", w=8,h=9)



myeloid_differential_violin<-plot_grid(plot_gene_UMAP(d10x.combined_myeloid,"LY96",0),
                                       plot_gene_UMAP(d10x.combined_myeloid,"CETP",0),
                                       plot_gene_UMAP(d10x.combined_myeloid,"CCL4",0),
                                       plot_gene_UMAP(d10x.combined_myeloid,"CTSB",0))
myeloid_differential_violin
save_plts(myeloid_differential_violin, "IFALD_myeloid_diff_genes_withIFALD_myeloidonly_umap", w=8,h=9)



myeloid_differential<-plot_grid(plot_gene_UMAP(d10x.combined_myeloid,"LYVE1", 0.8),
                                plot_gene_UMAP(d10x.combined_myeloid,"CD9", 0.8),
                                plot_gene_UMAP(d10x.combined_myeloid,"CLEC10A", 0.8),
                                plot_gene_UMAP(d10x.combined_myeloid,"HLA-DRB1", 0.8))
myeloid_differential
save_plts(myeloid_differential, "IFALD_myeloid_diff_genes_myeloidonly_tidy_MHCII_identity", w=7,h=5)



#############
## Heat Map of differential genes
#############
d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))

load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

d10x.combined_myeloid<-subset(d10x, subset = CellType_refined %in% c("RR Myeloid","Macrophage\n(MHCII high)","KC Like","Macrophage\n(CLEC9A high)","Cycling Myeloid","Myeloid Erythrocytes\n(phagocytosis)"))




# cluster0 age
kc_age<-c("ALB","SAA1","CCL3","CCL4","IL1B")
# cluster0 IFALD
kc_ifald<-c("CD5L","A2M","LY96")

# RR age
rr_age<-c("AREG","FOS")
# RR IFALD
rr_ifald<-c("S100A8","S100A9")

# MHCII age
mhcII_age<-c("HMOX1", "APOE")
# MHCII IFALD
mhcII_ifald<-c("VSIG4","CD163","HLA-DPA1","HLA-DRA")

### for adding *
sig_de_IFALD_RR<-read.csv(file=here("data","differential_IFALD_RR.csv"))
sig_de_IFALD_RR$age_condition<-"Ped\nIFALD"
sig_de_IFALD_RR$CellType_refined<-"RR Myeloid"

sig_de_IFALD_KC<-read.csv(file=here("data","differential_IFALD_KC.csv"))
sig_de_IFALD_KC$age_condition<-"Ped\nIFALD"
sig_de_IFALD_KC$CellType_refined<-"KC Like"

sig_de_IFALD_MHCII<-read.csv(file=here("data","differential_IFLAD_MHCII.csv"))
sig_de_IFALD_MHCII$age_condition<-"Ped\nIFALD"
sig_de_IFALD_MHCII$CellType_refined<-"Macrophage\n(MHCII high)"

sig_de_age_RR<-read.csv(file=here("data","differential_age_RR.csv"))
sig_de_age_RR$age_condition<-"Adult\nHealthy"
sig_de_age_RR$CellType_refined<-"RR Myeloid"

sig_de_age_KC<-read.csv(file=here("data","differential_age_KC.csv"))
sig_de_age_KC$age_condition<-"Adult\nHealthy"
sig_de_age_KC$CellType_refined<-"KC Like"

sig_de_age_MHCII<-read.csv(file=here("data","differential_age_MHCII.csv"))
sig_de_age_MHCII$age_condition<-"Adult\nHealthy"
sig_de_age_MHCII$CellType_refined<-"Macrophage\n(MHCII high)"

de<-rbind(sig_de_IFALD_RR, sig_de_IFALD_KC, sig_de_IFALD_MHCII, sig_de_age_RR, sig_de_age_KC, sig_de_age_MHCII)
colnames(de)[1]<-"variable"
de$label<-"*"


myeloid_age_heat<-plot_heat_map(d10x.combined_myeloid,c(kc_age, mhcII_age,rr_age), 
                                c("KC Like","Macrophage\n(MHCII high)","RR Myeloid"),T)
save_plts(myeloid_age_heat, "Myeloid_age_heat", w=7,h=5)

myeloid_IFALD_heat<-plot_heat_map(d10x.combined_myeloid,c(kc_ifald,mhcII_ifald,rr_ifald ), 
                                  c("KC Like","Macrophage\n(MHCII high)","RR Myeloid"),T)
save_plts(myeloid_IFALD_heat, "Myeloid_IFALD_heat", w=7,h=4)



