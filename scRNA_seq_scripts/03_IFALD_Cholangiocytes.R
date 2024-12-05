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


source("scRNA_seq_scripts/00_pretty_plots.R")
source("scRNA_seq_scripts/00_fanciest_UMAP.R")
source("scRNA_seq_scripts/00_plot_gene_exp.R")
source("scRNA_seq_scripts/00_entropy_d10x.R")



load(here("data","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

d10x.combined_cholangiocytes<-subset(d10x.combined, subset = CellType_refined %in% c("Cholangiocytes"))
rm(d10x.combined)
gc()
d10x.combined_cholangiocytes <- RunPCA(d10x.combined_cholangiocytes, npcs = 30, verbose = FALSE)
d10x.combined_cholangiocytes <- RunUMAP(d10x.combined_cholangiocytes, reduction = "pca", dims = 1:30)
d10x.combined_cholangiocytes <- FindNeighbors(d10x.combined_cholangiocytes, reduction = "pca", dims = 1:30)
d10x.combined_cholangiocytes <- FindClusters(d10x.combined_cholangiocytes, resolution = 0.3)

cholangiocytes_subtype<-DimPlot(d10x.combined_cholangiocytes, label=T)
cholangiocytes_subtype
save_plts(cholangiocytes_subtype, "IFALD_cholangiocytes_map_clusters", w=7,h=6)

cholangiocytes_cluster_umap<-DimPlot(d10x.combined_cholangiocytes, reduction = "umap", pt.size=0.25, label=T, group.by = "CellType_refined")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-7, y=-9, label=paste0("n = ",comma(ncol(d10x.combined_cholangiocytes))))
cholangiocytes_cluster_umap
save_plts(cholangiocytes_cluster_umap, "IFALD_cholangiocytes_map_celltype", w=7,h=6)

cell_num_cholangiocytes<-as.data.frame(table(d10x.combined_cholangiocytes$age_condition))
colnames(cell_num_cholangiocytes)<-c("age_condition","CellCount")
cholangiocytes_cluster_umap<-DimPlot(d10x.combined_cholangiocytes, reduction = "umap", pt.size=0.25, label=F,split.by = "age_condition", group.by = "CellType_refined", ncol=2)+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  geom_text(aes(y=-7, x=-9,label=paste0("n = ",comma(CellCount))),cell_num_cholangiocytes, hjust=-0.1, size=3)
cholangiocytes_cluster_umap
save_plts(cholangiocytes_cluster_umap, "IFALD_cholangiocytes_map", w=7,h=6)

cholangiocytes_cluster_umap<-DimPlot(d10x.combined_cholangiocytes, reduction = "umap", pt.size=0.25, label=F,split.by = "age_id", group.by = "CellType_refined", ncol=4)+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-14, y=-12, label=paste0("n = ",comma(ncol(d10x.combined_cholangiocytes))))
cholangiocytes_cluster_umap
save_plts(cholangiocytes_cluster_umap, "IFALD_cholangiocytes_map_individual", w=12,h=10)

plt_entropy_age_condition<-entropy_d10(d10x.combined_cholangiocytes, "age_condition")
entropy_cholangiocytesage<-entropy_plt(plt_entropy_age_condition, "age_condition", d10x.combined_cholangiocytes)
entropy_cholangiocytesage
save_plts(entropy_cholangiocytesage, "entropy_age_cholangiocytes", w=15,h=10)

##############
## Differential expression with age 
##############
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))

load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
d10x$Sex[which(d10x$individual%in%c("C113","C115"))]<-"M"

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")


d10x_raw_Cholangiocytes<-subset(d10x, subset = CellType_refined %in% c("Cholangiocytes"))

Idents(d10x_raw_Cholangiocytes)<-d10x_raw_Cholangiocytes$age_condition
table(d10x_raw_Cholangiocytes$age_condition)

## age differential
de<-FindMarkers(d10x_raw_Cholangiocytes, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]
sig_de[which(sig_de$avg_log2FC>0),]
sig_de[which(sig_de$avg_log2FC<0),]

### IFALD differential
de_IFALD<-FindMarkers(d10x_raw_Cholangiocytes, ident.1 = "Ped IFALD", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de_IFALD<-de_IFALD[which(de_IFALD$p_val_adj < 0.005 & abs(de_IFALD$avg_log2FC) > 1),]
sig_de_IFALD[which(sig_de_IFALD$avg_log2FC>0),]
sig_de_IFALD[which(sig_de_IFALD$avg_log2FC<0),]

write.csv(sig_de, file=here("data","differential_age_Cholangiocytes.csv"))
write.csv(sig_de_IFALD, file=here("data","differential_IFALD_Cholangiocytes.csv"))





###
## pathway adult/IFALD versus healthy ped
###
source("scRNA_seq_scripts/00_GSEA_function.R")
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

cholangiocytes_GSEA<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_hline(yintercept=16.5, color="grey")+scale_fill_manual(values=c("#D64A56","cornflowerblue"))
cholangiocytes_GSEA
save_plts(cholangiocytes_GSEA, "GSEA_Age_cholangiocytes", w=15,h=7)

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

cholangiocytes_GSEA<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_hline(yintercept=16.5, color="grey")+scale_fill_manual(values=c("cornflowerblue","#374eb8"))
cholangiocytes_GSEA
save_plts(cholangiocytes_GSEA, "GSEA_IFALD_cholangiocytes", w=15,h=7)
