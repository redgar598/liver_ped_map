#'---
#'title: scRNAseq signatures
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
library(ggsignif)
library(stats)
library(scales)



options(stringsAsFactors = FALSE)

source("R_functions/pretty_plots.R")



### load integrate for UMAP etc
load(here("data","adult_ped_integrated.rds"))


B_genes<-c("POU2F2","FCER2","MS4A1","LTB","CD37","CD79B","IGLC2","IGHG1","IGKC", "CD19")
T_genes<-c("CD3D","IL7R","CD8A","IL32")
NK_genes<-c("NKG7","CD7")
gd_genes<-c("GNLY")

d10x.combined_NK_T_B<-subset(d10x.combined, subset = CellType_rough %in% c("CD3_Tcell","nkTcell","gdTcell"))
rm(d10x.combined)
gc()
d10x.combined_NK_T_B <- RunPCA(d10x.combined_NK_T_B, npcs = 30, verbose = FALSE)
d10x.combined_NK_T_B <- RunUMAP(d10x.combined_NK_T_B, reduction = "pca", dims = 1:30)
d10x.combined_NK_T_B <- FindNeighbors(d10x.combined_NK_T_B, reduction = "pca", dims = 1:30)
d10x.combined_NK_T_B <- FindClusters(d10x.combined_NK_T_B, resolution = 0.1)
NK_T_B_umap<-DimPlot(d10x.combined_NK_T_B, label=T)
NK_T_B_umap
save_plts(NK_T_B_umap, "NK_T_B_umap", w=8,h=6)

NK_B_T_clusters<-d10x.combined_NK_T_B@meta.data



FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = NK_genes, ncol = 2)
FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = gd_genes)
FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = T_genes, ncol = 2)
FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = B_genes, ncol = 2)


## B cell hunt
FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = c("MS4A1"), ncol = 2)




#### DE for cluster spefici genes
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","d10x_adult_ped_raw.rds"))
d10x_raw_NK_T_B<-subset(d10x, cells = rownames(NK_B_T_clusters))
rm(d10x)
gc()

NK_B_T_clusters$index<-rownames(NK_B_T_clusters)
identical(colnames(d10x_raw_NK_T_B), NK_B_T_clusters$index)

d10x_raw_NK_T_B <- AddMetaData(d10x_raw_NK_T_B, metadata = NK_B_T_clusters)

d10x_raw_NK_T_B <- NormalizeData(d10x_raw_NK_T_B,scale.factor = 10000, normalization.method = "LogNormalize")



## testing factor
Idents(d10x_raw_NK_T_B) <- "seurat_clusters"
table(d10x_raw_NK_T_B$seurat_clusters)


#MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene,
#consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
#and linear regression for the continuous process (i.e., the expression level). 

NK_T_B_clusters<-unique(d10x_raw_NK_T_B$seurat_clusters)

diff_exp_all<-lapply(1:length(NK_T_B_clusters), function(x){
  de<-FindMarkers(d10x_raw_NK_T_B, ident.1 = NK_T_B_clusters[x], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
  de$gene<-rownames(de)
  rownames(de)<-NULL
  de<-de[,c(6,1:5)]
  de$cluster<-NK_T_B_clusters[x]
  de})


diff_exp_all<-do.call(rbind, diff_exp_all)

diff_exp_sig<-diff_exp_all[which(diff_exp_all$p_val_adj<0.005),]



# ### Top DE genes
diff_exp_sig %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = abs(avg_log2FC)) -> top10

top_DE<-as.data.frame(top10)

diff_exp_sig %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top_DE_up

top_DE_up<-as.data.frame(top_DE_up)

top_DE_up[which(top_DE_up$gene%in%c(T_genes,NK_genes,gd_genes)),]

###################################################
#0: CD3+ CD8+ Tcell
top_DE_up[which(top_DE_up$cluster=="0"),]
top_DE[which(top_DE$cluster=="0"),]
NK_T_N_0<-FeaturePlot(d10x.combined_NK_T_B, features = c("S100A4","IL7R","CD8A","CD3D"), min.cutoff = "q9", pt.size=1)
NK_T_N_0
save_plts(NK_T_N_0, "NK_T_N_cluster_0_markers", w=8,h=6)


#1: CD7+
top_DE[which(top_DE$cluster=="1"),]
top_DE_up[which(top_DE_up$cluster=="1"),]
NK_T_N_1<-FeaturePlot(d10x.combined_NK_T_B, features = c("CD7","FCER1G","IL2RB","CCL3"), min.cutoff = "q9", pt.size=1)
NK_T_N_1
save_plts(NK_T_N_1, "NK_T_N_cluster_1_markers", w=8,h=6)

#2: GNLY+
top_DE[which(top_DE$cluster=="2"),]
top_DE_up[which(top_DE_up$cluster=="2"),]
NK_T_N_2<-FeaturePlot(d10x.combined_NK_T_B, features = c("GNLY","CX3CR1","FCGR3A","FGFBP2"), min.cutoff = "q9", pt.size=1)
NK_T_N_2
save_plts(NK_T_N_2, "NK_T_N_cluster_2_markers", w=8,h=6)

#2: GNLY+
top_DE[which(top_DE$cluster=="2"),]
top_DE_up[which(top_DE_up$cluster=="2"),]
NK_T_N_2<-FeaturePlot(d10x.combined_NK_T_B, features = c("GNLY","CX3CR1","FCGR3A","FGFBP2"), min.cutoff = "q9", pt.size=1)
NK_T_N_2
save_plts(NK_T_N_2, "NK_T_N_cluster_2_markers", w=8,h=6)

#3: Soupy? ALB SERPINA1 SAA1
top_DE[which(top_DE$cluster=="3"),]
top_DE_up[which(top_DE_up$cluster=="3"),]
NK_T_N_3<-FeaturePlot(d10x.combined_NK_T_B, features = c("ALB","SERPINA1","HP","SAA1"), min.cutoff = "q9", pt.size=1)
NK_T_N_3
save_plts(NK_T_N_3, "NK_T_N_cluster_3_markers", w=8,h=6)

#4: CD8+ IL32+
top_DE[which(top_DE$cluster=="4"),]
top_DE_up[which(top_DE_up$cluster=="4"),]
NK_T_N_4<-FeaturePlot(d10x.combined_NK_T_B, features = c("CD8A","IL32","CD74","CD8B"), min.cutoff = "q9", pt.size=1)
NK_T_N_4
save_plts(NK_T_N_4, "NK_T_N_cluster_4_markers", w=8,h=6)

de_04<-FindMarkers(d10x_raw_NK_T_B, ident.1 = "0", ident.2= "4", test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
de_04_sig<-de_04[which(de_04$p_val_adj < 0.005 & abs(de_04$avg_log2FC) > 1),]
de_04_sig[which(de_04_sig$avg_log2FC<0),]

FeaturePlot(d10x.combined_NK_T_B, features = c("NKG7","CCL5","HLA-A","MT-ND4"), min.cutoff = "q9", pt.size=1)


