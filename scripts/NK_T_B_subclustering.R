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
library(colorspace)




options(stringsAsFactors = FALSE)

source("scripts/00_pretty_plots.R")



### load integrate for UMAP etc
load(here("data","adult_ped_integrated_refinedlabels.rds"))


B_genes<-c("POU2F2","FCER2","MS4A1","LTB","CD37","CD79B","IGLC2","IGHG1","IGKC", "CD19")
T_genes<-c("CD3D","IL7R","CD8A","IL32")
NK_genes<-c("NKG7","CD7")
gd_genes<-c("GNLY")

d10x.combined_NK_T<-subset(d10x.combined, subset = CellType_refined %in% c("CD3+ T-cells","NK-like cells","gd T-cells"))
rm(d10x.combined)
gc()
d10x.combined_NK_T <- RunPCA(d10x.combined_NK_T, npcs = 30, verbose = FALSE)
d10x.combined_NK_T <- RunUMAP(d10x.combined_NK_T, reduction = "pca", dims = 1:30)
d10x.combined_NK_T <- FindNeighbors(d10x.combined_NK_T, reduction = "pca", dims = 1:30)
d10x.combined_NK_T <- FindClusters(d10x.combined_NK_T, resolution = 0.1)
DimPlot(d10x.combined_NK_T)


FeaturePlot(d10x.combined_NK_T, reduction = "umap", features = NK_genes, ncol = 2)
FeaturePlot(d10x.combined_NK_T, reduction = "umap", features = gd_genes)
FeaturePlot(d10x.combined_NK_T, reduction = "umap", features = T_genes, ncol = 2)


# are adult T cells more hepatocyte
age_tcell<-DimPlot(d10x.combined_NK_T, split.by = "AgeGroup", group.by = "CellType_refined")+colscale_cellType

Hep_Tcell<-FeaturePlot(d10x.combined_NK_T, reduction = "umap", features = c("ALB","SERPINA1","HAMP","PCK1"), ncol = 2)
Tcell_hep<-plot_grid(age_tcell, Hep_Tcell, ncol=1, rel_heights = c(1.1,2))
save_plts(Tcell_hep, "age_tcell_hepmarker", w=8,h=10)



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
NK_T_N_0<-FeaturePlot(d10x.combined_NK_T, features = c("CD52","IL7R","CD8A","CD3D"), min.cutoff = "q9", pt.size=1)
NK_T_N_0
save_plts(NK_T_N_0, "NK_T_N_cluster_0_markers", w=8,h=6)


#1: Adult Hep NK? Soupy
top_DE[which(top_DE$cluster=="1"),]
top_DE_up[which(top_DE_up$cluster=="1"),]
NK_T_N_1<-FeaturePlot(d10x.combined_NK_T, features = c("ALB","APOA1","APOA2","SERPINA1"), min.cutoff = "q9", pt.size=1)
NK_T_N_1
save_plts(NK_T_N_1, "NK_T_N_cluster_1_markers", w=8,h=6)

#2: NK like
top_DE[which(top_DE$cluster=="2"),]
top_DE_up[which(top_DE_up$cluster=="2"),]
NK_T_N_2<-FeaturePlot(d10x.combined_NK_T, features = c("FCER1G","FOS","CCL3","IL2RB"), min.cutoff = "q9", pt.size=1)
NK_T_N_2
save_plts(NK_T_N_2, "NK_T_N_cluster_2_markers", w=8,h=6)

#3: GNLY+
top_DE[which(top_DE$cluster=="3"),]
top_DE_up[which(top_DE_up$cluster=="3"),]
NK_T_N_3<-FeaturePlot(d10x.combined_NK_T, features = c("GNLY","FCGR3A","GZMB","PRF1"), min.cutoff = "q9", pt.size=1)
NK_T_N_3
save_plts(NK_T_N_3, "NK_T_N_cluster_3_markers", w=8,h=6)

#4: Adult Hep NK? Soupy
top_DE[which(top_DE$cluster=="4"),]
top_DE_up[which(top_DE_up$cluster=="4"),]
NK_T_N_4<-FeaturePlot(d10x.combined_NK_T, features = c("ALB","SAA1","SAA2","CRP"), min.cutoff = "q9", pt.size=1)
NK_T_N_4
save_plts(NK_T_N_4, "NK_T_N_cluster_4_markers", w=8,h=6)


#5: CD8A TRAV17+ (alpha)
top_DE[which(top_DE$cluster=="5"),]
top_DE_up[which(top_DE_up$cluster=="5"),]
NK_T_N_5<-FeaturePlot(d10x.combined_NK_T, features = c("CD8A","CD74","TRAV17","IL32"), min.cutoff = "q9", pt.size=1)
NK_T_N_5
save_plts(NK_T_N_5, "NK_T_N_cluster_5_markers", w=8,h=6)

#5: CD8A TRDV2+ (delta)
top_DE[which(top_DE$cluster=="6"),]
top_DE_up[which(top_DE_up$cluster=="6"),]
NK_T_N_6<-FeaturePlot(d10x.combined_NK_T, features = c("TRDV2","TRGV9","IL7R","SLC4A10"), min.cutoff = "q9", pt.size=1)
NK_T_N_6
save_plts(NK_T_N_6, "NK_T_N_cluster_6_markers", w=8,h=6)


de_04<-FindMarkers(d10x_raw_NK_T_B, ident.1 = "1", ident.2= "4", test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
de_04_sig<-de_04[which(de_04$p_val_adj < 0.005 & abs(de_04$avg_log2FC) > 1),]
de_04_sig[which(de_04_sig$avg_log2FC<0),]
de_04_sig[which(de_04_sig$avg_log2FC>0),]


FeaturePlot(d10x.combined_NK_T, features = c("TRAC","CD3D","FCER1G","CD7"), min.cutoff = "q9", pt.size=1)

FeaturePlot(d10x.combined_NK_T, features = c("NKG7","CD8A","CD3D","CD4"), min.cutoff = "q9", pt.size=1)

## MacParland markers for tcell (fig9)
FeaturePlot(d10x.combined_NK_T, features = c("PTPRC","CD2","CD3D","TRAC"), min.cutoff = "q9", pt.size=1)
FeaturePlot(d10x.combined_NK_T, features = c("IL7R","KLRB1","NKG7","FCGR3A"), min.cutoff = "q9", pt.size=1)
FeaturePlot(d10x.combined_NK_T, features = c("GZMA","GZMB","GZMK","PRF1"), min.cutoff = "q9", pt.size=1)
FeaturePlot(d10x.combined_NK_T, features = c("CD79A","CD79B","CD27","IGHG1"), min.cutoff = "q9", pt.size=1)
FeaturePlot(d10x.combined_NK_T, features = c("MS4A1","LTB","CD52","IGHD"), min.cutoff = "q9", pt.size=1)



#######################
## Naive vs memory T cells
#######################

FeaturePlot(d10x.combined_NK_T, features = c("CD4","CD8A","FOXP3","PDCD1"), min.cutoff = "q9", pt.size=1)

d10x.combined_CD3<-subset(d10x.combined_NK_T, subset = CellType_rough %in% c("CD3_Tcell"))
d10x.combined_CD3 <- RunPCA(d10x.combined_CD3, npcs = 30, verbose = FALSE)
d10x.combined_CD3 <- RunUMAP(d10x.combined_CD3, reduction = "pca", dims = 1:30)
d10x.combined_CD3 <- FindNeighbors(d10x.combined_CD3, reduction = "pca", dims = 1:30)
d10x.combined_CD3 <- FindClusters(d10x.combined_CD3, resolution = 0.1)
CD3_umap<-DimPlot(d10x.combined_CD3, label=T)
CD3_umap
save_plts(CD3_umap, "CD3_umap", w=5,h=4)

DimPlot(d10x.combined_CD3, group.by = "individual", label=T)

FeaturePlot(d10x.combined_CD3, features = c("GZMB","GZMK","CXCR5","PDCD1"), min.cutoff = "q9", pt.size=1)
FeaturePlot(d10x.combined_CD3, features = c("TCF7","CTLA4","CD4","HAVCR2"), min.cutoff = "q9", pt.size=1)
FeaturePlot(d10x.combined_CD3, features = c("CD4","CD8A","SELL","CCR7"), min.cutoff = "q9", pt.size=1)

genes<-c("GZMB","GZMK","CXCR5","PDCD1","TCF7","CTLA4","CD4","CD8A", "HAVCR2")

FeaturePlot(d10x.combined_CD3, features = c("CD8A", "CD4"), blend = TRUE)


d10x.exp<-as.data.frame(d10x.combined_CD3[["RNA"]]@data)
d10x.exp.GOI<-d10x.exp[genes,]
d10x.exp.GOI$gene<-rownames(d10x.exp.GOI)
d10x.exp.GOI<-melt(d10x.exp.GOI)#

umap_mat<-as.data.frame(Embeddings(object = d10x.combined_CD3, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x.combined_CD3@meta.data
meta$cell<-rownames(meta)

plt<-merge(meta, umap_mat, by="cell")

plt<-merge(d10x.exp.GOI, plt,by.x="variable", by.y="cell")
plt$variable<-as.character(plt$variable)
plt$seurat_clusters<-as.character(plt$seurat_clusters)

plt_list<-lapply(genes, function(x){
  plt_gene<-plt[which(plt$gene==x),]
  ggplot(plt_gene, aes(UMAP_1, UMAP_2, color=value))+geom_point(size=0.75)+facet_wrap(~gene)+
    scale_color_continuous_sequential(palette = "Blues 2") +theme_bw() + th_present
})

plot_grid(plotlist = plt_list, align = 'hv', ncol = 3)





#######################
## Naive vs memory T cells (split by individual)
#######################
d10x.combined_CD3<-subset(d10x.combined_NK_T, subset = CellType_rough %in% c("CD3_Tcell"))
d10x.combined_CD3 <- RunPCA(d10x.combined_CD3, npcs = 30, verbose = FALSE)
d10x.combined_CD3 <- RunUMAP(d10x.combined_CD3, reduction = "pca", dims = 1:30)
d10x.combined_CD3 <- FindNeighbors(d10x.combined_CD3, reduction = "pca", dims = 1:30)
d10x.combined_CD3 <- FindClusters(d10x.combined_CD3, resolution = 0.1)
CD3_umap<-DimPlot(d10x.combined_CD3, label=T)
CD3_umap
save_plts(CD3_umap, "CD3_umap", w=5,h=4)

DimPlot(d10x.combined_CD3, split.by  = "individual", label=T)

d10x.combined_CD3$age_id<-paste(d10x.combined_CD3$individual, d10x.combined_CD3$AgeGroup)
DimPlot(d10x.combined_CD3, split.by  = "age_id", label=T)


FeaturePlot(d10x.combined_CD3, features = c("GZMB","GZMK","CXCR5","PDCD1"), min.cutoff = "q9", pt.size=1)
FeaturePlot(d10x.combined_CD3, features = c("TCF7","CTLA4","CD4","HAVCR2"), min.cutoff = "q9", pt.size=1)
FeaturePlot(d10x.combined_CD3, features = c("CD4","CD8A","SELL","CCR7"), min.cutoff = "q9", pt.size=1)

genes<-c("GZMB","GZMK","CXCR5","PDCD1","TCF7","CTLA4","CD4","CD8A", "HAVCR2")

FeaturePlot(d10x.combined_CD3, features = c("CD8A", "CD4"), blend = TRUE)





######
#TRM interesting markers
######

#klrg1/cd57
#cd69/cd103

#CD57 is encoded by  B3GAT1
#CD103 is encoded by ITGAE

TRM_markers_all<-FeaturePlot(d10x.combined, reduction = "umap", features = c("KLRG1", "B3GAT1", "CD69","ITGAE"), ncol = 2)
save_plts(TRM_markers_all, "TRM_markers", w=5,h=4)

d10x.combined_NK_T<-subset(d10x.combined, subset = CellType_rough %in% c("CD3+ T-cells","gd T-cells","NK-like cells"))
d10x.combined_NK_T <- RunPCA(d10x.combined_NK_T, npcs = 30, verbose = FALSE)
d10x.combined_NK_T <- RunUMAP(d10x.combined_NK_T, reduction = "pca", dims = 1:30)

FeaturePlot(d10x.combined_NK_T, reduction = "umap", features = c("KLRG1", "B3GAT1", "CD69","ITGAE"), ncol = 2)
TRM_markers_NKTB<-FeaturePlot(d10x.combined_NK_T, reduction = "umap", features = c("KLRG1", "B3GAT1", "CD69","ITGAE"), split="AgeGroup", ncol = 2)
save_plts(TRM_markers_NKTB, "TRM_markers_NK_T_B", w=5,h=10)

d10x.exp<-as.data.frame(d10x.combined_NK_T[["RNA"]]@data)
d10x.exp.GOI<-d10x.exp[c("KLRG1", "B3GAT1", "CD69","ITGAE"),]

d10x.exp.GOI$gene<-rownames(d10x.exp.GOI)
d10x.exp.GOI<-melt(d10x.exp.GOI)#

meta<-d10x.combined_NK_T@meta.data
meta$cell<-rownames(meta)

plt<-merge(d10x.exp.GOI, meta,by.x="variable", by.y="cell")

quant_markers<-ggplot(plt, aes(AgeGroup, value))+geom_violin(aes(fill=AgeGroup))+facet_wrap(~gene)+fillscale_age+theme_bw()+th
save_plts(quant_markers, "TRM_markers_NK_T_B_quantified", w=5,h=4)


