### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(reshape2)
library(gtools)


source("R_functions/pretty_plots.R")


load(here("data","adult_ped_integrated.rds"))

d10x.combined_myeloid<-subset(d10x.combined, subset = CellType_rough %in% c("Myeloid"))
rm(d10x.combined)
gc()
d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)
d10x.combined_myeloid <- RunUMAP(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindNeighbors(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindClusters(d10x.combined_myeloid, resolution = 0.1)
DimPlot(d10x.combined_myeloid, label=T)

myeliod_clusters<-d10x.combined_myeloid@meta.data

## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","d10x_adult_ped_raw.rds"))
d10x_raw_myeloid<-subset(d10x, cells = rownames(myeliod_clusters))
rm(d10x)
gc()

myeliod_clusters$index<-rownames(myeliod_clusters)
identical(colnames(d10x_raw_myeloid), myeliod_clusters$index)

d10x_raw_myeloid <- AddMetaData(d10x_raw_myeloid, metadata = myeliod_clusters)



d10x_raw_myeloid <- NormalizeData(d10x_raw_myeloid,scale.factor = 10000, normalization.method = "LogNormalize")


## testing factor
Idents(d10x_raw_myeloid) <- "seurat_clusters"
table(d10x_raw_myeloid$seurat_clusters)


#MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene,
#consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
#and linear regression for the continuous process (i.e., the expression level). 

myeloid_clusters<-unique(d10x_raw_myeloid$seurat_clusters)

diff_exp_all<-lapply(1:length(myeloid_clusters), function(x){
  de<-FindMarkers(d10x_raw_myeloid, ident.1 = myeloid_clusters[x], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
  de$gene<-rownames(de)
  rownames(de)<-NULL
  de<-de[,c(6,1:5)]
  de$cluster<-myeloid_clusters[x]
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

### plot markers
FeaturePlot(d10x.combined_myeloid, features = "DEFA3", min.cutoff = "q9", pt.size=1)

VlnPlot(d10x.combined_myeloid, features = "DEFA3", pt.size = 0, log=T)
VlnPlot(d10x.combined_myeloid, features = "DEFA3", pt.size = 0, log=F)

recent_recruit_myeloid<-c("S100A8","S100A9","CD68","LYZ")
kuffer_signature<-c("VSIG4","MARCO","CD5L","HMOX1")

top_DE_up[which(top_DE_up$gene%in%recent_recruit_myeloid),]
top_DE_up[which(top_DE_up$gene%in%kuffer_signature),]

FeaturePlot(d10x.combined_myeloid, features = kuffer_signature, min.cutoff = "q9", pt.size=1)
FeaturePlot(d10x.combined_myeloid, features = recent_recruit_myeloid, min.cutoff = "q9", pt.size=1)


###################################################
#0: Low quality, soupy, high MT
top_DE_up[which(top_DE_up$cluster=="0"),]
top_DE[which(top_DE$cluster=="0"),]
FeaturePlot(d10x.combined_myeloid, features = c("MT-CO2","MT-CO3","MT-ATP6","MT-ND4"), min.cutoff = "q9", pt.size=1)

#1: unclear
top_DE[which(top_DE$cluster=="1"),]
top_DE_up[which(top_DE_up$cluster=="1"),]
FeaturePlot(d10x.combined_myeloid, features = c("LYZ","CD300E","FCN1","CD5L"), min.cutoff = "q9", pt.size=1)


# 3,4,7 kupffer-like
top_DE_up[which(top_DE_up$gene%in%kuffer_signature),]
myeloid_347<-FeaturePlot(d10x.combined_myeloid, features = kuffer_signature, min.cutoff = "q9", pt.size=1)
myeloid_347
save_plts(myeloid_347, "myeloid_cluster_347_markers", w=8,h=6)

#3: maybe kupffer-like
top_DE[which(top_DE$cluster=="3"),]
top_DE_up[which(top_DE_up$cluster=="3"),]
FeaturePlot(d10x.combined_myeloid, features = c("C1QB","ARL4C","CCL3","CCL4"), min.cutoff = "q9", pt.size=1)

#4: maybe kupffer-like 
top_DE[which(top_DE$cluster=="4"),]
top_DE_up[which(top_DE_up$cluster=="4"),]
FeaturePlot(d10x.combined_myeloid, features = c("MARCO","SDC3","SLC40A1","VCAM1"), min.cutoff = "q9", pt.size=1)

#7: maybe kupffer-like (with hemaglobin?)
top_DE[which(top_DE$cluster=="7"),]
top_DE_up[which(top_DE_up$cluster=="7"),]
myeloid_7<-FeaturePlot(d10x.combined_myeloid, features = c("HBA2","HBA1","HBB","MARCO"), min.cutoff = "q9", pt.size=1)
myeloid_7
save_plts(myeloid_7, "myeloid_cluster_7_markers", w=8,h=6)



#5: B cells, almost all markers sig up
top_DE[which(top_DE$cluster=="5"),]
top_DE_up[which(top_DE_up$cluster=="5"),]
B_genes<-c("POU2F2","FCER2","MS4A1","LTB","CD37","CD79B","IGLC2","IGHG1","IGKC", "CD19")
myeloid_5<-FeaturePlot(d10x.combined_myeloid, features = c("POU2F2","CD37","CD19","CD79B"), min.cutoff = "q9", pt.size=1)
myeloid_5
save_plts(myeloid_5, "myeloid_cluster_5_markers", w=8,h=6)
FeaturePlot(d10x.combined_myeloid, features = B_genes, min.cutoff = "q9", pt.size=1)
diff_exp_sig[which(diff_exp_sig$gene%in%B_genes),]

#6: unclear maybe neutrophils too? (DEFA3, DEFA4 specific)
top_DE[which(top_DE$cluster=="6"),]
top_DE_up[which(top_DE_up$cluster=="6"),]
myeloid_6<-FeaturePlot(d10x.combined_myeloid, features = c("DEFA3","DEFA4","OLFM4","AZU1"), min.cutoff = "q9", pt.size=1)
myeloid_6
save_plts(myeloid_6, "myeloid_cluster_6_markers", w=8,h=6)

#2: Neutrophils
top_DE[which(top_DE$cluster=="2"),]
top_DE_up[which(top_DE_up$cluster=="2"),]
FeaturePlot(d10x.combined_myeloid, features = c("S100P","CMTM2","FCGR3B","IFITM2"), min.cutoff = "q9", pt.size=1)
#neutrophil makers from literature
neutro_gene<-c("CSF3R","FCGR3B","NAMPT","CXCR2")
diff_exp_sig[which(diff_exp_sig$gene%in%neutro_gene),]
myeloid_2<-FeaturePlot(d10x.combined_myeloid, features = neutro_gene, min.cutoff = "q9", pt.size=1)
myeloid_2
save_plts(myeloid_2, "myeloid_cluster_2_markers", w=8,h=6)



