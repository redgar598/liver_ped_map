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
myeloid_0<-FeaturePlot(d10x.combined_myeloid, features = c("MT-CO2","MT-CO3","MT-ATP6","MT-ND4"), min.cutoff = "q9", pt.size=1)
myeloid_0
save_plts(myeloid_0, "myeloid_cluster_0_markers", w=8,h=6)


#1: unclear
top_DE[which(top_DE$cluster=="1"),]
top_DE_up[which(top_DE_up$cluster=="1"),]
myeloid_1<-FeaturePlot(d10x.combined_myeloid, features = c("LYZ","CD300E","FCN1","EREG"), min.cutoff = "q9", pt.size=1)
myeloid_1
save_plts(myeloid_1, "myeloid_cluster_1_markers", w=8,h=6)


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


######
## Neutrophils have few genes expressed?
######
umap_nfeaturemyeloid<-FeaturePlot(d10x.combined_myeloid, features = "nFeature_RNA", min.cutoff = "q9", pt.size=1)
umap_nfeaturemyeloid + scale_color_continuous_sequential(palette = "ag_GrnYl") 

umap_ncountmyeloid<-FeaturePlot(d10x.combined_myeloid, features = "nCount_RNA", min.cutoff = "q9", pt.size=1)
umap_ncountmyeloid + scale_color_continuous_sequential(palette = "ag_GrnYl") 

#### Add back in cell cycle cause missing
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
d10x_raw_myeloid <- CellCycleScoring(d10x_raw_myeloid, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
identical(colnames(d10x_raw_myeloid), colnames(d10x.combined_myeloid))
d10x.combined_myeloid <- AddMetaData(d10x.combined_myeloid, metadata = d10x_raw_myeloid$Phase, col.name = 'Phase')


myeloid_cycle<-DimPlot(d10x.combined_myeloid, group.by = "Phase",pt.size=1) + scale_color_manual(values=c("#e6ab02","#386cb0","#1b9e77"))
myeloid_cycle
save_plts(myeloid_cycle, "myeloid_cell_cycle", w=8,h=6)



## RBC stuff
# Among the cognate receptors for these senescence signals that are expressed on RPMs and KCs, 
# TAM receptors AXL and MERTK, TIM4 and Fc receptor CD16 are the most abundant. A “don’t eat me” 
#signal that prevents the clearance of young and intact erythrocytes is provided by the interaction 
# between SIRPα at the surface of the macrophage with CD47 expressed by erythrocytes.
RBC_receptors<-FeaturePlot(d10x.combined_myeloid, features = c("TREM2","CD9","AXL","CD47","MERTK","TIMD4","FCGR3A"), min.cutoff = "q9", pt.size=1)
RBC_receptors
save_plts(RBC_receptors, "RBC_receptors", w=8,h=6)

FeaturePlot(d10x.combined_myeloid, features = c("HBA2", "AXL"), blend = TRUE)
RBC_Overlap<-FeaturePlot(d10x.combined_myeloid, features = c("HBA2", "FCGR3A"), blend = TRUE)
RBC_Overlap
save_plts(RBC_Overlap, "RBC_hem_receptor_Overlap", w=14,h=3)


d10x.combined_kupffer<-subset(d10x.combined_myeloid, subset = seurat_clusters %in% c(3,4,7))
rm(d10x.combined)
gc()
d10x.combined_kupffer <- RunPCA(d10x.combined_kupffer, npcs = 30, verbose = FALSE)
d10x.combined_kupffer <- RunUMAP(d10x.combined_kupffer, reduction = "pca", dims = 1:30)
d10x.combined_kupffer <- RunTSNE(d10x.combined_kupffer, reduction = "pca", dims = 1:30)
d10x.combined_kupffer <- FindNeighbors(d10x.combined_kupffer, reduction = "pca", dims = 1:30)
d10x.combined_kupffer <- FindClusters(d10x.combined_kupffer, resolution = 0.1)
DimPlot(d10x.combined_kupffer, label=T)
DimPlot(d10x.combined_kupffer, label=T, reduction="tsne")

RBC_receptors_KC<-FeaturePlot(d10x.combined_kupffer, features = c("TREM2","CD9","AXL","CD47","MERTK","TIMD4","FCGR3A"), min.cutoff = "q9", pt.size=1)
RBC_receptors_KC
RBC_hem_KC<-FeaturePlot(d10x.combined_kupffer, features = c("HBA2","HBA1","HBB","MARCO"), min.cutoff = "q9", pt.size=1)
RBC_hem_KC
RBC_Overlap<-FeaturePlot(d10x.combined_kupffer, features = c("HBA2", "FCGR3A"), blend = TRUE)
RBC_Overlap



######
## Kupffer-like cells
######
# Guilliams KC markers
KC_markers<-FeaturePlot(d10x.combined_myeloid, features = c("SLC16A9","CD5L","SLC40A1","VSIG4"), min.cutoff = "q9", pt.size=1)
KC_markers

# Guilliams Bile-duct associated lipid associated macropahges (BC-LAM) used human orthologos
LAM_markers<-FeaturePlot(d10x.combined_myeloid, features = c("TREM2","CLEC4A","GPNMB","SPP1","TREM2","LGALS1","FTL"), min.cutoff = "q9", pt.size=1)
LAM_markers



###########
## Label KC-like and recent recruited mono
###########
d10x_raw_myeloid@meta.data$CellType_rough[which(d10x_raw_myeloid@meta.data$seurat_clusters=="1")]<-"RR_myeloid"
d10x_raw_myeloid@meta.data$CellType_rough[which(d10x_raw_myeloid@meta.data$seurat_clusters%in%c("3","4"))]<-"KC_like"
d10x_raw_myeloid@meta.data$CellType_rough[which(d10x_raw_myeloid@meta.data$seurat_clusters=="7")]<-"macro_RBC"
d10x_raw_myeloid@meta.data$CellType_rough[which(d10x_raw_myeloid@meta.data$seurat_clusters=="5")]<-"bcell"
d10x_raw_myeloid@meta.data$CellType_rough[which(d10x_raw_myeloid@meta.data$seurat_clusters%in%c("2","6"))]<-"Neutrophil"
d10x_raw_myeloid@meta.data$CellType_rough[which(d10x_raw_myeloid@meta.data$seurat_clusters=="0")]<-"badmyeloid"


Idents(d10x_raw_myeloid) <- "CellType_rough"

de<-FindMarkers(d10x_raw_myeloid, ident.1 = "RR_myeloid", ident.2 = "KC_like", test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]


#########
## pathway KC versus recently recruited
#########
source("R_functions/GSEA_function_tmp.R")

GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

gene_list = de$avg_log2FC
names(gene_list) = rownames(de)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

FeaturePlot(d10x.combined_myeloid, features = c("VCAN","LYZ","CD5L","MARCO"), min.cutoff = "q9", pt.size=1)

res = GSEA(gene_list, GO_file, pval = 0.05)
top10_pathways<-data.frame(pathway=sapply(1:10, function(x) strsplit(res$Results$pathway[x], "%")[[1]][1]), direction=sapply(1:10, function(x) res$Results$Enrichment[x]))

res$Plot
View(top10_pathways)

plt_path<-res$Results
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Enrichment_Cell<-"Up-regulated in \nrecently recruited myeloid"
plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \nKupffer-like cells"

plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))

plt_path$direction_label<-as.factor(plt_path$Enrichment)
levels(plt_path$direction_label)<-c(0.1,-0.1)
plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))

GO_KC_vs_recentrcruti<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_vline(xintercept = 0, color="grey40")+scale_fill_manual(values=c("#fd8d3c","#6baed6"))+ 
  guides(fill = guide_legend(override.aes = list(size=5)))
save_plts(GO_KC_vs_recentrcruti, "GO_KC_vs_recently_recruited", w=20,h=10)


d10x.combined_myeloid@meta.data$CellType_rough[which(d10x.combined_myeloid@meta.data$seurat_clusters=="1")]<-"Recently recruited myeloid"
d10x.combined_myeloid@meta.data$CellType_rough[which(d10x.combined_myeloid@meta.data$seurat_clusters%in%c("3","4"))]<-"Kupffer-like"
d10x.combined_myeloid@meta.data$CellType_rough[which(d10x.combined_myeloid@meta.data$seurat_clusters=="7")]<-"Kupffer-like RBC"
d10x.combined_myeloid@meta.data$CellType_rough[which(d10x.combined_myeloid@meta.data$seurat_clusters=="5")]<-"B cell"
d10x.combined_myeloid@meta.data$CellType_rough[which(d10x.combined_myeloid@meta.data$seurat_clusters%in%c("2","6"))]<-"Neutrophil"
d10x.combined_myeloid@meta.data$CellType_rough[which(d10x.combined_myeloid@meta.data$seurat_clusters=="0")]<-"Low Quality Myeloid"

de_groups<-DimPlot(d10x.combined_myeloid, label=T, group.by="CellType_rough")+scale_color_manual(values=c("#c994c7","#fd8d3c","#8c2d04","grey","#1a9850","#6baed6"))
de_groups
save_plts(de_groups, "myeloid_DE_groups", w=10,h=6)

