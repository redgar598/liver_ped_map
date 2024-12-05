## Load Libraries
library(here)
library(Seurat)

#library(SCINA)
library(reshape2)
library(dplyr)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

library(sp)
#library(rgeos)
library(scales)
library(gridExtra)
#library(rhdf5)

source("scripts/00_pretty_plots.R")
source("scripts/00_long_functions.R")




#only healthy
#d10x<-readRDS(file = here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_d10x_adult_ped_raw.rds"))
d10x<-readRDS(file = here("../../../projects/macparland/RE/PediatricAdult/processed_data/","IFALD_d10x_adult_ped_raw.rds"))

#load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
load(here("../../../projects/macparland/RE/PediatricAdult/processed_data/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))


cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
gc()

d10x <- subset(d10x, subset = Treatment == "Healthy")
gc()

levels(d10x$CellType_refined)[which(levels(d10x$CellType_refined)=="Myeloid Erythrocytes\n(phagocytosis)" )]<-"Myeloid Erythrocytes"
levels(d10x$CellType_refined)[which(levels(d10x$CellType_refined)=="Macrophage\n(MHCII high)" )]<-"Macrophage MHCII high"

d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")
gc()



################
## load and merge all xenium BIDCell data
################
#
print("Load xenium BIDCell Data")

# count_files<-c("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_2/cell_gene_matrices/2024_03_21_16_22_10/expr_mat.csv",
#                "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_3/cell_gene_matrices/2024_03_21_16_38_48/expr_mat.csv",
#                "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_4/cell_gene_matrices/2024_03_21_17_03_05/expr_mat.csv",
#                "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C105/cell_gene_matrices/2024_05_17_13_32_03/expr_mat.csv",
#                "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C85/cell_gene_matrices/2024_05_17_13_29_40/expr_mat.csv",
#                "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C95/cell_gene_matrices/2024_07_04_12_35_29/expr_mat.csv",
#                "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C101/cell_gene_matrices/2024_06_28_09_23_40/expr_mat.csv")

count_files<-c("/cluster/projects/macparland/RE/xenium_liver/C94_2/cell_gene_matrices/2024_03_21_16_22_10/expr_mat.csv",
               "/cluster/projects/macparland/RE/xenium_liver/C94_3/cell_gene_matrices/2024_03_21_16_38_48/expr_mat.csv",
               "/cluster/projects/macparland/RE/xenium_liver/C94_4/cell_gene_matrices/2024_03_21_17_03_05/expr_mat.csv",
               "/cluster/projects/macparland/RE/xenium_liver/output-XETG00082__0016948__C105_C2__20231128__212539/cell_gene_matrices/2024_05_17_13_32_03/expr_mat.csv",
               "/cluster/projects/macparland/RE/xenium_liver/output-XETG00082__0016948__C85__20231128__212539/cell_gene_matrices/2024_05_17_13_29_40/expr_mat.csv",
               "/cluster/projects/macparland/RE/xenium_liver/output-XETG00082__0018061__C95_1__20240221__212447/cell_gene_matrices/2024_07_04_12_35_29/expr_mat.csv",
               "/cluster/projects/macparland/RE/xenium_liver/output-XETG00082__0016940__C101A__20231128__212539/cell_gene_matrices/2024_06_28_09_23_40/expr_mat.csv")



d10x.list <- sapply(count_files, function(file_path){
  counts<-read.csv(file_path)

  counts$X<-NULL
  rownames(counts) <- counts$cell_id
  counts$cell_id <- NULL

  # Transpose the data and convert to sparse matrix.
  mat <- as(t(as.matrix(counts)), "sparseMatrix")

  seu <- CreateSeuratObject(counts=mat)
  seu$sample<-strsplit(file_path,"/")[[1]][6]
  seu
})

d10x.list

xenium.obj <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "xenium_liver")
rm(d10x.list)
gc()

save(xenium.obj, file="/cluster/projects/macparland/RE/xenium_liver/BIDCell_liver_seurat.RData")
#save(xenium.obj, file="/media/redgar/Seagate Portable Drive/liver_BIDCell_output/BIDCell_liver_seurat.RData")


metadata_add<-xenium.obj@meta.data


## 477 genes in both
reference_subsample <- subset(d10x, features = intersect(rownames(d10x), rownames(xenium.obj)))
reference_subsample
rm(d10x)
gc()

## make seurat object for merging from xenium data
subset.matrix <- xenium.obj@assays$RNA
rm(xenium.obj)
gc()

xenium.obj_RNA <- CreateSeuratObject(subset.matrix) # Create a new Seurat object with just the genes of interest
identical(colnames(xenium.obj_RNA), rownames(metadata_add))
xenium.obj_RNA<-AddMetaData(xenium.obj_RNA, metadata_add)

xenium.obj_RNA<-JoinLayers(xenium.obj_RNA)


###############
## Integrate
###############
objetc.list<-list(xenium.obj_RNA, reference_subsample)

objetc.list <- lapply(X = objetc.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = objetc.list)
objetc.list <- lapply(X = objetc.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


## Identify anchors
chem.anchors <- FindIntegrationAnchors(object.list = objetc.list, anchor.features = features, reduction = "rpca")
ref_xenium_integrated <- IntegrateData(anchorset = chem.anchors)

DefaultAssay(ref_xenium_integrated) <- "integrated"



# Run the standard workflow for visualization and clustering
ref_xenium_integrated <- ScaleData(ref_xenium_integrated, verbose = FALSE)
ref_xenium_integrated <- RunPCA(ref_xenium_integrated, npcs = 30, verbose = FALSE)
ref_xenium_integrated <- RunUMAP(ref_xenium_integrated, reduction = "pca", dims = 1:30)
ref_xenium_integrated <- FindNeighbors(ref_xenium_integrated, reduction = "pca", dims = 1:30)
ref_xenium_integrated <- FindClusters(ref_xenium_integrated, resolution = 0.5)

ref_xenium_integrated$dataset<-"xenium"
ref_xenium_integrated$dataset[which(!(is.na(ref_xenium_integrated$CellType_refined)))]<-"scRNAseq_reference"


DimPlot(ref_xenium_integrated, group.by = "dataset", raster=FALSE)
DimPlot(ref_xenium_integrated, group.by = "CellType_refined",raster=FALSE)+colscale_cellType
DimPlot(ref_xenium_integrated, group.by = "seurat_clusters",raster=FALSE)


table(ref_xenium_integrated$CellType_refined,ref_xenium_integrated$seurat_clusters)


umap_mat<-as.data.frame(Embeddings(object = ref_xenium_integrated, reduction = "umap"))
umap_mat$cell<-rownames(umap_mat)
meta<-ref_xenium_integrated@meta.data
meta$cell<-rownames(meta)

plt_umap<-merge(meta, umap_mat, by="cell")

save(plt_umap, file="/cluster/projects/macparland/RE/xenium_liver/liver_reference_integration_labelling_BIDcell.RData")










####################
## Label by cluster
####################

load(here("data/liver_reference_integration_labelling_BIDcell.RData"))


len_x_bar<-((range(plt_umap$umap_1))[2]-(range(plt_umap$umap_1))[1])/10
len_y_bar<-((range(plt_umap$umap_2))[2]-(range(plt_umap$umap_2))[1])/10
arr <- list(x = min(plt_umap$umap_1), y = min(plt_umap$umap_2), x_len = len_x_bar, y_len = len_y_bar)


celltype_label<- plt_umap %>% group_by(CellType_refined) %>% summarise(mean_umap1=mean(umap_1),mean_umap2=mean(umap_2))
scRNAseq_reference_clusters<-ggplot()+
  geom_point(aes(umap_1, umap_2, fill=seurat_clusters),plt_umap, size=1.5, shape=21,stroke = NA)+
  geom_point(aes(umap_1, umap_2),plt_umap[which(plt_umap$dataset=="scRNAseq_reference"),], color="black", size=1, stroke = 1)+
  geom_point(aes(umap_1, umap_2, color=CellType_refined),plt_umap[which(plt_umap$dataset=="scRNAseq_reference"),], size=0.75)+
  colscale_cellType_ped+scale_fill_manual(values=sample(grep("^gr[ea]y", colours(), value = TRUE)[140:180],36))+
  geom_label(aes(mean_umap1, mean_umap2, label=CellType_refined), data=celltype_label,size=5, color="black", fill="white")+
  theme_void()+theme(legend.text=element_text(size=10),
                     legend.title=element_text(size=8),
                     plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=7,hjust = 0.05, vjust = 10),
                     axis.title.y = element_text(size=7,hjust = 0.05, vjust=-9,angle = 90))+
  annotate("segment",
           x = arr$x, xend = arr$x + c(arr$x_len, 0),
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt')))+
  guides(colour = guide_legend(override.aes = list(size=3)))
scRNAseq_reference_clusters

ggsave(scRNAseq_reference_clusters, file=here("figures/jpeg/sc_reference_clusters_BIDCell.jpeg"), h=20, w=25, bg = "white")




### Label based on scRNAseq clustering
cluster_label<- plt_umap %>% group_by(seurat_clusters) %>% summarise(mean_umap1=mean(umap_1),mean_umap2=mean(umap_2))

seurat_clusters<-ggplot()+
  geom_point(aes(umap_1, umap_2, color=seurat_clusters),plt_umap, size=1.5)+
  geom_point(aes(umap_1, umap_2),plt_umap[which(plt_umap$dataset=="scRNAseq_reference"),], color="black", size=1, stroke = 1)+
  geom_point(aes(umap_1, umap_2, color=seurat_clusters),plt_umap[which(plt_umap$dataset=="scRNAseq_reference"),], size=0.75)+
  geom_label(aes(mean_umap1, mean_umap2, label=seurat_clusters), data=cluster_label,size=5, color="black", fill="white")+
  theme_void()+theme(legend.text=element_text(size=10),
                     legend.title=element_text(size=8),
                     plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=7,hjust = 0.05, vjust = 10),
                     axis.title.y = element_text(size=7,hjust = 0.05, vjust=-9,angle = 90))+
  annotate("segment",
           x = arr$x, xend = arr$x + c(arr$x_len, 0),
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt')))+
  guides(colour = guide_legend(override.aes = list(size=3)))
seurat_clusters
ggsave(seurat_clusters, file=here("figures/jpeg/seuratclusters_BIDCell.jpeg"), h=20, w=25, bg = "white")



plt_bar_scRef<-plt_umap[c("seurat_clusters","CellType_refined")] %>% count(seurat_clusters, CellType_refined)
plt_bar_all<-plt_umap[c("seurat_clusters","CellType_refined")] %>% count(seurat_clusters)

top_vote<-do.call(rbind,lapply(1:nrow(plt_bar_all), function(x){
  cluster_count<-plt_bar_scRef[which(plt_bar_scRef$seurat_clusters==plt_bar_all$seurat_clusters[x]),]
  if(length(which(!(is.na(cluster_count$CellType_refined))))==0){
    data.frame(cluster=plt_bar_all$seurat_clusters[x],max_cell_allen=NA,percent_allen=NA)
  }else{
  cluster_count<-cluster_count[which(!(is.na(cluster_count$CellType_refined))),]
  max_allen<-max((cluster_count$n/plt_bar_all$n[x])*100)
  data.frame(cluster=plt_bar_all$seurat_clusters[x],
             max_cell_allen=cluster_count$CellType_refined[which(((cluster_count$n/plt_bar_all$n[x])*100)==max_allen)],
             percent_allen=max_allen)}}))
top_vote


barplots<-ggplot(plt_bar_scRef, aes(seurat_clusters, n, fill=CellType_refined))+geom_bar(stat="identity")+fillscale_cellType_ped
barplots
ggsave(barplots, file=here("figures/jpeg/seuratclusters_barplot_BIDCell.jpeg"), h=10, w=15, bg = "white")


ggplot()+
  geom_point(aes(umap_1, umap_2, color=seurat_clusters),plt_umap, size=0.1)+
  scale_color_manual(values=c(rep("grey",10),"red",rep("grey",30)))


#scRNAseq confident
c(5,12,16,17, 18, 19,23)

#Hep subtype
c(1:4, 14,22)

#HSC subtype
c(8,9)

#LSEC subtype
c(7,11)

#Colangiocyte subtype
c(15, 21)

#Refine
c(6,10,13,20)

####################
##check markers in xenium data
#####################
cell_markers_general<-read.csv(here("data/Xenium_CombinedPanel.csv"))
custom_markers<-read.csv(here("data/Xenium_MacParlandGeneListUpdated.csv"))
custom_markers$Gene[grep("HAL", custom_markers$Gene)]<-"HAL"


load(here("/cluster/projects/macparland/RE/xenium_liver/BIDCell_liver_seurat.RData"))
#load(here("data/BIDCell_liver_seurat.RData"))



xenium.obj <- NormalizeData(xenium.obj)
xenium.obj <- FindVariableFeatures(xenium.obj, selection.method = "vst", nfeatures = 2000)
xenium.obj <- ScaleData(xenium.obj)

xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.4)

BIDCell_UMAP<-DimPlot(xenium.obj,raster=FALSE, label=T)
save_plts(BIDCell_UMAP, "BIDCell_UMAP", w=10, h=10)


xenium.obj$sample<-sapply(1:ncol(xenium.obj), function(x) strsplit(colnames(xenium.obj)[x],"_")[[1]][2])
xenium.obj$sample<-as.factor(xenium.obj$sample)
levels(xenium.obj$sample)<-c("C94_2", "C94_3","C94_4","C105","C85","C95","C101")

BIDCell_UMAP<-DimPlot(xenium.obj,raster=FALSE, label=F, group.by = "sample")
save_plts(BIDCell_UMAP, "BIDCell_UMAP_sample", w=10, h=10)

FeaturePlot(xenium.obj, features=c("CYP2A6", "CYP1A2", "CYP2A7","CYP2E1"), raster=F)


#########

plt_umap_xenium<-plt_umap[which(plt_umap$dataset=="xenium"),]
plt_umap_xenium<-plt_umap_xenium[,-grep("integration.|_RNA|_color|_label|_id|_order|_name",colnames(plt_umap_xenium))]
plt_umap_xenium$reference_cluster<-plt_umap_xenium$seurat_clusters

plt_umap_xenium<-plt_umap_xenium[!duplicated(plt_umap_xenium),]
rownames(plt_umap_xenium)<-plt_umap_xenium$cell

plt_umap_xenium<-plt_umap_xenium[match(colnames(xenium.obj), rownames(plt_umap_xenium)),]
identical(colnames(xenium.obj), rownames(plt_umap_xenium))
plt_umap_xenium$cell<-NULL

xenium.obj<-AddMetaData(xenium.obj, plt_umap_xenium)

BIDCell_UMAP<-DimPlot(xenium.obj,raster=FALSE, group.by = "reference_cluster", label=T)
save_plts(BIDCell_UMAP, "BIDCell_UMAP_reference_clusters", w=10, h=10)

Idents(xenium.obj)<-"reference_cluster"


pdf(file = here("figures/custom_marker_dot_plot_BIDCell_integratedclusters.pdf"), w=15, h=10)
lapply(unique(custom_markers$Cell.Type), function(type){
  DotPlot(object = xenium.obj, features = custom_markers$Gene[which(custom_markers$Cell.Type==type)])+xlab(type)})
dev.off()

pdf(file = here("figures/custom_marker_dot_plot_BIDCell_integratedclusters_other_label.pdf"), w=15, h=10)
lapply(unique(custom_markers$X), function(type){
  DotPlot(object = xenium.obj, features = custom_markers$Gene[which(custom_markers$X==type)])+xlab(type)})
dev.off()


rownames(xenium.obj)[grep("HLA",rownames(xenium.obj))]

# MHCII high
FeaturePlot(xenium.obj, features=c("HLA.DQA1", "HLA.DQB1", "HLA.DQB2"), raster=F) #MHCII
DotPlot(object = xenium.obj, features = c("HLA.DQA1", "HLA.DQB1", "HLA.DQB2"))

# portal vs central
custom_markers[grep("CYP",custom_markers$Gene),]
DotPlot(object = xenium.obj, features = c("CYP1A2", "CYP2A6"))
FeaturePlot(xenium.obj, features=c("CYP2A6", "CYP1A2", "CYP2A7","CYP2E1"), raster=F)

## cycling
DotPlot(object = xenium.obj, features = c("TOP2A", "MKI67"))


FeaturePlot(xenium.obj, features=c("nCount_RNA"), raster=F) #QC
FeaturePlot(xenium.obj, features=c("nFeature_RNA"), raster=F) #QC

## cdc
custom_markers[grep("CLEC",custom_markers$Gene),]
DotPlot(object = xenium.obj, features = c("CLEC10A", "CLEC9A"))



### two types of LSEC, why different
Idents(xenium.obj)<-"reference_cluster"
xenium.obj<-JoinLayers(xenium.obj)
de_LSEC<-FindMarkers(xenium.obj,ident.1 = 7, ident.2 = 11,  min.pct = 0.25)

head(de_LSEC)
tail(de_LSEC)

cell_markers_general[which(cell_markers_general$Gene%in%rownames(tail(de_LSEC))),]
custom_markers[which(custom_markers$Gene%in%rownames(de_LSEC)[(nrow(de_LSEC)-10):nrow(de_LSEC)]),]

cell_markers_general[which(cell_markers_general$Gene%in%rownames(head(de_LSEC))),]
custom_markers[which(custom_markers$Gene%in%rownames(de_LSEC)[1:10]),]






###########################
## Label cell types
###########################

#scRNAseq confident
c(5,12,16,17, 18, 19,23)

#Hep subtype
c(1:4, 14,22)

#HSC subtype
c(8,9)

#LSEC subtype
c(7,11)

#Colangiocyte subtype
c(15, 21)

#Refine
c(6,10,13,20)

###############
## relabel based on reference
###############
plt_umap_xenium<-plt_umap[which(plt_umap$dataset=="xenium"),]
plt_umap_xenium<-plt_umap_xenium[,-grep("integration.|_RNA|_color|_label|_id|_order|_name",colnames(plt_umap_xenium))]
plt_umap_xenium$reference_cluster<-plt_umap_xenium$seurat_clusters

plt_umap_xenium<-merge(plt_umap_xenium, top_vote, by.x="reference_cluster", by.y="cluster")

plt_umap_xenium$CellType<-as.character(plt_umap_xenium$reference_cluster)
plt_umap_xenium$CellType[which(plt_umap_xenium$reference_cluster%in%c(5,12,15,16,17, 18, 19,23))]<-plt_umap_xenium$max_cell_allen[which(plt_umap_xenium$reference_cluster%in%c(5,12,15,16,17, 18, 19,23))]
plt_umap_xenium$CellType[which(plt_umap_xenium$reference_cluster%in%c(21))]<-"Cholangiocyte (Biliary)"

plt_umap_xenium$CellType[which(plt_umap_xenium$reference_cluster%in%c(0:4, 14,22))]<-"Hepatocyte"

plt_umap_xenium$CellType[which(plt_umap_xenium$reference_cluster%in%c(7))]<-"LSEC"
plt_umap_xenium$CellType[which(plt_umap_xenium$reference_cluster%in%c(11))]<-"VEC"

plt_umap_xenium$CellType[which(plt_umap_xenium$reference_cluster%in%c(9))]<-"HSC (Quiescent)"
plt_umap_xenium$CellType[which(plt_umap_xenium$reference_cluster%in%c(8))]<-"HSC (Activated)"

plt_umap_xenium$CellType[which(plt_umap_xenium$reference_cluster%in%c(6))]<-"CD3+ T-cells"
plt_umap_xenium$CellType[which(plt_umap_xenium$reference_cluster%in%c(20))]<-"B cells"
plt_umap_xenium$CellType[which(plt_umap_xenium$reference_cluster%in%c(10))]<-"Myeloid"
plt_umap_xenium$CellType[which(plt_umap_xenium$reference_cluster%in%c(13))]<-"Hepatocyte (Cycling)"

# 
# 
# 
# 
# 
# 
# 
# 
# 
# ggplot(plt_umap_xenium, aes(umap_1,umap_2))+
#   geom_point(size = 0.06, colour= "black", stroke = 1)+
#   geom_point(aes(color=CellType),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
#   colscale_cellType+
#   annotate("segment",
#            x = arr$x, xend = arr$x + c(arr$x_len, 0),
#            y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
#            arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
#   theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
#                      axis.title.x = element_text(size=5,hjust = 0.05),
#                      axis.title.y = element_text(size=5,hjust = 0.05,angle = 90))+
#   guides(colour = guide_legend(override.aes = list(size=3)))
# 
# 
# 
# plt_umap_xenium<-plt_umap_xenium[,c("cell","sample","CellType")]
# plt_umap_xenium<-plt_umap_xenium[!duplicated(plt_umap_xenium$cell),]
# 
# plt_umap_xenium$sample<-sapply(1:nrow(plt_umap_xenium), function(x) strsplit(plt_umap_xenium$cell[x],"_")[[1]][2])
# plt_umap_xenium$sample<-as.factor(plt_umap_xenium$sample)
# levels(plt_umap_xenium$sample)<-c("C94_2", "C94_3","C94_4","C105","C85","C95","C101")
# 
# 
# save(plt_umap_xenium, file=here("data/cell_type_labels_BIDCell_rough.RData"))


######################
#### Refine Hepatocytes
######################
#custom_markers[grep("CYP",custom_markers$Gene),]


load(here("data/cell_type_labels_BIDCell_rough.RData"))

#load(here("data/BIDCell_liver_seurat.RData"))
load(here("/cluster/projects/macparland/RE/xenium_liver/BIDCell_liver_seurat.RData"))


xenium.obj <- NormalizeData(xenium.obj)
xenium.obj <- FindVariableFeatures(xenium.obj, selection.method = "vst", nfeatures = 2000)
xenium.obj <- ScaleData(xenium.obj)

xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.4)

BIDCell_UMAP<-DimPlot(xenium.obj,raster=FALSE, label=T)
save_plts(BIDCell_UMAP, "BIDCell_UMAP", w=10, h=10)


xenium.obj$sample<-sapply(1:ncol(xenium.obj), function(x) strsplit(colnames(xenium.obj)[x],"_")[[1]][2])
xenium.obj$sample<-as.factor(xenium.obj$sample)
levels(xenium.obj$sample)<-c("C94_2", "C94_3","C94_4","C105","C85", "C95", "C101")

BIDCell_UMAP<-DimPlot(xenium.obj,raster=FALSE, label=F, group.by = "sample")
save_plts(BIDCell_UMAP, "BIDCell_UMAP_sample", w=10, h=10)

FeaturePlot(xenium.obj, features=c("CYP2A6", "CYP1A2", "CYP2A7","CYP2E1"), raster=F)



plt_umap_xenium<-plt_umap_xenium[match(colnames(xenium.obj), plt_umap_xenium$cell),]
identical(colnames(xenium.obj),plt_umap_xenium$cell)
rownames(plt_umap_xenium)<-plt_umap_xenium$cell
plt_umap_xenium$cell<-NULL

xenium.obj<-AddMetaData(xenium.obj, plt_umap_xenium)



xenium.obj_hep<-subset(xenium.obj, subset = CellType %in% c("Hepatocyte"))

xenium.obj_hep_C94<-subset(xenium.obj_hep, subset = sample %in% c("C94_2","C94_3","C94_4"))
xenium.obj_hep_C85<-subset(xenium.obj_hep, subset = sample %in% c("C85"))
xenium.obj_hep_C105<-subset(xenium.obj_hep, subset = sample %in% c("C105"))
xenium.obj_hep_C95<-subset(xenium.obj_hep, subset = sample %in% c("C95"))
xenium.obj_hep_C101<-subset(xenium.obj_hep, subset = sample %in% c("C101"))



cluster_hep<-function(d10x){
  d10x <- NormalizeData(d10x)
  d10x <- FindVariableFeatures(d10x, selection.method = "vst", nfeatures = 2000)
  d10x <- ScaleData(d10x)

  d10x <- RunPCA(d10x, npcs = 30, features = rownames(d10x))
  d10x <- RunUMAP(d10x, dims = 1:30)
  d10x <- FindNeighbors(d10x, reduction = "pca", dims = 1:30)
  d10x <- FindClusters(d10x, resolution = 0.6)
  d10x
}

xenium.obj_hep_C94<-cluster_hep(xenium.obj_hep_C94)
xenium.obj_hep_C85<-cluster_hep(xenium.obj_hep_C85)
xenium.obj_hep_C105<-cluster_hep(xenium.obj_hep_C105)
xenium.obj_hep_C95<-cluster_hep(xenium.obj_hep_C95)
xenium.obj_hep_C101<-cluster_hep(xenium.obj_hep_C101)

save(xenium.obj_hep_C94, xenium.obj_hep_C85, xenium.obj_hep_C105, xenium.obj_hep_C95, xenium.obj_hep_C101, file=here("/cluster/projects/macparland/RE/xenium_liver","hepatocyte_each_sample.RData"))

load(here("/media/redgar/Seagate Portable Drive/xenium_liver/","hepatocyte_each_sample.RData"))


plot_grid(plot_grid(
  DimPlot(xenium.obj_hep_C94,raster=FALSE, label=T),
  DotPlot(xenium.obj_hep_C94, features=c("CYP2A6", "CYP1A2","CYP2E1", "CYP2A7", "MKI67","TOP2A"))),
  FeaturePlot(xenium.obj_hep_C94, features=c("CYP2A7", "CYP1A2", "MKI67","TOP2A"), raster=F), ncol=1, rel_heights = c(1,2) )


Idents(xenium.obj_hep_C94)<-"seurat_clusters"
xenium.obj_hep_C94<-JoinLayers(xenium.obj_hep_C94)

C94_labels<-xenium.obj_hep_C94@meta.data
C94_labels$CellType[which(C94_labels$seurat_clusters%in%c(0,1,8,5))]<-"Hepatocyte (Pericentral)"
C94_labels$CellType[which(C94_labels$seurat_clusters%in%c(2,3,4,6,7))]<-"Hepatocyte (Periportal)"




plot_grid(plot_grid(
  DimPlot(xenium.obj_hep_C85,raster=FALSE, label=T),
  DotPlot(xenium.obj_hep_C85, features=c("CYP2A6", "CYP1A2","CYP2E1", "CYP2A7", "MKI67","TOP2A"))),
  FeaturePlot(xenium.obj_hep_C85, features=c("CYP2A7", "CYP1A2", "MKI67","TOP2A"), raster=F), ncol=1, rel_heights = c(1,2) )

plot_grid(DimPlot(xenium.obj_hep_C85,raster=FALSE, label=T),
          FeaturePlot(xenium.obj_hep_C85, features=c("CYP3A4", "CYP2E1", "CYP1A2","CYP2A7","CYP2A6"), raster=F))


C85_labels<-xenium.obj_hep_C85@meta.data
C85_labels$CellType[which(C85_labels$seurat_clusters%in%c(1,2,8,5,9,4,7,10:13))]<-"Hepatocyte (Pericentral)"
C85_labels$CellType[which(C85_labels$seurat_clusters%in%c(0,3,6))]<-"Hepatocyte (Periportal)"



plot_grid(plot_grid(
  DimPlot(xenium.obj_hep_C105,raster=FALSE, label=T),
  DotPlot(xenium.obj_hep_C105, features=c("CYP2A6","CYP2A7", "CYP1A2", "CYP3A4","CYP2E1", "MKI67","TOP2A"))),
  FeaturePlot(xenium.obj_hep_C105, features=c("CYP2A7", "CYP1A2", "MKI67","TOP2A"), raster=F), ncol=1, rel_heights = c(1,2) )

FeaturePlot(xenium.obj_hep_C105, features=c("CYP3A4", "CYP2E1", "CYP1A2","CYP2A7","CYP2A6"), raster=F)

C105_labels<-xenium.obj_hep_C105@meta.data
C105_labels$CellType[which(C105_labels$seurat_clusters%in%c(1,2,3,4,5,6,8,9))]<-"Hepatocyte (Pericentral)"
C105_labels$CellType[which(C105_labels$seurat_clusters%in%c(0,7))]<-"Hepatocyte (Periportal)"


plot_grid(plot_grid(
  DimPlot(xenium.obj_hep_C95,raster=FALSE, label=T),
  DotPlot(xenium.obj_hep_C95, features=c("CYP2A6","CYP2A7", "CYP1A2", "CYP3A4","CYP2E1", "MKI67","TOP2A"))),
  FeaturePlot(xenium.obj_hep_C95, features=c("CYP2A7", "CYP1A2", "MKI67","TOP2A"), raster=F), ncol=1, rel_heights = c(1,2) )

FeaturePlot(xenium.obj_hep_C95, features=c("CYP3A4", "CYP2E1", "CYP1A2","CYP2A7","CYP2A6"), raster=F)

C95_labels<-xenium.obj_hep_C95@meta.data
C95_labels$CellType[which(C95_labels$seurat_clusters%in%c(1,2,3,6,10,11,15))]<-"Hepatocyte (Pericentral)"
C95_labels$CellType[which(C95_labels$seurat_clusters%in%c(0,4,5,7,8,9,12,13,14,16))]<-"Hepatocyte (Periportal)"


plot_grid(plot_grid(
  DimPlot(xenium.obj_hep_C101,raster=FALSE, label=T),
  DotPlot(xenium.obj_hep_C101, features=c("CYP2A6","CYP2A7", "CYP1A2", "CYP3A4","CYP2E1", "MKI67","TOP2A"))),
  FeaturePlot(xenium.obj_hep_C101, features=c("CYP2A7", "CYP1A2", "MKI67","TOP2A"), raster=F), ncol=1, rel_heights = c(1,2) )

FeaturePlot(xenium.obj_hep_C101, features=c("CYP3A4", "CYP2E1", "CYP1A2","CYP2A7","CYP2A6"), raster=F)

C101_labels<-xenium.obj_hep_C101@meta.data
C101_labels$CellType[which(C101_labels$seurat_clusters%in%c(0,2,6,4,5,7,8))]<-"Hepatocyte (Pericentral)"
C101_labels$CellType[which(C101_labels$seurat_clusters%in%c(1,3,9:16))]<-"Hepatocyte (Periportal)"

rm(xenium.obj_hep_C101,xenium.obj_hep_C105, xenium.obj_hep_C85, xenium.obj_hep_C94, xenium.obj_hep_C95)
gc()


######################
#### Refine Myeloid
######################
load(here("data/cell_type_labels_BIDCell_rough.RData"))
load(here("data/BIDCell_liver_seurat.RData"))



plt_umap_xenium<-plt_umap_xenium[match(colnames(xenium.obj), plt_umap_xenium$cell),]
identical(colnames(xenium.obj),plt_umap_xenium$cell)
rownames(plt_umap_xenium)<-plt_umap_xenium$cell
plt_umap_xenium$cell<-NULL

xenium.obj<-AddMetaData(xenium.obj, plt_umap_xenium)

xenium.obj$sample<-sapply(1:ncol(xenium.obj), function(x) strsplit(colnames(xenium.obj)[x],"_")[[1]][2])
xenium.obj$sample<-as.factor(xenium.obj$sample)
levels(xenium.obj$sample)<-c("C94_2", "C94_3","C94_4","C105","C85","C95","C101")

xenium.obj_myeloid<-subset(xenium.obj, subset = CellType %in% c("Myeloid","Cycling Myeloid","KC Like"))

xenium.obj_myeloid <- NormalizeData(xenium.obj_myeloid)
xenium.obj_myeloid <- FindVariableFeatures(xenium.obj_myeloid, selection.method = "vst", nfeatures = 2000)
xenium.obj_myeloid <- ScaleData(xenium.obj_myeloid)
xenium.obj_myeloid <- RunPCA(xenium.obj_myeloid, npcs = 30, features = rownames(xenium.obj_myeloid))
xenium.obj_myeloid <- RunUMAP(xenium.obj_myeloid, dims = 1:30)
xenium.obj_myeloid <- FindNeighbors(xenium.obj_myeloid, reduction = "pca", dims = 1:30)
xenium.obj_myeloid <- FindClusters(xenium.obj_myeloid, resolution = 0.4)

xenium.obj_myeloid$sample<-sapply(1:ncol(xenium.obj_myeloid), function(x) strsplit(colnames(xenium.obj_myeloid)[x],"_")[[1]][2])
xenium.obj_myeloid$sample<-as.factor(xenium.obj_myeloid$sample)
levels(xenium.obj_myeloid$sample)<-c("C94_2", "C94_3","C94_4","C105","C85","C95","C101")

save(xenium.obj_myeloid, file=here("/media/redgar/Seagate Portable Drive/xenium_liver/","myeloid_cells.RData"))


plot_grid(plot_grid(DimPlot(xenium.obj_myeloid,raster=FALSE, label=T),
          DimPlot(xenium.obj_myeloid,raster=FALSE, label=F, group.by = "sample"), ncol=1),
          DotPlot(xenium.obj_myeloid, features=c("MARCO", "CD68","CD5L", "LYZ","S100A9","S100A8","HLA.DQA1", "HLA.DQB1", "HLA.DQB2","TOP2A", "MKI67")), ncol=2)

DimPlot(xenium.obj_myeloid,raster=FALSE, label=F, group.by = "CellType")


plot_grid(plot_grid(DimPlot(xenium.obj_myeloid,raster=FALSE, label=T),
                    DimPlot(xenium.obj_myeloid,raster=FALSE, label=F, group.by = "sample"), ncol=2),
          FeaturePlot(xenium.obj_myeloid, features=c("MARCO", "CD68","CD5L", "LYZ","S100A9","S100A8","HLA.DQA1", "HLA.DQB1", "HLA.DQB2","TOP2A", "MKI67"), raster=F), ncol=1, rel_heights = c(1,3))

plot_grid(plot_grid(DimPlot(xenium.obj_myeloid,raster=FALSE, label=T),
                    DimPlot(xenium.obj_myeloid,raster=FALSE, label=F, group.by = "sample"), ncol=2),
          FeaturePlot(xenium.obj_myeloid, features=c("CLEC10A", "CLEC9A"), raster=F), ncol=1)

plot_grid(plot_grid(DimPlot(xenium.obj_myeloid,raster=FALSE, label=T),
                    DimPlot(xenium.obj_myeloid,raster=FALSE, label=F, group.by = "sample"), ncol=2),
          FeaturePlot(xenium.obj_myeloid, features=c("PECAM1", "CD36","PTPRC"), raster=F), ncol=1)

plot_grid(plot_grid(DimPlot(xenium.obj_myeloid,raster=FALSE, label=T),
                    DimPlot(xenium.obj_myeloid,raster=FALSE, label=F, group.by = "sample"), ncol=2),
          FeaturePlot(xenium.obj_myeloid, features=c("HLA.DQA1", "HLA.DQB1", "HLA.DQB2"), raster=F), ncol=1, rel_heights = c(1,2))


myeloid_labels<-xenium.obj_myeloid@meta.data
myeloid_labels$CellType[which(myeloid_labels$seurat_clusters%in%c(4))]<-"Macrophage MHCII High"
myeloid_labels$CellType[which(myeloid_labels$seurat_clusters%in%c(2))]<-"Mono-Mac"
myeloid_labels$CellType[which(myeloid_labels$seurat_clusters%in%c(1,0,3,5,6,7))]<-"KC Like"
myeloid_labels$CellType[which(myeloid_labels$seurat_clusters%in%c(8))]<-"Cycling Myeloid"

######################
#### Refine B cells
######################
load(here("data/cell_type_labels_BIDCell_rough.RData"))
load(here("data/BIDCell_liver_seurat.RData"))

plt_umap_xenium<-plt_umap_xenium[match(colnames(xenium.obj), plt_umap_xenium$cell),]
identical(colnames(xenium.obj),plt_umap_xenium$cell)
rownames(plt_umap_xenium)<-plt_umap_xenium$cell
plt_umap_xenium$cell<-NULL

xenium.obj<-AddMetaData(xenium.obj, plt_umap_xenium)

xenium.obj$sample<-sapply(1:ncol(xenium.obj), function(x) strsplit(colnames(xenium.obj)[x],"_")[[1]][2])
xenium.obj$sample<-as.factor(xenium.obj$sample)
levels(xenium.obj$sample)<-c("C94_2", "C94_3","C94_4","C105","C85","C95","C101")

xenium.obj_bcells<-subset(xenium.obj, subset = reference_cluster %in% c(19,20))

xenium.obj_bcells <- NormalizeData(xenium.obj_bcells)
xenium.obj_bcells <- FindVariableFeatures(xenium.obj_bcells, selection.method = "vst", nfeatures = 2000)
xenium.obj_bcells <- ScaleData(xenium.obj_bcells)
xenium.obj_bcells <- RunPCA(xenium.obj_bcells, npcs = 30, features = rownames(xenium.obj_bcells))
xenium.obj_bcells <- RunUMAP(xenium.obj_bcells, dims = 1:30)
xenium.obj_bcells <- FindNeighbors(xenium.obj_bcells, reduction = "pca", dims = 1:30)
xenium.obj_bcells <- FindClusters(xenium.obj_bcells, resolution = 0.4)


plot_grid(plot_grid(DimPlot(xenium.obj_bcells,raster=FALSE, label=T),
                    DimPlot(xenium.obj_bcells,raster=FALSE, label=F, group.by = "sample"), ncol=1),
          DotPlot(xenium.obj_bcells, features=c("CD19", "IGHG1","MS4A1","TOP2A", "MKI67")), ncol=2)
FeaturePlot(xenium.obj_bcells, features=c("CD19", "IGHG1","MS4A1","TOP2A", "MKI67"))

b_labels<-xenium.obj_bcells@meta.data
b_labels$CellType[which(b_labels$seurat_clusters%in%c(0,2,5))]<-"Plasma Cells"
b_labels$CellType[which(b_labels$seurat_clusters%in%c(1,3,4))]<-"Mature B-cells"



rm(xenium.obj)
gc()

######################
#### Add Refined Hepatocyte and Myeloid labels
######################
plt_umap_xenium$CellType<-sapply(1:nrow(plt_umap_xenium), function(x){
  if(rownames(plt_umap_xenium)[x]%in%rownames(C94_labels)){
    C94_labels$CellType[which(rownames(C94_labels)==rownames(plt_umap_xenium)[x])]
  }else{
    if(rownames(plt_umap_xenium)[x]%in%rownames(C85_labels)){
      C85_labels$CellType[which(rownames(C85_labels)==rownames(plt_umap_xenium)[x])]
    }else{
        if(rownames(plt_umap_xenium)[x]%in%rownames(C105_labels)){
          C105_labels$CellType[which(rownames(C105_labels)==rownames(plt_umap_xenium)[x])]
        }else{
          if(rownames(plt_umap_xenium)[x]%in%rownames(C101_labels)){
            C101_labels$CellType[which(rownames(C101_labels)==rownames(plt_umap_xenium)[x])]
            }else{
              if(rownames(plt_umap_xenium)[x]%in%rownames(C95_labels)){
                C95_labels$CellType[which(rownames(C95_labels)==rownames(plt_umap_xenium)[x])]
                }else{
                  if(rownames(plt_umap_xenium)[x]%in%rownames(myeloid_labels)){
                    myeloid_labels$CellType[which(rownames(myeloid_labels)==rownames(plt_umap_xenium)[x])]
                    }else{
                      if(rownames(plt_umap_xenium)[x]%in%rownames(b_labels)){
                        b_labels$CellType[which(rownames(b_labels)==rownames(plt_umap_xenium)[x])]
                        }else{plt_umap_xenium$CellType[x]}
          }}}}}}})

plt_umap_xenium$cell<-rownames(plt_umap_xenium)

# ## from spatial exploring
# plt_umap_xenium$CellType[which(plt_umap_xenium$CellType=="LSEC (Periportal)")]<-"VEC"

plt_umap_xenium$sample<-sapply(1:nrow(plt_umap_xenium), function(x) strsplit(plt_umap_xenium$cell[x],"_")[[1]][2])
plt_umap_xenium$sample<-as.factor(plt_umap_xenium$sample)
levels(plt_umap_xenium$sample)<-c("C94_2", "C94_3","C94_4","C105","C85","C95","C101")

save(plt_umap_xenium, file=here("data/cell_type_labels_BIDCell.RData"))
