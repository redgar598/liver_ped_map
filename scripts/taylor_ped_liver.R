### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(reshape2)
library(gtools)
library(SoupX)
library(colorspace)
library(cowplot)
library(DropletQC)
library(SCINA)



source("scripts/00_pretty_plots.R")
source("scripts/00_entropy_d10x.R")
source("scripts/00_fanciest_UMAP.R")
source("scripts/00_plot_gene_exp.R")


taylor_liver_h5ad <- Read10X_h5(here("data/GSM4983457_NC_filtered_feature_bc_matrix.h5"))

taylor_liver <- CreateSeuratObject(
  counts = taylor_liver_h5ad,
  min.cells = 0, min.features = 0
)

head(taylor_liver)


######################
## cell cycle gene expression
######################
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

taylor_liver <- CellCycleScoring(taylor_liver, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

taylor_liver[["percent.mt"]] <- PercentageFeatureSet(taylor_liver, pattern = "^MT-")
taylor_liver <- subset(taylor_liver, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 25)


################
## Normalize scale and UMAP
################
taylor_liver <- SCTransform(taylor_liver, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), verbose = FALSE)

# dimension reduction
taylor_liver <- RunPCA(taylor_liver, verbose = FALSE)
taylor_liver <- RunUMAP(taylor_liver, dims = 1:30)
taylor_liver <- RunTSNE(taylor_liver, dims = 1:30)

# cluster
taylor_liver <- FindNeighbors(taylor_liver, reduction = "pca", dims = 1:20)
taylor_liver <- FindClusters(taylor_liver, resolution = 0.5)


DimPlot(taylor_liver, reduction = "umap", pt.size=0.25, label=T)
DimPlot(taylor_liver, reduction = "tsne", pt.size=0.25, label=T)


###################
## rough cell labels
###################
Macrophage_genes<-c( "PTPRC", "CD68", "MARCO","CD5L","VSIG4", "MAF", "LYZ", "CSTA", "S100A8", "S10049",
                     "CD14", "CD74", "GPBAR1", "ID3")

B_genes<-c("POU2F2","FCER2","MS4A1","LTB","CD37","CD79B","IGLC2","IGHG1","IGKC", "CD19")
NK_T_genes<-c("CD3D","TRDC","SELL","IL7R","CCR7","S100A4","CD8A","GNLY","NKG7")


LEC_genes<-c("CALCRL", "VWF", "RAMP2", "STAB2", "LYVE1", "PECAM1", "ENG", "FCGR2B", "F8", "SPARCL1",
             "ID1", "SOX18", "CD32B", "ID3")
Hepatocyte_genes<-c("ALB", "HAMP", "ARG1", "PCK1", "AFP", "BCHE", "HAL", "SCD", "CPS1", "CYP3A4",
                    "ELF3", "CRP", "GSTA2", "AKR1C1", "MGST1", "CYP3A5", "ALDH1A1", "ADH1A", "CYP2E1",
                    "GLS2", "SDS", "GLUL", "AKR1D1", "HPR",
                    "HMGCS1", "IGSF23", "ACSS2", "G6PC", "ID3")
Cholangiocytes_genes<-c( "EPCAM", "SOX9", "KRT1", "KRT7", "ANXA4", "KRT18", "ID3")

HSCs_genes<-c( "RBP1", "LRAT", "PDE3B", "ACTA2", "AOX1", "PDE3D", "PDE4D", "SPARC", "TAGLN", "COL1A1", "COL1A2", "COL3A1",
               "TIMP1", "DCN", "MYL9", "TPM2", "MEG3", "BGN", "IGFBP7", "IGFBP3", "CYR61", "IGFBP6", "CCL2", "COLEC11",
               "CTGF", "HGF", "ID3")

T_genes<-c("CD3D","IL7R","CD8A","IL32")
NK_genes<-c("NKG7","CD7")
gd_genes<-c("GNLY")

genes<-unique(c(Macrophage_genes, NK_T_genes, B_genes,LEC_genes,Hepatocyte_genes,Cholangiocytes_genes,HSCs_genes))

taylor_liver.exp<-as.data.frame(taylor_liver[["RNA"]]@data)
taylor_liver.exp.GOI<-taylor_liver.exp[genes,]
taylor_liver.exp.GOI$gene<-rownames(taylor_liver.exp.GOI)
taylor_liver.exp.GOI<-melt(taylor_liver.exp.GOI)#

meta<-taylor_liver@meta.data
meta$cell<-rownames(meta)

plt<-merge(taylor_liver.exp.GOI, meta,by.x="variable", by.y="cell")

plt$variable<-as.character(plt$variable)

## possible cells types
cluster_marker_mean<-function(gene_list, type){
  plt_epi<-plt[which(plt$gene%in%gene_list),]
  mean_type<-as.data.frame(tapply(plt_epi$value, plt_epi$seurat_clusters, mean))
  colnames(mean_type)<-type
  mean_type
}

cell_rough<-cbind(cluster_marker_mean(Macrophage_genes, "Myeloid"),
                  cluster_marker_mean(NK_T_genes, "NK_T"),
                  cluster_marker_mean(LEC_genes, "LEC"),
                  cluster_marker_mean(Hepatocyte_genes, "Hepatocyte"),
                  cluster_marker_mean(Cholangiocytes_genes, "Cholangiocytes"),
                  cluster_marker_mean(HSCs_genes, "HSC"),
                  cluster_marker_mean(B_genes, "B_cell"))


cell_rough$CellType_rough<-sapply(1:nrow(cell_rough), function(x) {
  compart<-colnames(cell_rough)[which(cell_rough[x,] == max(cell_rough[x,]))]
  if(length(compart)==1){compart}else{"Unclear"}
})

cell_rough$seurat_clusters<-rownames(cell_rough)




meta<-taylor_liver@meta.data
meta$cell<-rownames(meta)
plt_summary<-merge(meta, cell_rough[,c("seurat_clusters","CellType_rough")], by="seurat_clusters")
plt_summary<-plt_summary[match(rownames(taylor_liver@meta.data),plt_summary$cell),]
identical(plt_summary$cell, rownames(taylor_liver@meta.data))

rownames(plt_summary)<-plt_summary$cell

taylor_liver<- AddMetaData(taylor_liver, plt_summary)
taylor_liver
table(taylor_liver$CellType_rough)

taylor_liver@meta.data$CellType_rough<-as.factor(taylor_liver@meta.data$CellType_rough)
levels(taylor_liver@meta.data$CellType_rough)<-c("B-cells","Myeloid cells","NK and T cells")

DimPlot(taylor_liver, reduction = "umap",group.by="CellType_rough", pt.size=0.15, label=T)+
  colscale_cellType+ggtitle("")+annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(taylor_liver))))


#################
## SCINA
#################
DefaultAssay(taylor_liver)<-"RNA"
taylor_liver <- NormalizeData(taylor_liver,scale.factor = 10000, normalization.method = "LogNormalize")

signatures<-read.csv(here("data/Liver_Markers - Human_for_SCINA.csv"))

taylor_liver_exp <- GetAssayData(taylor_liver)
results = SCINA(taylor_liver_exp, signatures, max_iter = 100, convergence_n = 10,
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
taylor_liver$SCINA_broad<-results$cell_labels
gc()

no_refined<-data.frame(cell=colnames(taylor_liver), SCINA_broad=results$cell_labels, SCINA_refined=NA)


### Cell subsets
RBCsignatures<-read.csv(here("data/Liver_Markers - Erythrocytes.csv"))
Neutrosignatures<-read.csv(here("data/Liver_Markers - Neurtophil.csv"))
Tsignatures<-read.csv(here("data/Liver_Markers - Tcell.csv"))
Bsignatures<-read.csv(here("data/Liver_Markers - Bcell.csv"))
myeloidsignatures<-read.csv(here("data/Liver_Markers - Myeloid.csv"))


taylor_liver.combined_myeloid<-subset(taylor_liver, subset = SCINA_broad == "Myeloid")
taylor_liver_exp <- GetAssayData(taylor_liver.combined_myeloid)
results = SCINA(taylor_liver_exp, myeloidsignatures, max_iter = 100, convergence_n = 10,
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
taylor_liver.combined_myeloid$SCINA_refined<-results$cell_labels


taylor_liver.combined_tcell<-subset(taylor_liver, subset = SCINA_broad == "T_cell")
taylor_liver_exp <- GetAssayData(taylor_liver.combined_tcell)
results = SCINA(taylor_liver_exp, Tsignatures, max_iter = 100, convergence_n = 10,
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
taylor_liver.combined_tcell$SCINA_refined<-results$cell_labels

taylor_liver.combined_neutro<-subset(taylor_liver, subset = SCINA_broad == "Neutrophil")
taylor_liver_exp <- GetAssayData(taylor_liver.combined_neutro)
results = SCINA(taylor_liver_exp, Neutrosignatures, max_iter = 100, convergence_n = 10,
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
taylor_liver.combined_neutro$SCINA_refined<-results$cell_labels



SCINA_cell_labels<-rbind(#taylor_liver.combined_RBC@meta.data[,c("cell","SCINA_broad","SCINA_refined")],taylor_liver.combined_neutro@meta.data[,c("cell","SCINA_broad","SCINA_refined")]
                         taylor_liver.combined_myeloid@meta.data[,c("cell","SCINA_broad","SCINA_refined")],
                         taylor_liver.combined_neutro@meta.data[,c("cell","SCINA_broad","SCINA_refined")],
                         taylor_liver.combined_tcell@meta.data[,c("cell","SCINA_broad","SCINA_refined")])
SCINA_cell_labels$cell<-rownames(SCINA_cell_labels)

SCINA_cell_labels<-rbind(SCINA_cell_labels, no_refined[which(!(no_refined$cell%in%SCINA_cell_labels$cell)),])



length(which(SCINA_cell_labels$cell%in%colnames(taylor_liver)))

SCINA_cell_labels<-SCINA_cell_labels[match(colnames(taylor_liver), SCINA_cell_labels$cell),]
identical(colnames(taylor_liver), SCINA_cell_labels$cell)
rownames(SCINA_cell_labels)<-SCINA_cell_labels$cell
taylor_liver <- AddMetaData(taylor_liver, metadata = SCINA_cell_labels)

taylor_liver$SCINA_refined[which(is.na(taylor_liver$SCINA_refined))]<-taylor_liver$SCINA_broad[which(is.na(taylor_liver$SCINA_refined))]

DimPlot(taylor_liver, reduction = "umap",group.by="SCINA_refined", pt.size=0.15, label=T)+
  colscale_cellType+ggtitle("")+
  annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(taylor_liver))))


taylor_liver$SCINA_refined<-as.factor(taylor_liver$SCINA_refined)
levels(taylor_liver$SCINA_refined)<-c("B-cells","CD3+ T-cells","Cholangiocytes","Erythrocytes",
                                      "gd T-cells","Hepatocytes",
                                       "HSC","KC Like","LSEC",
                                       "Macrophage\n(MHCII high)","Neutrophil","NK-like cells","Platelets",
                                       "RR Myeloid","Unknown")

SCINA_cellUMAP<-DimPlot(taylor_liver, reduction = "umap",group.by="SCINA_refined", pt.size=0.15, label=F)+
  colscale_cellType+ggtitle("")+
  annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(taylor_liver))))
SCINA_cellUMAP



#############
## Annotate clusters by majority cells in cluster
#############
cluster_level<-table(taylor_liver$seurat_clusters, taylor_liver$SCINA_refined)
cluster_level
cluster_level_percent<-as.data.frame.matrix(signif((cluster_level/rowSums(cluster_level))*100,2))
cluster_level_percent$Unknown<-NULL

cluster_consensus<-data.frame(cluster=rownames(cluster_level_percent),cluster_consensus=sapply(1:nrow(cluster_level_percent), function(x) {
  celltype<-colnames(cluster_level_percent)[which(cluster_level_percent[x,]==max(cluster_level_percent[x,]) & which(cluster_level_percent[x,]>20))]
  if(length(celltype)==0){"unknown"}else{celltype}
}))

meta_add<-merge(taylor_liver@meta.data, cluster_consensus, by.x="seurat_clusters", by.y="cluster")

meta_add<-meta_add[match(colnames(taylor_liver), meta_add$cell),]
identical(colnames(taylor_liver), meta_add$cell)
rownames(meta_add)<-meta_add$cell
taylor_liver <- AddMetaData(taylor_liver, metadata = meta_add)

SCINA_cellUMAP<-DimPlot(taylor_liver, reduction = "umap",group.by="cluster_consensus", pt.size=0.15, label=F)+
  colscale_cellType+ggtitle("")+
  annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(taylor_liver))))
SCINA_cellUMAP


DefaultAssay(taylor_liver)<-"SCT"
FeaturePlot(taylor_liver, features = Macrophage_genes, min.cutoff = "q9", pt.size=0.15)
FeaturePlot(taylor_liver, features = LEC_genes, min.cutoff = "q9", pt.size=0.15)
FeaturePlot(taylor_liver, features = c("HBB","HBA2","HBA1","FCGR3A"), min.cutoff = "q9", pt.size=0.15)
FeaturePlot(taylor_liver, features = c("COL1A1","HBA2","HBA1","FCGR3A"), min.cutoff = "q9", pt.size=0.15)
FeaturePlot(taylor_liver, features = c("MKI67","PTPRC","IGFBP7","FCGR2B"), min.cutoff = "q9", pt.size=0.15)
FeaturePlot(taylor_liver, features = Cholangiocytes_genes, min.cutoff = "q9", pt.size=0.15)


DefaultAssay(taylor_liver) <- "RNA"
DotPlot(object = taylor_liver, features = B_genes)+xlab("B Cell Marker")
DotPlot(object = taylor_liver, features = T_genes)+xlab("T Cell Marker")
DotPlot(object = taylor_liver, features = NK_genes)+xlab("NK Cell Marker")
DotPlot(object = taylor_liver, features = LEC_genes)+xlab("LSEC Marker")
DotPlot(object = taylor_liver, features = Hepatocyte_genes[1:15])+xlab("Hepatocyte Marker")
DotPlot(object = taylor_liver, features = Hepatocyte_genes[16:29])+xlab("Hepatocyte Marker")
DotPlot(object = taylor_liver, features = Cholangiocytes_genes)+xlab("Cholangiocyte Marker")
DotPlot(object = taylor_liver, features = HSCs_genes[1:14])+xlab("HSC Marker")
DotPlot(object = taylor_liver, features = HSCs_genes[15:27])+xlab("HSC Marker")
DotPlot(object = taylor_liver, features = Macrophage_genes)+xlab("Macrophage Marker")


## unclear clusters
DefaultAssay(taylor_liver) <- "RNA"

cluster4.markers <-  FindMarkers(taylor_liver, ident.1 = 4,ident.2 = 3, min.pct = 0.25)
head(cluster4.markers, n = 10)
head(cluster4.markers[which(cluster4.markers$avg_log2FC>0),], n = 20)
head(cluster4.markers[which(cluster4.markers$avg_log2FC<0),], n = 20)

FeaturePlot(taylor_liver, features = c("XIST","NEAT1","TPT1","CCNL1"), min.cutoff = "q9", pt.size=1)
FeaturePlot(taylor_liver, features = c("CALCRL","LYVE1","ID1","MEG3"), min.cutoff = "q9", pt.size=1)

cluster12.markers <-  FindMarkers(taylor_liver, ident.1 = 12, min.pct = 0.25)
head(cluster12.markers, n = 10)
head(cluster12.markers[which(cluster12.markers$avg_log2FC>0),], n = 20)
FeaturePlot(taylor_liver, features = c("NEAT1","S100A16","IGFBP7","SEPP1"), min.cutoff = "q9", pt.size=1)
FeaturePlot(taylor_liver, features = c("CFH","LINC00299","KRT86","SPINK2"), min.cutoff = "q9", pt.size=1)


## relabel some
taylor_liver@meta.data$cluster_consensus[which(taylor_liver@meta.data$seurat_clusters%in%c("3"))]<-"LSEC"

SCINA_cellUMAP<-DimPlot(taylor_liver, reduction = "umap",group.by="cluster_consensus", pt.size=0.15, label=F)+
  colscale_cellType+ggtitle("")+
  annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(taylor_liver))))
SCINA_cellUMAP
