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


dataset_loc <- here("/media/redgar/Seagate Portable Drive/Caudate_Flush")


caudate_10x<- Read10X(file.path(dataset_loc,"Caudate","filtered_gene_bc_matrices/GRCh38"))
colnames(caudate_10x) <- paste(sapply(strsplit(colnames(caudate_10x),split="-"),'[[',1L),"caudate",sep="-")
caudate_10x<-CreateSeuratObject(counts = caudate_10x, project = "caudate", min.cells = 0, min.features = 0)

flush_10x<- Read10X(file.path(dataset_loc,"Flush","filtered_gene_bc_matrices/GRCh38"))
colnames(flush_10x) <- paste(sapply(strsplit(colnames(flush_10x),split="-"),'[[',1L),"flush",sep="-")
flush_10x<-CreateSeuratObject(counts = flush_10x, project = "flush", min.cells = 0, min.features = 0)

d10x.list<-list(flush_10x, caudate_10x)

## cell counts
plt_count_raw<-lapply(1:length(d10x.list), function(x) {
  df<-data.frame(raw_cell_count=nrow(d10x.list[[x]]@meta.data),sample=unique(d10x.list[[x]]@meta.data$orig.ident))
  df})
plt_count_raw<-do.call(rbind, plt_count_raw)
print(plt_count_raw)

#'## QC
#'The percentage of reads that map to the mitochondrial genome
#'Low-quality / dying cells often exhibit extensive mitochondrial contamination
#'We calculate mitochondrial QC metrics with the PercentageFeatureSet function, which calculates the percentage of counts originating from a set of features
#'We use the set of all genes starting with MT- as a set of mitochondrial genes

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
invisible(lapply(1:length(d10x.list), function(x){
  d10x.list[[x]][["percent.mt"]] <<- PercentageFeatureSet(d10x.list[[x]], pattern = "^MT-")}))

# Show QC metrics for the first 5 cells
print(head(d10x.list[[2]]@meta.data, 5))

#'Low-quality cells or empty droplets will often have very few genes
#'Cell doublets or multiplets may exhibit an aberrantly high gene count

# Visualize QC metrics
#nFeature number of unique genes
#nCount number of total molecules
plt_QC_data<-do.call(rbind, lapply(1:length(d10x.list), function(x) d10x.list[[x]]@meta.data))

qc_plts<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
  geom_point() +
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th

qc_plts_sample<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
  geom_point() + facet_wrap(~orig.ident)+
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th
save_plts(qc_plts_sample, "flush_intital_QC_plts_sample", w=12,h=4)

MT_plt_orig.ident<-ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
  facet_wrap(~orig.ident, scales="free_y")+
  geom_vline(xintercept = 25)+ theme_bw()+xlab("Percent Mitochondrial")+th
save_plts(MT_plt_orig.ident, "flush_percentMT_plt_sample", w=8,h=4)


#'We filter cells that have unique feature counts over 6,000 or less than 500
#'We filter cells that have >10% mitochondrial counts
#'we will also filter doublets as called by scrublet

invisible(lapply(1:length(d10x.list), function(x){
  d10x.list[[x]] <<- subset(d10x.list[[x]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 25)
}))

d10x.list

## cell counts after QC
plt_count_QC<-lapply(1:length(d10x.list), function(x) {
  df<-data.frame(qc_cell_count=nrow(d10x.list[[x]]@meta.data),orig.ident=unique(d10x.list[[x]]@meta.data$orig.ident))
  df})
plt_count_QC<-do.call(rbind, plt_count_QC)
print(plt_count_QC)

meta<-merge(plt_count_raw, plt_count_QC, by.x="sample", by.y="orig.ident")

cell_count<-grid.arrange(ggplot(meta, aes(sample, raw_cell_count))+
                           geom_boxplot()+geom_point()+
                           theme_bw()+
                           ylab("Total Cell Number")+th+fillscale_age+ylim(0,60000)+
                           theme(legend.position = "none")+ggtitle("Before Quality Control"),
                         ggplot(meta, aes(sample, qc_cell_count))+
                           geom_boxplot()+geom_point()+
                           theme_bw()+
                           ylab("Total Cell Number")+th+fillscale_age+ylim(0,60000)+
                           theme(legend.position = "none")+ggtitle("After Quality Control"), ncol=2)

save_plts(cell_count, "flush_QC_cellcount_age", w=8,h=4)


d10x <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "caudate_flush")#add.cell.ids = alldata_names2,

d10x


################
## Normalize scale and UMAP
################
d10x <- NormalizeData(d10x)
d10x <- FindVariableFeatures(d10x, selection.method = "vst", nfeatures = 2000)
d10x <- ScaleData(d10x) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))

# dimension reduction
d10x <- RunPCA(d10x, ndims.print = 1:10, nfeatures.print = 10)
d10x <- RunUMAP(d10x, dims = 1:30)
d10x <- RunTSNE(d10x, dims = 1:30)



######################
## cell cycle gene expression
######################
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

d10x <- CellCycleScoring(d10x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

pca_cellcycle<-DimPlot(d10x, reduction="pca",  group.by = "Phase")
pca_nfeature<-FeaturePlot(d10x, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)



## regress out cell cycle and other covariates
#Transformed data will be available in the SCT assay, which is set as the default after running sctransform
#By default, sctransform accounts for cellular sequencing depth, or nUMIs.
d10x <- SCTransform(d10x, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), verbose = FALSE)

# dimension reduction
d10x <- RunPCA(d10x, verbose = FALSE)
d10x <- RunUMAP(d10x, dims = 1:30)
d10x <- RunTSNE(d10x, dims = 1:30)

# cluster
d10x <- FindNeighbors(d10x, reduction = "pca", dims = 1:20)
d10x <- FindClusters(d10x, resolution = 0.5)



###############
## visualize
###############
SCT_cluster_umap<-DimPlot(d10x, reduction = "umap", pt.size=0.25, label=T)
save_plts(SCT_cluster_umap, "flush_SCT_cluster_umap", w=6,h=4)

cell_pca_SCT<-DimPlot(d10x, reduction="pca", group.by="Phase")
save_plts(cell_pca_SCT, "flush_cell_PCA_afterSCT", w=6,h=4)

nFeature_UMAP_SCT<-FeaturePlot(d10x, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
save_plts(nFeature_UMAP_SCT, "flush_nfeature_UMAP_afterSCT", w=6,h=4)

chem_umap_sct<-DimPlot(d10x, reduction = "umap", group.by = "orig.ident", pt.size=0.25)
save_plts(chem_umap_sct, "flush_sample_SCT_umap", w=6,h=4)



## Low quality cluster?
MT_umap_SCT<-FeaturePlot(d10x, features = "percent.mt", min.cutoff = "q9", pt.size=1, split.by="orig.ident")
save_plts(pca_nfeature, "flush_pca_nfeature", w=6,h=4)

ncount_umap_SCT<-FeaturePlot(d10x, features = "nCount_RNA", min.cutoff = "q9", pt.size=1)
save_plts(ncount_umap_SCT, "flush_ncount_umap_SCT", w=6,h=4)

nfeature_umap_SCT<-FeaturePlot(d10x, features = "nFeature_RNA", min.cutoff = "q9", pt.size=1)
save_plts(nfeature_umap_SCT, "flush_nfeature_umap_SCT", w=6,h=4)



















#########################
## that is a very bad caudate, super high MT, So I am just going to look at the flush
#########################
flush_10x<- Read10X(file.path(dataset_loc,"Flush","filtered_gene_bc_matrices/GRCh38"))
colnames(flush_10x) <- paste(sapply(strsplit(colnames(flush_10x),split="-"),'[[',1L),"flush",sep="-")
flush_10x<-CreateSeuratObject(counts = flush_10x, project = "flush", min.cells = 0, min.features = 0)

flush_10x[["percent.mt"]] <- PercentageFeatureSet(flush_10x, pattern = "^MT-")
flush_10x_raw <- subset(flush_10x, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 25)

# Run the standard workflow for visualization and clustering
flush_10x <- NormalizeData(flush_10x_raw)
flush_10x <- FindVariableFeatures(flush_10x, selection.method = "vst", nfeatures = 2000)
flush_10x <- ScaleData(flush_10x) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))

flush_10x <- RunPCA(flush_10x, npcs = 30, verbose = FALSE)
flush_10x <- RunUMAP(flush_10x, reduction = "pca", dims = 1:30)
flush_10x <- RunTSNE(flush_10x, dims = 1:30)

flush_10x <- FindNeighbors(flush_10x, reduction = "pca", dims = 1:30)
flush_10x <- FindClusters(flush_10x, resolution = 0.5)

flush_10x

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
flush_10x <- CellCycleScoring(flush_10x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

###########
## Visualize integration
###########
SCT_cluster_umap<-DimPlot(flush_10x, reduction = "umap", pt.size=0.25, label=T)
save_plts(SCT_cluster_umap, "flush_only_cluster_umap", w=6,h=4)

MT_umap_sct<-FeaturePlot(flush_10x, reduction = "umap", features = "percent.mt", pt.size=0.25)
save_plts(MT_umap_sct, "flush_MT_only_umap", w=5,h=4)



############################
##### Example markers to plot for the liver [From Diana]
############################

#' #################
#' ## Rough annotation
#' #################
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


#################
## SCINA
#################
d10_flush_raw <- NormalizeData(d10_flush_raw,scale.factor = 10000, normalization.method = "LogNormalize")

signatures<-read.csv(here("data/Liver_Markers_with_citations - Human_for_SCINA.csv"))

d10x_exp <- GetAssayData(d10_flush_raw)
results = SCINA(d10x_exp, signatures, max_iter = 100, convergence_n = 10,
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
d10_flush_raw$SCINA_broad<-results$cell_labels
gc()

no_refined<-data.frame(cell=colnames(d10_flush_raw), SCINA_broad=results$cell_labels, SCINA_refined=NA)


### Cell subsets
## Cell subsets
RBCsignatures<-read.csv(here("data/Liver_Markers_with_citations - Erythrocytes.csv"))
Neutrosignatures<-read.csv(here("data/Liver_Markers_with_citations - Neutrophil.csv"))
Tsignatures<-read.csv(here("data/Liver_Markers_with_citations - Tcell.csv"))
Bsignatures<-read.csv(here("data/Liver_Markers_with_citations - Bcell.csv"))
myeloidsignatures<-read.csv(here("data/Liver_Markers_with_citations - Myeloid.csv"))

d10_flush_raw$cell<-rownames(d10_flush_raw@meta.data)

# d10_flush_RBC<-subset(d10_flush, subset = SCINA_broad == "Erythrocytes")
# d10x_exp <- GetAssayData(d10_flush_RBC)
# results = SCINA(d10x_exp, RBCsignatures, max_iter = 100, convergence_n = 10,
#                 convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
# d10_flush_RBC$SCINA_refined<-results$cell_labels

d10_flush_myeloid<-subset(d10_flush_raw, subset = SCINA_broad == "Myeloid")
d10x_exp <- GetAssayData(d10_flush_myeloid)
results = SCINA(d10x_exp, myeloidsignatures, max_iter = 100, convergence_n = 10,
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
d10_flush_myeloid$SCINA_refined<-results$cell_labels

d10_flush_bcell<-subset(d10_flush_raw, subset = SCINA_broad == "B_cell")
d10x_exp <- GetAssayData(d10_flush_bcell)
results = SCINA(d10x_exp, Bsignatures, max_iter = 100, convergence_n = 10,
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
d10_flush_bcell$SCINA_refined<-results$cell_labels

d10_flush_tcell<-subset(d10_flush_raw, subset = SCINA_broad == "T_cell")
d10x_exp <- GetAssayData(d10_flush_tcell)
results = SCINA(d10x_exp, Tsignatures, max_iter = 100, convergence_n = 10,
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
d10_flush_tcell$SCINA_refined<-results$cell_labels

d10_flush_neutro<-subset(d10_flush_raw, subset = SCINA_broad == "Neutrophil")
d10x_exp <- GetAssayData(d10_flush_neutro)
results = SCINA(d10x_exp, Neutrosignatures, max_iter = 100, convergence_n = 10,
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
d10_flush_neutro$SCINA_refined<-results$cell_labels



SCINA_cell_labels<-rbind(#d10_flush_RBC@meta.data[,c("cell","SCINA_broad","SCINA_refined")],
                         d10_flush_neutro@meta.data[,c("cell","SCINA_broad","SCINA_refined")],
                         d10_flush_myeloid@meta.data[,c("cell","SCINA_broad","SCINA_refined")],
                         d10_flush_bcell@meta.data[,c("cell","SCINA_broad","SCINA_refined")],
                         d10_flush_tcell@meta.data[,c("cell","SCINA_broad","SCINA_refined")])
SCINA_cell_labels$cell<-rownames(SCINA_cell_labels)

SCINA_cell_labels<-rbind(SCINA_cell_labels, no_refined[which(!(no_refined$cell%in%SCINA_cell_labels$cell)),])


length(which(SCINA_cell_labels$cell%in%colnames(flush_10x)))

SCINA_cell_labels<-SCINA_cell_labels[match(colnames(flush_10x), SCINA_cell_labels$cell),]
identical(colnames(flush_10x), SCINA_cell_labels$cell)
rownames(SCINA_cell_labels)<-SCINA_cell_labels$cell
flush_10x <- AddMetaData(flush_10x, metadata = SCINA_cell_labels)

flush_10x$SCINA_refined[which(is.na(flush_10x$SCINA_refined))]<-flush_10x$SCINA_broad[which(is.na(flush_10x$SCINA_refined))]

DimPlot(flush_10x, reduction = "umap",group.by="SCINA_refined", pt.size=0.15, label=T)+
  ggtitle("")+
  annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(flush_10x))))



flush_10x$SCINA_refined<-as.factor(flush_10x$SCINA_refined)
levels(flush_10x$SCINA_refined)<-c("CD3+ T-cells","Cholangiocytes","Neutrophil","Erythrocytes","gd T-cells","Hepatocytes",
                                   "HSC","KC Like","LSEC",
                                   "Mature B-cells","Macrophage\n(MHCII high)","Neutrophil","NK-like cells",
                                   "RR Myeloid","Unknown")

SCINA_cellUMAP<-DimPlot(flush_10x, reduction = "umap",group.by="SCINA_refined", pt.size=0.15, label=T)+
  colscale_cellType+ggtitle("")+
  annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(flush_10x))))
SCINA_cellUMAP
save_plts(SCINA_cellUMAP, "flush_rPCA_SCINA_umap", w=6,h=4)


#############
## Annotate clusters by majority cells in cluster
#############
cluster_level<-table(flush_10x$seurat_clusters, flush_10x$SCINA_refined)
cluster_level
cluster_level_percent<-as.data.frame.matrix(signif((cluster_level/rowSums(cluster_level))*100,2))
cluster_level_percent$Unknown<-NULL

cluster_consensus<-data.frame(cluster=rownames(cluster_level_percent),cluster_consensus=sapply(1:nrow(cluster_level_percent), function(x) {
  celltype<-colnames(cluster_level_percent)[which(cluster_level_percent[x,]==max(cluster_level_percent[x,]) & which(cluster_level_percent[x,]>20))]
  if(length(celltype)==0){"unknown"}else{celltype}
}))

meta_add<-merge(flush_10x@meta.data, cluster_consensus, by.x="seurat_clusters", by.y="cluster")

meta_add<-meta_add[match(colnames(flush_10x), meta_add$cell),]
identical(colnames(flush_10x), meta_add$cell)
rownames(meta_add)<-meta_add$cell
flush_10x <- AddMetaData(flush_10x, metadata = meta_add)




########
## refining cell labels
########
Idents(flush_10x)<-"seurat_clusters"
DefaultAssay(flush_10x) <- "RNA"

pdf(file = here("figures/flush_dot_plots.pdf"), w=10, h=10)
DotPlot(object = flush_10x, features = B_genes)+xlab("B Cell Marker")
DotPlot(object = flush_10x, features = T_genes)+xlab("T Cell Marker")
DotPlot(object = flush_10x, features = NK_genes)+xlab("NK Cell Marker")
DotPlot(object = flush_10x, features = LEC_genes)+xlab("LSEC Marker")
DotPlot(object = flush_10x, features = Hepatocyte_genes[1:15])+xlab("Hepatocyte Marker")
DotPlot(object = flush_10x, features = Hepatocyte_genes[16:29])+xlab("Hepatocyte Marker")
DotPlot(object = flush_10x, features = Cholangiocytes_genes)+xlab("Cholangiocyte Marker")
DotPlot(object = flush_10x, features = HSCs_genes[1:14])+xlab("HSC Marker")
DotPlot(object = flush_10x, features = HSCs_genes[15:27])+xlab("HSC Marker")
DotPlot(object = flush_10x, features = Macrophage_genes)+xlab("Macrophage Marker")
dev.off()



DimPlot(flush_10x, reduction = "umap")+scale_color_manual(values=c(rep("grey",13),"red",rep("grey", 19)))

FeaturePlot(flush_10x, features = gd_genes)
DimPlot(flush_10x, reduction = "umap", group.by = "Phase")


#############
# #relabel some clusters
#############
flush_10x@meta.data$CellType_refined<-as.character(flush_10x@meta.data$cluster_consensus)
flush_10x@meta.data$CellType_refined[which(flush_10x@meta.data$seurat_clusters=="10")]<-"Plasma cells"
flush_10x@meta.data$CellType_refined[which(flush_10x@meta.data$seurat_clusters=="2")]<-"Naive CD4 T-cells"
flush_10x@meta.data$CellType_refined[which(flush_10x@meta.data$seurat_clusters%in%c(5,6,7))]<-"Cycling T-cells"


DimPlot(flush_10x, reduction = "umap",group.by="CellType_refined", pt.size=0.15, label=T)+colscale_cellType+ggtitle("")+
  annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(flush_10x))))


## what is 13
cluster.markers <-  FindMarkers(flush_10x, ident.1 = 13, min.pct = 0.25)
head(cluster.markers, n = 10)
head(cluster.markers[which(cluster.markers$avg_log2FC>0),], n = 20)
FeaturePlot(flush_10x, features = c("CD1C","IDO1","FCER1A","CLEC9A"), min.cutoff = "q9", pt.size=1)








##############
## Myeloid clustering
##############
d10_flush_myeloid<-subset(flush_10x, subset = CellType_refined %in% c("Macrophage\n(MHCII high)"))
d10_flush_myeloid <- RunPCA(d10_flush_myeloid, npcs = 30, verbose = FALSE)
d10_flush_myeloid <- RunUMAP(d10_flush_myeloid, reduction = "pca", dims = 1:30)
d10_flush_myeloid <- FindNeighbors(d10_flush_myeloid, reduction = "pca", dims = 1:30)
d10_flush_myeloid <- FindClusters(d10_flush_myeloid, resolution = 0.5)

myeloid_cluster_umap<-DimPlot(d10_flush_myeloid, reduction = "umap", pt.size=0.25, label=T)
myeloid_cluster_umap
save_plts(myeloid_cluster_umap, "flush_myeloid_cluster_umap", w=5,h=4)

recent_recruit_myeloid<-c("S100A8","S100A9","CD68","LYZ")
kuffer_signature<-c("VSIG4","MARCO","CD5L","HMOX1")
neutro_gene<-c("CSF3R","FCGR3B","NAMPT","CXCR2","DEFA3","DEFA4")
MHCII<-c("HLA-DRA","CD74","HLA-DPB1","HLA-DQB1")
LSEC<-c("CALCRL","STAB2","FCN2","FCN3")


FeaturePlot(d10_flush_myeloid, features = recent_recruit_myeloid, min.cutoff = "q9", pt.size=0.25)
FeaturePlot(d10_flush_myeloid, features = kuffer_signature, min.cutoff = "q9", pt.size=0.25)
FeaturePlot(d10_flush_myeloid, features = neutro_gene, min.cutoff = "q9", pt.size=0.25)
FeaturePlot(d10_flush_myeloid, features = MHCII, min.cutoff = "q9", pt.size=0.25)
FeaturePlot(d10_flush_myeloid, features = LSEC, min.cutoff = "q9", pt.size=0.25)



## one unclear cluster
cluster.markers <-  FindMarkers(d10_flush_myeloid, ident.1 = 2, min.pct = 0.25)
head(cluster.markers, n = 10)
head(cluster.markers[which(cluster.markers$avg_log2FC>0),], n = 20)
FeaturePlot(d10_flush_myeloid, features = c("CLEC9A","CD1C","FCER1A"), min.cutoff = "q9", pt.size=1)


d10_flush_myeloid@meta.data$CellType_refined<-as.character(d10_flush_myeloid@meta.data$CellType_refined)
d10_flush_myeloid@meta.data$CellType_refined[which(d10_flush_myeloid@meta.data$seurat_clusters%in%c("0"))]<-"RR Myeloid"
d10_flush_myeloid@meta.data$CellType_refined[which(d10_flush_myeloid@meta.data$seurat_clusters%in%c("1"))]<-"Macrophage (MHCII high)"
d10_flush_myeloid@meta.data$CellType_refined[which(d10_flush_myeloid@meta.data$seurat_clusters%in%c("2"))]<-"DC"



##############
## Relabel subtypes
##############

flush_10x@meta.data$CellType_refined<-sapply(1:nrow(flush_10x@meta.data), function(x){
  if(rownames(flush_10x@meta.data)[x]%in%rownames(d10_flush_myeloid@meta.data)){
    d10_flush_myeloid@meta.data$CellType_refined[which(rownames(d10_flush_myeloid@meta.data)==rownames(flush_10x@meta.data)[x])]
        }else{flush_10x@meta.data$CellType_refined[x]}
      })

flush_10x@meta.data$CellType_refined<-as.factor(flush_10x@meta.data$CellType_refined)
print(levels(flush_10x@meta.data$CellType_refined))
levels(flush_10x@meta.data$CellType_refined)[which(levels(flush_10x@meta.data$CellType_refined)=="Macrophage (MHCII high)")]<-"Macrophage\n(MHCII high)"

all_refined_cluster_umap<-DimPlot(flush_10x, reduction = "umap", pt.size=0.25, label=T,label.size = 3, group.by = "CellType_refined")+
  colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-13, y=-13, label=paste0("n = ",comma(ncol(flush_10x))))
all_refined_cluster_umap
save_plts(all_refined_cluster_umap, "flush_refined_cellType_map", w=12,h=8)


all_refined_cluster_umap_nolab<-DimPlot(flush_10x, reduction = "umap", pt.size=0.25, group.by = "CellType_refined")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-13, y=-13, label=paste0("n = ",comma(ncol(flush_10x))))
all_refined_cluster_umap_nolab
save_plts(all_refined_cluster_umap_nolab, "flush_refined_cellType_map_nolabel", w=12,h=8)

fancyUMAP_all<-fanciest_UMAP(flush_10x,NA,F)
save_plts(fancyUMAP_all, "flush_refined_cellType_umpa_fancy", w=6,h=4)



##############
## Save integrated with refined cluster labels
##############
save(flush_10x, file=paste(here("/media/redgar/Seagate Portable Drive/processed_data/"),"flush_refinedlabels.rds", sep=""))
cell_label<-flush_10x@meta.data
save(cell_label, file=paste(here("/media/redgar/Seagate Portable Drive/processed_data/"),"flush_cell_labels_Refined.rds", sep=""))



print(sessionInfo())




