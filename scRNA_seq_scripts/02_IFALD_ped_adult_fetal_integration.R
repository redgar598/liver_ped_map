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
library(anndata)
library(RColorBrewer)

source("scRNA_seq_scripts/00_pretty_plots.R")
source("scRNA_seq_scripts/00_entropy_d10x.R")

# #################
# ## Load raw QC'ed data
# #################
# d10x_fetal<-readRDS(here("../../../projects/macparland/RE/PediatricAdult/processed_data","d10x_fetal_raw.rds"))
# #d10x_fetal<-readRDS("/media/redgar/Seagate Portable Drive/fetal_liver/d10x_fetal_raw.rds")
# 
# d10x_ped_IFALD<-readRDS(file = here("../../../projects/macparland/RE/PediatricAdult/processed_data","IFALD_d10x_adult_ped_raw.rds"))
# ## add cell type labels
# load(here("../../../projects/macparland/RE/PediatricAdult/processed_data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
# cell_label$index<-rownames(cell_label)
# cell_label<-cell_label[match(colnames(d10x_ped_IFALD), cell_label$index),]
# identical(colnames(d10x_ped_IFALD), cell_label$index)
# d10x_ped_IFALD <- AddMetaData(d10x_ped_IFALD, metadata = cell_label)
# 
# 
# #################
# ## Match meta data
# #################
# d10x_fetal@meta.data$Barcodes<-NULL
# colnames(d10x_fetal@meta.data)[which(colnames(d10x_fetal@meta.data)=="Extract.Name")]<-"Sample"
# colnames(d10x_fetal@meta.data)[which(colnames(d10x_fetal@meta.data)=="Characteristics.individual.")]<-"individual"
# d10x_fetal@meta.data$Treatment<-"Healthy"
# d10x_fetal@meta.data$Tissue<-"TLH"
# d10x_fetal@meta.data$chemistry<-"3pr"
# colnames(d10x_fetal@meta.data)[which(colnames(d10x_fetal@meta.data)=="Characteristics.sex.")]<-"Sex"
# colnames(d10x_fetal@meta.data)[which(colnames(d10x_fetal@meta.data)=="Characteristics.age.")]<-"Age"
# d10x_fetal@meta.data$AgeGroup<-"Fetal"
# d10x_fetal@meta.data$FreshorFrozen<-"fresh"
# d10x_fetal@meta.data$BMI<-NA
# d10x_fetal@meta.data$relALBChange<-NA
# d10x_fetal@meta.data$nuclear_fraction<-NA
# d10x_fetal@meta.data$cell_status<-NA
# d10x_fetal@meta.data$Perfused<-NA
# colnames(d10x_fetal@meta.data)[which(colnames(d10x_fetal@meta.data)=="Cell.Labels")]<-"CellType_refined"
# d10x_fetal@meta.data$age_condition<-paste(d10x_fetal$AgeGroup, d10x_fetal$Treatment, sep=" ")
# 
# 
# d10x_ped_IFALD@meta.data$file<-NULL
# d10x_ped_IFALD@meta.data$Approx_bam_GB<-NULL
# d10x_ped_IFALD@meta.data$Characteristics.facs.sorting.<-NA
# d10x_ped_IFALD@meta.data$Sample<-d10x_ped_IFALD@meta.data$individual
# d10x_ped_IFALD@meta.data$individual<-sapply(1:nrow(d10x_ped_IFALD@meta.data), function(x) strsplit(d10x_ped_IFALD@meta.data$individual[x],"_")[[1]][1])
# d10x_ped_IFALD@meta.data$CellType_rough <-NULL
# d10x_ped_IFALD@meta.data$second_best_cell  <-NULL
# d10x_ped_IFALD@meta.data$S.Score   <-NULL
# d10x_ped_IFALD@meta.data$G2M.Score <-NULL
# d10x_ped_IFALD@meta.data$Phase     <-NULL
# d10x_ped_IFALD@meta.data$old.ident<-NULL
# d10x_ped_IFALD@meta.data$integrated_snn_res.0.5 <-NULL
# d10x_ped_IFALD@meta.data$seurat_clusters<-NULL
# d10x_ped_IFALD@meta.data$age_id<-NULL
# d10x_ped_IFALD@meta.data$index<-NULL
# d10x_ped_IFALD@meta.data$SCINA_broad<-NULL
# d10x_ped_IFALD@meta.data$SCINA_refined<-NULL
# d10x_ped_IFALD@meta.data$cluster_consensus<-NULL
# d10x_ped_IFALD@meta.data$Characteristics.clinical.information.<-NA
# 
# 
# head(d10x_ped_IFALD@meta.data)
# head(d10x_fetal@meta.data)
# 
# d10x_fetal@meta.data<-d10x_fetal@meta.data[,colnames(d10x_ped_IFALD@meta.data)]
# 
# head(d10x_fetal@meta.data)
# head(d10x_ped_IFALD@meta.data)
# 
# 
# 
# 
# 
# 
# #################
# ## Merge
# #################
# d10x <- merge(d10x_ped_IFALD,d10x_fetal, merge.data=TRUE, project = "IFALD_fetal_adult_ped_map")
# saveRDS(d10x, file = here("../../../projects/macparland/RE/PediatricAdult/processed_data","Fetal_IFALD_d10x_adult_ped_raw.rds"))

d10x<-readRDS(here("../../../projects/macparland/RE/PediatricAdult/processed_data","Fetal_IFALD_d10x_adult_ped_raw.rds"))
unique(d10x$individual)


###############
## Integrate by donor
###############
#https://satijalab.org/seurat/articles/integration_rpca.html
print("RUNNING INTEGRATION")

## run integration across donor and hopefully that will also smooth out differences with chemistry and batch?
d10x.list<- SplitObject(d10x, split.by = "Sample")

# normalize, identify variable features and score cell cycle for each dataset independently
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

d10x.list <- lapply(X = d10x.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = d10x.list)
d10x.list <- lapply(X = d10x.list, FUN = function(x) {
  #x <- ScaleData(x, features = features, verbose = FALSE)
  x <- ScaleData(x, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

## Identify anchors
chem.anchors <- FindIntegrationAnchors(object.list = d10x.list, anchor.features = features, reduction = "rpca")
d10x.fetal_ped_IFALD <- IntegrateData(anchorset = chem.anchors)

DefaultAssay(d10x.fetal_ped_IFALD) <- "integrated"

print("INTEGRATED")
unique(d10x.fetal_ped_IFALD$individual)


# Run the standard workflow for visualization and clustering
d10x.fetal_ped_IFALD <- ScaleData(d10x.fetal_ped_IFALD, verbose = FALSE)
d10x.fetal_ped_IFALD <- RunPCA(d10x.fetal_ped_IFALD, npcs = 30, verbose = FALSE)
d10x.fetal_ped_IFALD <- RunUMAP(d10x.fetal_ped_IFALD, reduction = "pca", dims = 1:30)
d10x.fetal_ped_IFALD <- RunTSNE(d10x.fetal_ped_IFALD, dims = 1:30)

d10x.fetal_ped_IFALD <- FindNeighbors(d10x.fetal_ped_IFALD, reduction = "pca", dims = 1:30)
d10x.fetal_ped_IFALD <- FindClusters(d10x.fetal_ped_IFALD, resolution = 0.5)

d10x.fetal_ped_IFALD


##############
## Save integrated to look local
##############
save(d10x.fetal_ped_IFALD, file=paste(here("../../../projects/macparland/RE/PediatricAdult/processed_data/"),"Fetal_IFALD_adult_ped_integrated.rds", sep=""))
cell_label<-d10x.fetal_ped_IFALD@meta.data
save(cell_label, file=paste(here("../../../projects/macparland/RE/PediatricAdult/processed_data/"),"Fetal_IFALD_adult_ped_cellRough.rds", sep=""))


# load(file="/cluster/projects/macparland/RE/PediatricAdult/processed_data/Fetal_IFALD_adult_ped_integrated.rds")
# d10x<-readRDS(file="/cluster/projects/macparland/RE/PediatricAdult/processed_data/Fetal_IFALD_d10x_adult_ped_raw.rds")



###################
### save for local plotting
###################

Macrophage_genes<-c( "PTPRC", "MARCO","CD74")
LEC_genes<-c("CALCRL","RAMP2")
Hepatocyte_genes<-c("ALB", "CYP3A4")
Cholangiocytes_genes<-c( "EPCAM", "KRT7")
HSCs_genes<-c( "IGFBP7",  "SPARC")
T_genes<-c("CD3D","CD8A")
NK_genes<-c("NKG7","CD7")
gd_genes<-c("GNLY")
RBC<-c("HBB","HBA2","HBA1","FCGR3A")
MAST<-c("TPSAB1", "AREG")
recent_recruit_myeloid<-c("S100A8","S100A9","CD68","LYZ")
kuffer_signature<-c("VSIG4","CD5L")
neutro_gene<-c("CSF3R","FCGR3B")
MHCII<-c("HLA-DRA","HLA-DPB1")
b_genes_noIG<-c("MS4A1", "CD79B")
immunoglobins<-c("IGKC","IGHG1")

head(d10x.fetal_ped_IFALD)

umap_mat<-as.data.frame(Embeddings(object = d10x.fetal_ped_IFALD, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x.fetal_ped_IFALD@meta.data
meta$cell<-rownames(meta)

plt<-merge(meta, umap_mat, by="cell")

## raw expression values
gene_exp<-FetchData(d10x, vars=c(Macrophage_genes,LEC_genes,Hepatocyte_genes,Cholangiocytes_genes,HSCs_genes,T_genes,NK_genes,gd_genes,RBC,
                                 MAST, recent_recruit_myeloid, kuffer_signature, neutro_gene, MHCII, b_genes_noIG, immunoglobins))
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)

plt<-merge(gene_exp,plt, by="cell")


head(plt)

save(plt, file=here("../../../projects/macparland/RE/PediatricAdult/processed_data","Fetal_IFALD_adult_ped_pltData.RData"))



#############
## PCAs for plotting
#############
d10x$CellType_harmonized<-d10x$CellType_refined

levels(d10x$CellType_harmonized)[which(levels(d10x$CellType_harmonized)=="Kupffer Cell")]<-"KC Like"
levels(d10x$CellType_harmonized)[which(levels(d10x$CellType_harmonized)%in%c("Mono-Mac","Monocyte-DC precursor","Monocyte" ))]<-"Mono-Mac"
levels(d10x$CellType_harmonized)[which(levels(d10x$CellType_harmonized)%in%c("DC2","DC1" ))]<-"Macrophage\n(MHCII high)"
levels(d10x$CellType_harmonized)[which(levels(d10x$CellType_harmonized)%in%c("NK-like cells","NK" ))]<-"NK cell"

d10x.combined_myeloid<-subset(d10x, subset = CellType_harmonized %in% c("Mono-Mac","KC Like","Macrophage\n(MHCII high)","CDC1"))


d10x.combined_myeloid <- NormalizeData(d10x.combined_myeloid)
d10x.combined_myeloid <- FindVariableFeatures(d10x.combined_myeloid, selection.method = "vst", nfeatures = 2000)
d10x.combined_myeloid <- ScaleData(d10x.combined_myeloid)
d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)

d10x.combined_myeloid[["pca"]]

head(Embeddings(d10x.combined_myeloid, reduction = "pca")[, 1:5])
head(Loadings(d10x.combined_myeloid, reduction = "pca")[, 1:5])
head(Stdev(d10x.combined_myeloid, reduction = "pca"))

#' ## PCA for batch effect
Loadings<-as.data.frame(Loadings(d10x.combined_myeloid, reduction = "pca"))
embed<-as.data.frame(Embeddings(d10x.combined_myeloid, reduction = "pca"))
vars <- Stdev(d10x.combined_myeloid, reduction = "pca")^2
Importance<-vars/sum(vars)
print(Importance[1:10])

meta_categorical <- d10x.combined_myeloid@meta.data[, c("CellType_refined","CellType_harmonized","age_condition")]  # input column numbers in meta that contain categorical variables
meta_continuous <- d10x.combined_myeloid@meta.data[, c("percent.mt","nFeature_RNA","Age")]  # input column numbers in meta that contain continuous variables

save(embed,vars, Importance, meta_categorical, meta_continuous,Loadings, file=here("../../../projects/macparland/RE/PediatricAdult/processed_data","Fetal_ped_IFALD_adult_PCA_myeloid.RData"))


###########
# subset objects for local plotting
###########
#load(here("../../../projects/macparland/RE/PediatricAdult/processed_data","Fetal_IFALD_adult_ped_integrated.rds"))

load(here("../../../projects/macparland/RE/PediatricAdult/processed_dataFetal_IFALD_adult_ped_integrated.rds"))


d10x.fetal_ped_IFALD$CellType_harmonized<-d10x.fetal_ped_IFALD$CellType_refined
levels(d10x.fetal_ped_IFALD$CellType_harmonized)[which(levels(d10x.fetal_ped_IFALD$CellType_harmonized)=="Kupffer Cell")]<-"KC Like"
levels(d10x.fetal_ped_IFALD$CellType_harmonized)[which(levels(d10x.fetal_ped_IFALD$CellType_harmonized)%in%c("Mono-Mac","Monocyte-DC precursor","Monocyte" ))]<-"Mono-Mac"
levels(d10x.fetal_ped_IFALD$CellType_harmonized)[which(levels(d10x.fetal_ped_IFALD$CellType_harmonized)%in%c("DC2","DC1" ))]<-"Macrophage\n(MHCII high)"
levels(d10x.fetal_ped_IFALD$CellType_harmonized)[which(levels(d10x.fetal_ped_IFALD$CellType_harmonized)%in%c("NK-like cells","NK" ))]<-"NK cell"

d10x.combined_myeloid<-subset(d10x.fetal_ped_IFALD, subset = CellType_refined %in% c("RR Myeloid","KC Like","Macrophage\n(MHCII high)","Cycling Myeloid",
                                                                        "CDC1","VCAM1+ Erythroblastic Island Macrophage",
                                                                        "pDC precursor","Neutrophil-myeloid progenitor","Mono-NK",
                                                                        "Myeloid Erythrocytes\n(phagocytosis)","HSC/MPP","Monocyte-DC precursor","Mono-Mac",
                                                                        "Kupffer Cell","Neutrophil-myeloid progenitor","Mono-NK",
                                                                        "VCAM1+ Erythroblastic Island Macrophage","Monocyte","Erythroblastic Island Macrophage"))
d10x.combined_bcell<-subset(d10x.fetal_ped_IFALD, subset = CellType_harmonized %in% c("Mature B-cells","pro B cell","pre pro B cell","B cell",
                                                                      "pre B cell","Plasma cells"))
d10x.combined_tcell<-subset(d10x.fetal_ped_IFALD, subset = CellType_harmonized %in% c("NK cell","CD3+ T-cells","CLNK T-cells","Cycling T-cells",
                                                                      "ILC precursor","Early lymphoid/T lymphocyte","Mono-NK",
                                                                      "gd T-cells",""))
d10x.combined_HSC<-subset(d10x.fetal_ped_IFALD, subset = CellType_harmonized %in% c("HSC"))

save(d10x.combined_bcell, file=here("../../../projects/macparland/RE/PediatricAdult/processed_data/Fetal_IFALD_adult_ped_integrated_bcell_only.RData"))
save(d10x.combined_tcell, file=here("../../../projects/macparland/RE/PediatricAdult/processed_data/Fetal_IFALD_adult_ped_integrated_tcell_only.RData"))
save(d10x.combined_HSC, file=here("../../../projects/macparland/RE/PediatricAdult/processed_data/Fetal_IFALD_adult_ped_integrated_HSC_only.RData"))


########
## Myeloid overlapping
########

d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)
d10x.combined_myeloid <- RunUMAP(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindNeighbors(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindClusters(d10x.combined_myeloid, resolution = 0.2)


myeloid_cluster_umap<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=T)
myeloid_cluster_umap
save_plts(myeloid_cluster_umap, "IFALD_fetal_myeloid_cluster_umap", w=5,h=4)

myeloid_cluster_umap_individual<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=T, split.by = "age_condition", ncol=2)
myeloid_cluster_umap_individual
save_plts(myeloid_cluster_umap_individual, "IFALD_fetal_myeloid_cluster_umap_individual", w=10,h=6)

DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=T,
        group.by = "CellType_refined",split.by = "age_condition", ncol=2)+colscale_cellType_fetal

DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=T,
        group.by = "CellType_refined")+colscale_cellType_fetal


d10x.combined_myeloid$CellType_harmonized<-as.factor(d10x.combined_myeloid$CellType_refined)
levels(d10x.combined_myeloid$CellType_harmonized)[which(levels(d10x.combined_myeloid$CellType_harmonized)=="Kupffer Cell")]<-"KC Like"
# levels(d10x.combined_myeloid$CellType_harmonized)[which(levels(d10x.combined_myeloid$CellType_harmonized)%in%c("Mono-Mac"))]<-"Mono-Mac"
# levels(d10x.combined_myeloid$CellType_harmonized)[which(levels(d10x.combined_myeloid$CellType_harmonized)%in%c("Monocyte" ))]<-"Monocyte"


myeloid_cluster_umap<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=F,
                                         group.by = "CellType_harmonized",split.by = "age_condition", ncol=2)+colscale_cellType_fetal_combo
myeloid_cluster_umap
save_plts(myeloid_cluster_umap, "IFALD_fetal_myeloid_cluster_umap_groups", w=10,h=6)

######## Subset to just the overlapping cell types
d10x.combined_myeloid<-subset(d10x.combined_myeloid, subset = CellType_harmonized %in% c("Mono-Mac","KC Like","Macrophage\n(MHCII high)","CDC1","Monocyte"))
d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)
d10x.combined_myeloid <- RunUMAP(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindNeighbors(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindClusters(d10x.combined_myeloid, resolution = 0.2)

save(d10x.combined_myeloid, file=here("../../../projects/macparland/RE/PediatricAdult/processed_data/Fetal_IFALD_adult_ped_integrated_myeloid_only.RData"))

