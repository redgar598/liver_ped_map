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




source("R_functions/pretty_plots.R")


dataset_loc <- here("../../../projects/macparland/RE/PediatricAdult")

samples<-list.files(dataset_loc)
samples<-samples[-grep("meta",samples)]
print(samples)

#meta<-read.table(here("data/input_metadata.txt"), header=T)
meta<-read.table(here(dataset_loc,"input_metadata.txt"), header=T)
meta$Sample_ID[which(meta$Sample_ID=="C85_caud3pr")]<-"C85_caud5pr"


d10x.list <- sapply(1:length(samples), function(y){
  print(file.path(dataset_loc,paste(samples[y], sep=""),"filtered_feature_bc_matrix"))
  d10x <- Read10X(file.path(dataset_loc,paste(samples[y], sep=""),"filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),samples[y],sep="-")
  # print(dim(d10x))
  #' Initialize the Seurat object with the raw (non-normalized data).
  d10x<-CreateSeuratObject(counts = d10x, project = "ped_adult_map", min.cells = 0, min.features = 0)
  
  ## SoupX needs clusters so quickly make clusters for each sample
  d10x    <- SCTransform(d10x, verbose = F)
  d10x    <- RunPCA(d10x, verbose = F)
  d10x    <- RunUMAP(d10x, dims = 1:30, verbose = F)
  d10x    <- FindNeighbors(d10x, dims = 1:30, verbose = F)
  d10x    <- FindClusters(d10x, verbose = T)
  meta_clusters    <- d10x@meta.data
  
  sc = load10X(file.path(dataset_loc,paste(samples[y], sep="")))
  sc = setClusters(sc, setNames(meta_clusters$seurat_clusters, rownames(meta_clusters)))
  
  ######
  ## Load data and estimate soup profile
  ######
  # Estimate rho
  sc = autoEstCont(sc)
  #Genes with highest expression in background. These are often enriched for ribosomal proteins.
  print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20))
  # Clean the data
  sc = adjustCounts(sc)
  
  d10x = CreateSeuratObject(sc)
  
  #add meta data to each seurat object
  meta_cell<-data.frame(cell=colnames(d10x), individual=samples[y])
  meta_cell_add<-merge(meta_cell, meta, by.x="individual", by.y="Sample_ID")
  meta_cell_add<-meta_cell_add[match(colnames(d10x), meta_cell_add$cell),]
  print(identical(meta_cell_add$cell, colnames(d10x)))
  rownames(meta_cell_add)<-meta_cell_add$cell
  d10x<- AddMetaData(d10x, meta_cell_add)
  d10x
})

d10x.list



## cell counts
plt_count_raw<-lapply(1:length(d10x.list), function(x) {
  df<-data.frame(raw_cell_count=nrow(d10x.list[[x]]@meta.data),individual=unique(d10x.list[[x]]@meta.data$individual))
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
save(plt_QC_data, file=here("data","QC_metrics.Rdata"))

#load(here("data","QC_metrics.Rdata"))



qc_plts<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
  geom_point() +
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th
save_plts(qc_plts, "intital_QC_plts", w=6,h=4)

qc_plts_chem<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
  geom_point() + facet_wrap(~Chemistry)+
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th
save_plts(qc_plts_chem, "intital_QC_plts_chemistry", w=12,h=4)

qc_plts_individual<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
  geom_point() + facet_wrap(~individual)+
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th
save_plts(qc_plts_chem, "intital_QC_plts_individual", w=12,h=4)

MT_plt<-ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
  geom_vline(xintercept = 25)+ theme_bw()+xlab("Percent Mitochondrial")+th
save_plts(MT_plt, "percentMT_plt", w=6,h=4)

MT_plt_individual<-ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
  facet_wrap(~individual, scales="free_y")+
  geom_vline(xintercept = 25)+ theme_bw()+xlab("Percent Mitochondrial")+th
save_plts(MT_plt_individual, "percentMT_plt_individual", w=8,h=4)



#'We filter cells that have unique feature counts over 6,000 or less than 500
#'We filter cells that have >10% mitochondrial counts
#'we will also filter doublets as called by scrublet
d10x.list.raw<-d10x.list

invisible(lapply(1:length(d10x.list), function(x){
  d10x.list[[x]] <<- subset(d10x.list[[x]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 25)
}))

d10x.list

## cell counts after QC
plt_count_QC<-lapply(1:length(d10x.list), function(x) {
  df<-data.frame(qc_cell_count=nrow(d10x.list[[x]]@meta.data),individual=unique(d10x.list[[x]]@meta.data$individual))
  df})
plt_count_QC<-do.call(rbind, plt_count_QC)
print(plt_count_QC)

counts<-merge(plt_count_raw, plt_count_QC, by="individual")
meta<-merge(meta,counts,by.x="Sample_ID", by.y="individual")

cell_count<-grid.arrange(ggplot(meta, aes(AgeGroup, raw_cell_count,fill=AgeGroup))+
                           geom_boxplot()+geom_point()+
                           theme_bw()+geom_text(aes(label=Sample_ID), hjust=-0.25, size=3)+xlab("Age Group")+
                           ylab("Total Cell Number")+th+fillscale_age+ylim(0,60000)+
                           theme(legend.position = "none")+ggtitle("Before Quality Control"),
                         ggplot(meta, aes(AgeGroup, qc_cell_count,fill=AgeGroup))+
                           geom_boxplot()+geom_point()+
                           theme_bw()+geom_text(aes(label=Sample_ID), hjust=-0.25, size=3)+xlab("Age Group")+
                           ylab("Total Cell Number")+th+fillscale_age+ylim(0,60000)+
                           theme(legend.position = "none")+ggtitle("After Quality Control"), ncol=2)

save_plts(cell_count, "QC_cellcount_age", w=8,h=4)


d10x <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "adult_ped_map")#add.cell.ids = alldata_names2,

d10x

saveRDS(d10x, file = here("data","d10x_adult_ped_raw.rds"))


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
save_plts(pca_cellcycle, "pca_cellcycle", w=6,h=4)

pca_nfeature<-FeaturePlot(d10x, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
save_plts(pca_nfeature, "pca_nfeature", w=6,h=4)



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
save_plts(SCT_cluster_umap, "SCT_cluster_umap", w=6,h=4)

SCT_cluster_tsne<-DimPlot(d10x, reduction = "tsne", pt.size=0.25, label=T)
save_plts(SCT_cluster_tsne, "SCT_cluster_tsne", w=6,h=4)


cell_pca_SCT<-DimPlot(d10x, reduction="pca", group.by="Phase")
save_plts(cell_pca_SCT, "cell_PCA_afterSCT", w=6,h=4)

nFeature_UMAP_SCT<-FeaturePlot(d10x, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
save_plts(nFeature_UMAP_SCT, "nfeature_UMAP_afterSCT", w=6,h=4)



chem_umap_sct<-DimPlot(d10x, reduction = "umap", group.by = "Chemistry", pt.size=0.25)
save_plts(chem_umap_sct, "chem_SCT_umap", w=6,h=4)

age_umap_sct<-DimPlot(d10x, reduction = "umap", group.by = "AgeGroup", pt.size=0.25)+fillscale_age
save_plts(age_umap_sct, "age_SCT_umap", w=6,h=4)

individual_umap_sct<-DimPlot(d10x, reduction = "umap", group.by = "individual", pt.size=1)
save_plts(individual_umap_sct, "individual_SCT_UMAP", w=6,h=4)


## Low quality cluster?
MT_umap_SCT<-FeaturePlot(d10x, features = "percent.mt", min.cutoff = "q9", pt.size=1)
save_plts(pca_nfeature, "pca_nfeature", w=6,h=4)

ncount_umap_SCT<-FeaturePlot(d10x, features = "nCount_RNA", min.cutoff = "q9", pt.size=1)
save_plts(ncount_umap_SCT, "ncount_umap_SCT", w=6,h=4)

nfeature_umap_SCT<-FeaturePlot(d10x, features = "nFeature_RNA", min.cutoff = "q9", pt.size=1)
save_plts(nfeature_umap_SCT, "nfeature_umap_SCT", w=6,h=4)





###############
## Integration
###############
#https://satijalab.org/seurat/articles/integration_rpca.html
print("RUNNING INTEGRATION")
## Start back with raw data
d10x <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "adult_ped_map")#add.cell.ids = alldata_names2,
d10x

# split the dataset into a list of two seurat objects (3' and 5')
d10x.list.chem <- SplitObject(d10x, split.by = "Chemistry")

# normalize, identify variable features and score cell cycle for each dataset independently
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

d10x.list.chem <- lapply(X = d10x.list.chem, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = d10x.list.chem)
d10x.list.chem <- lapply(X = d10x.list.chem, FUN = function(x) {
  #x <- ScaleData(x, features = features, verbose = FALSE)
  x <- ScaleData(x, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})



## Identify anchors
chem.anchors <- FindIntegrationAnchors(object.list = d10x.list.chem, anchor.features = features, reduction = "rpca")
d10x.combined <- IntegrateData(anchorset = chem.anchors)

DefaultAssay(d10x.combined) <- "integrated"

print("INTEGRATED")

# Run the standard workflow for visualization and clustering
d10x.combined <- ScaleData(d10x.combined, verbose = FALSE)
d10x.combined <- RunPCA(d10x.combined, npcs = 30, verbose = FALSE)
d10x.combined <- RunUMAP(d10x.combined, reduction = "pca", dims = 1:30)
d10x.combined <- RunTSNE(d10x.combined, dims = 1:30)

d10x.combined <- FindNeighbors(d10x.combined, reduction = "pca", dims = 1:30)
d10x.combined <- FindClusters(d10x.combined, resolution = 0.5)

d10x.combined


###########
## Visualize integration
###########
SCT_cluster_umap<-DimPlot(d10x.combined, reduction = "umap", pt.size=0.25, label=T)
save_plts(SCT_cluster_umap, "rPCA_cluster_umap", w=6,h=4)

SCT_cluster_tsne<-DimPlot(d10x.combined, reduction = "tsne", pt.size=0.25, label=T)
save_plts(SCT_cluster_tsne, "rPCA_cluster_tsne", w=6,h=4)

chem_umap_sct<-DimPlot(d10x.combined, reduction = "umap", group.by = "Chemistry", pt.size=0.25)
save_plts(chem_umap_sct, "chem_rPCA_umap", w=6,h=4)

age_umap_sct<-DimPlot(d10x.combined, reduction = "umap", group.by = "AgeGroup", pt.size=0.25)+fillscale_age
save_plts(age_umap_sct, "age_rPCA_umap", w=6,h=4)

MT_umap_sct<-FeaturePlot(d10x.combined, reduction = "umap", features = "percent.mt", pt.size=0.25)
save_plts(MT_umap_sct, "MT_rPCA_umap", w=5,h=4)

individual_umap_sct<-DimPlot(d10x.combined, reduction = "umap", group.by = "individual", pt.size=0.25)
save_plts(individual_umap_sct, "individual_rPCA_UMAP", w=6,h=4)

individual_split<-DimPlot(d10x.combined, reduction = "umap", group.by = "individual", split.by="individual",pt.size=0.25)
save_plts(individual_split, "individual_facet_rPCA_UMAP", w=20,h=4)

age_split<-DimPlot(d10x.combined, reduction = "umap", group.by = "individual", split.by="AgeGroup",pt.size=0.25)
save_plts(age_split, "age_facet_rPCA_UMAP", w=10,h=5)


############################
##### Example markers to plot for the liver [From Diana]
############################

######
#Immune cells
######

## Macrophages
DotPlot(d10x.combined, features = c( "PTPRC", "CD68", "MARCO","CD5L","VSIG4", "MAF", "LYZ", "CSTA", "S100A8", "S10049", "CD14", "CD74", "GPBAR1", "ID3"), cols=c("blue", "red")) + RotatedAxis()

## NK/T/B cells
DotPlot(d10x.combined, features = c("PTPRC", "CD2", "CD3E", "IL7R", "KLRB1",
                                     "NKG7", "GZMA", "GZMB", "GZMK" , "PRF1","CD4", "CD8A","CD247", "TRAC","TRDC", "TRGC1", "TRGC2", "TRBC1",
                                     "TRBC2", "S1PR1", "CD28", "CD27", "SELL", "CCR7", "CXCR4","CCR4","FAS",
                                     "FOXP3", "CTLA4", "LAG3", "TNFRSF4","TNFRSF18", "ICOS" ,"CD69", "CD79A", "CD79B", "IGHG1", "MS4A1",
                                     "LTB", "CD52", "IGHD", "CD19", "ID3"), cols=c("blue", "red")) + RotatedAxis()

# LEC and LSEC
DotPlot(d10x.combined, features = c("CALCRL", "VWF", "RAMP2", "STAB2", "LYVE1", "PECAM1", "ENG", "FCGR2B", "F8", "SPARCL1", "ID1", "SOX18", "CD32B", "ID3"), cols=c("blue", "red")) + RotatedAxis()

######
#Hepatocytes
######
DotPlot(d10x.combined, features=c("ALB", "HAMP", "ARG1", "PCK1", "AFP", "BCHE", "HAL", "SCD", "CPS1", "CYP3A4",
                                   "ELF3", "CRP", "GSTA2", "AKR1C1", "MGST1", "CYP3A5", "ALDH1A1", "ADH1A", "CYP2E1",
                                   "GLS2", "SDS", "GLUL", "AKR1D1", "HPR",
                                   "HMGCS1", "IGSF23", "ACSS2", "G6PC", "ID3"),
        cols=c("blue", "red")) + RotatedAxis()



######
#Cholangiocytes
######
DotPlot(d10x.combined,features = c( "EPCAM", "SOX9", "KRT1", "KRT7", "ANXA4", "KRT18", "ID3"), cols=c("blue", "red")) + RotatedAxis()

######
#HSCs
######
DotPlot(d10x.combined,features = c( "RBP1", "LRAT", "PDE3B", "ACTA2", "AOX1", "PDE3D", "PDE4D", "SPARC", "TAGLN", "COL1A1", "COL1A2", "COL3A1", "TIMP1", "DCN", "MYL9", "TPM2", "MEG3", "BGN", "IGFBP7", "IGFBP3", "CYR61", "IGFBP6", "CCL2", "COLEC11", "CTGF", "HGF", "ID3"), cols=c("blue", "red")) + RotatedAxis()
FeaturePlot(d10x.combined, reduction = "umap", features = c("PTPRC", "CD3D", "CD68", "CD79A","TRDC", "NKG7", "KRT7", "CALCRL", "ACTA2", "MS4A1", "CYP3A4", "SCD", "FCN2", "CD4", "CD8A", "FCER1A", "MARCO", "LYZ", "VSIG4", "FOLR2", "ID3"), ncol = 4)
key_markers<-FeaturePlot(d10x.combined, reduction = "umap", features = c("EPCAM", "SOX9", "SELL", "PTPRC",
                                                            "TRDC", "NKG7", "CALCRL", "VWF",
                                                            "MARCO", "LYZ","COL1A1","IGFBP3"), ncol = 4)
save_plts(key_markers, "markers_rPCA_UMAP", w=25,h=20)



#################
## Rough annotation
#################
Macrophage_genes<-c( "PTPRC", "CD68", "MARCO","CD5L","VSIG4", "MAF", "LYZ", "CSTA", "S100A8", "S10049",
                     "CD14", "CD74", "GPBAR1", "ID3")
NK_T_B_genes<-c("PTPRC", "CD2", "CD3E", "IL7R", "KLRB1","NKG7", "GZMA", "GZMB", "GZMK" , "PRF1","CD4",
                "CD8A","CD247", "TRAC","TRDC", "TRGC1", "TRGC2", "TRBC1","TRBC2", "S1PR1", "CD28", "CD27",
                "SELL", "CCR7", "CXCR4","CCR4","FAS",  "FOXP3", "CTLA4", "LAG3", "TNFRSF4","TNFRSF18",
                "ICOS" ,"CD69", "CD79A", "CD79B", "IGHG1", "MS4A1","LTB", "CD52", "IGHD", "CD19", "ID3")

B_genes<-c("POU2F2","FCER2","MS4A1")
T_genes<-c("CD3D","TRDC","SELL","IL7R","CCR7","S100A4","CD8A")
NK_genes<-c("GNLY","NKG7")


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

genes<-unique(c(Macrophage_genes, NK_T_B_genes,LEC_genes,Hepatocyte_genes,Cholangiocytes_genes,HSCs_genes))

d10x.exp<-as.data.frame(d10x.combined[["RNA"]]@data)
d10x.exp.GOI<-d10x.exp[genes,]


# GENES<-rownames(d10x.combined)
# indMyGene<-which(GENES%in%genes)
#
# CountMyGeneMyCell<-GetAssayData(object = d10x.combined, slot = 'counts', assay='RNA')[which(rownames(d10x.combined)%in%genes), ]
# d10x.exp.GOI<-as.data.frame(as.matrix(CountMyGeneMyCell))

d10x.exp.GOI$gene<-rownames(d10x.exp.GOI)
d10x.exp.GOI<-melt(d10x.exp.GOI)#

meta<-d10x.combined@meta.data
meta$cell<-rownames(meta)

plt<-merge(d10x.exp.GOI, meta,by.x="variable", by.y="cell")

plt$variable<-as.character(plt$variable)

## possible cells types
cluster_marker_mean<-function(gene_list, type){
  plt_epi<-plt[which(plt$gene%in%gene_list),]
  mean_type<-as.data.frame(tapply(plt_epi$value, plt_epi$seurat_clusters, mean))
  colnames(mean_type)<-type
  mean_type
}

cell_rough<-cbind(cluster_marker_mean(Macrophage_genes, "Myeloid"),
      cluster_marker_mean(NK_T_B_genes, "NK_T_B"),
      cluster_marker_mean(LEC_genes, "LEC"),
      cluster_marker_mean(Hepatocyte_genes, "Hepatocyte"),
      cluster_marker_mean(Cholangiocytes_genes, "Cholangiocytes"),
      cluster_marker_mean(HSCs_genes, "HSC"))


cell_rough$CellType_rough<-sapply(1:nrow(cell_rough), function(x) {
  compart<-colnames(cell_rough)[which(cell_rough[x,] == max(cell_rough[x,]))]
  if(length(compart)==1){compart}else{"Unclear"}
})

cell_rough$seurat_clusters<-rownames(cell_rough)

meta<-d10x.combined@meta.data
meta$cell<-rownames(meta)
plt_summary<-merge(meta, cell_rough[,c("seurat_clusters","CellType_rough")], by="seurat_clusters")
plt_summary<-plt_summary[match(rownames(d10x.combined@meta.data),plt_summary$cell),]
identical(plt_summary$cell, rownames(d10x.combined@meta.data))

rownames(plt_summary)<-plt_summary$cell

d10x.combined<- AddMetaData(d10x.combined, plt_summary)

###########
## Sub cluster NK T B
###########

B_genes<-c("POU2F2","FCER2","MS4A1","LTB","CD37","CD79B","IGLC2","IGHG1","IGKC", "CD19")
T_genes<-c("CD3D","IL7R","CD8A","IL32")
NK_genes<-c("NKG7","CD7")
gd_genes<-c("GNLY")

d10x.combined_NK_T_B<-subset(d10x.combined, subset = CellType_rough %in% c("NK_T_B"))
d10x.combined_NK_T_B <- RunPCA(d10x.combined_NK_T_B, npcs = 30, verbose = FALSE)
d10x.combined_NK_T_B <- RunUMAP(d10x.combined_NK_T_B, reduction = "pca", dims = 1:30)

DimPlot(d10x.combined_NK_T_B, reduction = "umap", pt.size=0.25)

FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = NK_genes, ncol = 2)
FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = gd_genes)
FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = T_genes, ncol = 2)
FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = B_genes, ncol = 2)

## same method on sub cluster
genes<-unique(c(T_genes, NK_genes,gd_genes))

d10x.exp<-as.data.frame(d10x.combined_NK_T_B[["RNA"]]@data)
d10x.exp.GOI<-d10x.exp[genes,]
d10x.exp.GOI$gene<-rownames(d10x.exp.GOI)
d10x.exp.GOI<-melt(d10x.exp.GOI)#

meta<-d10x.combined_NK_T_B@meta.data
meta$cell<-rownames(meta)

plt<-merge(d10x.exp.GOI, meta,by.x="variable", by.y="cell")
plt$variable<-as.character(plt$variable)
plt$seurat_clusters<-as.character(plt$seurat_clusters)


## possible cells types
cluster_marker_mean<-function(gene_list, type){
  plt_epi<-plt[which(plt$gene%in%gene_list),]
  mean_type<-as.data.frame(tapply(plt_epi$value, plt_epi$seurat_clusters, mean))
  colnames(mean_type)<-type
  mean_type
}

cell_rough<-cbind(cluster_marker_mean(T_genes, "CD3_Tcell"),
                  cluster_marker_mean(NK_genes, "nkTcell"),
                  cluster_marker_mean(gd_genes, "gdTcell"))


cell_rough$CellType_rough_sub<-sapply(1:nrow(cell_rough), function(x) {
  compart<-colnames(cell_rough)[which(cell_rough[x,] == max(cell_rough[x,]))]
  if(length(compart)==1){compart}else{"Unclear"}
})

cell_rough$seurat_clusters<-rownames(cell_rough)

lapply(1:nrow(cell_rough),function(x){
  d10x.combined@meta.data[which(d10x.combined@meta.data$seurat_clusters==cell_rough$seurat_clusters[x]),]$CellType_rough<<-cell_rough$CellType_rough_sub[x]})
table(d10x.combined@meta.data$CellType_rough)


##############
## Save integrated to look local
##############
save(d10x.combined, file=paste(here("data/"),"adult_ped_integrated.rds", sep=""))
cell_label<-d10x.combined@meta.data
save(cell_label, file=paste(here("data/"),"adult_ped_cellRough.rds", sep=""))


load(here("data","adult_ped_integrated.rds"))


#####
## plot cell types
#####
d10x.combined@meta.data$CellType_rough<-as.factor(d10x.combined@meta.data$CellType_rough)
levels(d10x.combined@meta.data$CellType_rough)<-c("CD3+ T-cells","Cholangiocytes",
                                                  "gd T-cells","Hepatocytes",
                                                  "HSC","LSEC","Myeloid cells","NK-like cells")

roughcell_cluster_umap<-DimPlot(d10x.combined, reduction = "umap",group.by="CellType_rough", pt.size=0.15, label=T)+colscale_cellType+ggtitle("")+
  annotate("text", x = -11, y = -12, label = paste0("n = ",comma(ncol(d10x.combined))))
roughcell_cluster_umap
save_plts(roughcell_cluster_umap, "rPCA_roughcellType_cluster_umap", w=6,h=4)
roughcell_cluster_umap<-DimPlot(d10x.combined, reduction = "umap",group.by="CellType_rough", pt.size=0.15)+colscale_cellType+ggtitle("")+
  annotate("text", x = -11, y = -12, label = paste0("n = ",comma(ncol(d10x.combined))))
roughcell_cluster_umap
save_plts(roughcell_cluster_umap, "rPCA_roughcellType_cluster_umap_nolab", w=6,h=4)

d10x.combined$age_id<-paste(d10x.combined$individual, d10x.combined$AgeGroup)
d10x.combined$age_id[which(d10x.combined$age_id=="C85_caud5pr Ped")]<-"C85_caud5pr Ped (Frozen)"
d10x.combined$age_id<-as.factor(d10x.combined$age_id)
#order by % cell passing QC
d10x.combined$age_id<-factor(d10x.combined$age_id, c("C86_caud3pr Adult","C92_caud3pr Adult", "C63_caud5pr Adult","C70_caud5pr Adult",
                                                     "C61_caud5pr Adult","C82_caud3pr Adult","C98_caud3pr Adult",
                                                     "C85_caud5pr Ped (Frozen)","C96_caud3pr Ped","C93_caud3pr Ped", "C64_caud5pr Ped"))
levels(d10x.combined$age_id)<-gsub(" ","\n", levels(d10x.combined$age_id))

individual_split<-DimPlot(d10x.combined, reduction = "umap", group.by = "CellType_rough", split.by="age_id",pt.size=0.25, ncol=5)+colscale_cellType+ggtitle("")
save_plts(individual_split, "individual_roughCell_facet_rPCA_UMAP", w=22,h=8)

cell_num_all<-as.data.frame(table(d10x.combined@meta.data$AgeGroup))
colnames(cell_num_all)<-c("AgeGroup","CellCount")
age_split<-DimPlot(d10x.combined, reduction = "umap", group.by = "CellType_rough", split.by="AgeGroup",pt.size=0.25)+colscale_cellType+ggtitle("")+
  geom_text(aes(x=-11, y=-12, label=paste0("n = ",comma(CellCount))), cell_num_all)
save_plts(age_split, "age_roughCell_facet_rPCA_UMAP", w=10,h=5)



########
## mean genes per cell
########
d10x.combined@meta.data %>%
  group_by(individual) %>%
  summarise(max = mean(nFeature_RNA, na.rm=TRUE))


######
#TRM interesting markers
######

#klrg1/cd57
#cd69/cd103

#CD57 is encoded by  B3GAT1
#CD103 is encoded by ITGAE

TRM_markers_all<-FeaturePlot(d10x.combined, reduction = "umap", features = c("KLRG1", "B3GAT1", "CD69","ITGAE"), ncol = 2)
save_plts(TRM_markers_all, "TRM_markers", w=5,h=4)

d10x.combined_NK_T_B<-subset(d10x.combined, subset = CellType_rough %in% c("CD3+ T-cells","gd T-cells","NK-like cells"))
d10x.combined_NK_T_B <- RunPCA(d10x.combined_NK_T_B, npcs = 30, verbose = FALSE)
d10x.combined_NK_T_B <- RunUMAP(d10x.combined_NK_T_B, reduction = "pca", dims = 1:30)

FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = c("KLRG1", "B3GAT1", "CD69","ITGAE"), ncol = 2)
TRM_markers_NKTB<-FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = c("KLRG1", "B3GAT1", "CD69","ITGAE"), split="AgeGroup", ncol = 2)
save_plts(TRM_markers_NKTB, "TRM_markers_NK_T_B", w=5,h=10)

d10x.exp<-as.data.frame(d10x.combined_NK_T_B[["RNA"]]@data)
d10x.exp.GOI<-d10x.exp[c("KLRG1", "B3GAT1", "CD69","ITGAE"),]

d10x.exp.GOI$gene<-rownames(d10x.exp.GOI)
d10x.exp.GOI<-melt(d10x.exp.GOI)#

meta<-d10x.combined_NK_T_B@meta.data
meta$cell<-rownames(meta)

plt<-merge(d10x.exp.GOI, meta,by.x="variable", by.y="cell")

quant_markers<-ggplot(plt, aes(AgeGroup, value))+geom_violin(aes(fill=AgeGroup))+facet_wrap(~gene)+fillscale_age+theme_bw()+th
save_plts(quant_markers, "TRM_markers_NK_T_B_quantified", w=5,h=4)


##############
## Really no hepatocytes in ped samples
##############
## diff in ped (none in ped?)
hep_age<-FeaturePlot(d10x.combined, reduction = "umap", features = c("ALB", "CPS1", "CYP3A4",
                                                            "MGST1", "CYP2E1"),
            ncol = 2, split.by='AgeGroup')
save_plts(hep_age, "hep_markers_age_rPCA_UMAP", w=6,h=14)

##############
## Myeloid clustering
##############
d10x.combined_myeloid<-subset(d10x.combined, subset = CellType_rough %in% c("Myeloid cells"))
d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)
d10x.combined_myeloid <- RunUMAP(d10x.combined_myeloid, reduction = "pca", dims = 1:30)

d10x.combined_myeloid <- FindNeighbors(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindClusters(d10x.combined_myeloid, resolution = 0.1)


myeloid_cluster_umap<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=T)
myeloid_cluster_umap
save_plts(myeloid_cluster_umap, "myeloid_cluster_umap", w=5,h=4)

myeloid_cluster_umap_individual<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=T, split.by = "age_id", ncol=5)
myeloid_cluster_umap_individual
save_plts(myeloid_cluster_umap_individual, "myeloid_cluster_umap_individual", w=10,h=6)


cell_cluster_count<-d10x.combined_myeloid@meta.data %>%  group_by(age_id, seurat_clusters) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

cell_cluster_count<-as.data.frame(cell_cluster_count)
cell_cluster_count<-merge(cell_cluster_count, d10x.combined_myeloid@meta.data[,c("age_id","AgeGroup")], by="age_id")
cell_cluster_count<-cell_cluster_count[!duplicated(cell_cluster_count),]

bar_individual<-ggplot(cell_cluster_count, aes(fill=seurat_clusters, y=n, x=age_id)) + 
  geom_bar(position="fill", stat="identity", color="black")+theme_bw()+th+
  facet_wrap(~AgeGroup, scale="free_x")
save_plts(bar_individual, "bar_individual_myeloid", w=14,h=6)

## low quality clusters?
umap_MTmyeloid<-FeaturePlot(d10x.combined_myeloid, features = "percent.mt", min.cutoff = "q9", pt.size=1)
save_plts(umap_MTmyeloid, "umap_MTmyeloid", w=5,h=4)

cluster_MT<-d10x.combined_myeloid@meta.data %>%  group_by(seurat_clusters) %>%
  dplyr::summarize(Mean = mean(percent.mt, na.rm=TRUE))
cell_cluster_count<-as.data.frame(cluster_MT)

box_MT<-ggplot(d10x.combined_myeloid@meta.data, aes(seurat_clusters, percent.mt)) + 
  geom_violin(fill="lightgrey", color="lightgrey")+xlab("Myeloid Cluster")+ylab("Percent MT")+
  geom_boxplot(width=0.1, outlier.size = 0.05)+theme_bw()+th
box_MT
save_plts(box_MT, "box_MT", w=8,h=3)

umap_nfeaturemyeloid<-FeaturePlot(d10x.combined_myeloid, features = "nFeature_RNA", min.cutoff = "q9", pt.size=1)
save_plts(umap_nfeaturemyeloid, "umap_nfeaturemyeloid", w=5,h=4)

cluster_nFeature<-d10x.combined_myeloid@meta.data %>%  group_by(seurat_clusters) %>%
  dplyr::summarize(Mean = mean(nFeature_RNA, na.rm=TRUE))
cell_cluster_count<-as.data.frame(cluster_nFeature)

box_nfeature<-ggplot(d10x.combined_myeloid@meta.data, aes(seurat_clusters, nFeature_RNA)) + 
  geom_violin(fill="lightgrey", color="lightgrey")+xlab("Myeloid Cluster")+ylab("Number of Feature \n(nFeature_RNA)")+
  geom_boxplot(width=0.1, outlier.size = 0.05)+theme_bw()+th
box_nfeature
save_plts(box_nfeature, "box_nfeature", w=8,h=3)

print(sessionInfo())




