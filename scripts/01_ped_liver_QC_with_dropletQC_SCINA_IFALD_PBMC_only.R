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

#' 
#' dataset_loc <- here("/media/redgar/Seagate Portable Drive/ped_liver_map_raw")
#' #dataset_loc <- here("../../../projects/macparland/RE/PediatricAdult/ped_liver_map_raw")
#' 
#' samples<-list.files(dataset_loc)
#' #samples<-samples[-grep("data_transfer",samples)]
#' print(samples)
#' 
#' meta<-read.table(here("data/data_transfer_updated_may15_2023_IFALD_PBMC.csv"), header=T, sep=",")
#' #meta<-read.table(here(dataset_loc,"data_transfer_updated_may15_2023_IFALD_PBMC.csv"), header=T, sep=",")
#' 
#' samples<-samples[which(samples%in%meta$file)]
#' 
#' samples<-samples[grep("PBMC",samples)]
#' 
#' d10x.list <- sapply(1:length(samples), function(y){
#'   caud<-meta$Sample_ID[which(meta$file == samples[y])]
#'   print(caud)
#'   print(file.path(dataset_loc,paste(samples[y],"/outs", sep=""),"filtered_feature_bc_matrix"))
#'   d10x <- Read10X(file.path(dataset_loc,paste(samples[y],"/outs", sep=""),"filtered_feature_bc_matrix"))
#'   colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),caud,sep="-")
#'   # print(dim(d10x))
#'   #' Initialize the Seurat object with the raw (non-normalized data).
#'   d10x<-CreateSeuratObject(counts = d10x, project = "ped_adult_map", min.cells = 0, min.features = 0)
#'   
#'   if(dir.exists(file.path(dataset_loc,paste(samples[y],"/outs/raw_feature_bc_matrix", sep="")))) {
#'     ## SoupX needs clusters so quickly make clusters for each sample
#'     d10x    <- SCTransform(d10x, verbose = F)
#'     d10x    <- RunPCA(d10x, verbose = F)
#'     d10x    <- RunUMAP(d10x, dims = 1:30, verbose = F)
#'     d10x    <- FindNeighbors(d10x, dims = 1:30, verbose = F)
#'     d10x    <- FindClusters(d10x, verbose = T)
#'     meta_clusters    <- d10x@meta.data
#'     
#'     sc = load10X(file.path(dataset_loc,paste(samples[y],"/outs", sep="")))
#'     sc = setClusters(sc, setNames(meta_clusters$seurat_clusters, rownames(meta_clusters)))
#'     
#'     ######
#'     ## Load data and estimate soup profile
#'     ######
#'     # Estimate rho
#'     sc = autoEstCont(sc, forceAccept=TRUE)
#'     #Genes with highest expression in background. These are often enriched for ribosomal proteins.
#'     print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20))
#'     print(unique(sc$metaData$rho))
#'     # Clean the data
#'     out = adjustCounts(sc)
#'     
#'     ## Save a metric of soupness (ALB change after soupX)
#'     DR = sc$metaData[,sc$DR]
#'     df = DR
#'     old = colSums(sc$toc["ALB",rownames(df),drop=FALSE])
#'     new = colSums(out["ALB",rownames(df),drop=FALSE])
#'     relChange = (old-new)/old
#'     df$old = old
#'     df$new = new
#'     df$relALBChange=relChange
#'     df$cell<-rownames(df)
#'     
#'     ## make seurat object of adjusted counts
#'     d10x = CreateSeuratObject(out)
#'     
#'     #add meta data to each seurat object
#'     meta_cell<-data.frame(cell=colnames(d10x), individual=caud)
#'     meta_cell_add<-merge(meta_cell, meta, by.x="individual", by.y="Sample_ID")
#'     meta_cell_add<-merge(meta_cell_add, df[,c("cell","relALBChange")], by="cell")
#'     meta_cell_add<-meta_cell_add[match(colnames(d10x), meta_cell_add$cell),]
#'     print(identical(meta_cell_add$cell, colnames(d10x)))
#'     rownames(meta_cell_add)<-meta_cell_add$cell
#'     d10x<- AddMetaData(d10x, meta_cell_add)
#'   }else{
#'     #add meta data to each seurat object
#'     meta_cell<-data.frame(cell=colnames(d10x), individual=caud)
#'     meta_cell_add<-merge(meta_cell, meta, by.x="individual", by.y="Sample_ID")
#'     meta_cell_add$relALBChange<-NA
#'     meta_cell_add<-meta_cell_add[match(colnames(d10x), meta_cell_add$cell),]
#'     print(identical(meta_cell_add$cell, colnames(d10x)))
#'     rownames(meta_cell_add)<-meta_cell_add$cell
#'     d10x<- AddMetaData(d10x, meta_cell_add)
#'   }
#'   
#'   ######
#'   ## dropletQC
#'   ######
#'   if(file.exists(file.path(dataset_loc,paste(samples[y],"/outs", sep=""),"possorted_genome_bam.bam"))){
#'     nf1 <- nuclear_fraction_tags(
#'       outs = file.path(dataset_loc,paste(samples[y],"/outs", sep="")),
#'       tiles = 1, cores = 1, verbose = FALSE)
#'     head(nf1)
#'     
#'     print(identical(rownames(nf1), colnames(d10x)))
#'     d10x<- AddMetaData(d10x, nf1)
#'     d10x
#'     
#'     nf.umi <- data.frame(nf=d10x$nuclear_fraction,
#'                          umi=d10x$nCount_RNA)
#'     
#'     # Run identify_empty_drops
#'     empty_drop <- identify_empty_drops(nf_umi=nf.umi)
#'     empty_drop$individual<-d10x$individual
#'     empty_drop_damagedcell <- identify_damaged_cells(empty_drop, verbose = FALSE, output_plots = F)
#'     
#'     head(empty_drop_damagedcell[[1]])
#'     table(empty_drop_damagedcell[[1]]$cell_status)
#'     
#'     print(identical(rownames(empty_drop_damagedcell[[1]]), colnames(d10x)))
#'     d10x<- AddMetaData(d10x, empty_drop_damagedcell[[1]])
#'     d10x$nf<-NULL
#'     d10x$umi<-NULL
#'     print(d10x)
#'     d10x}else{
#'       d10x$nuclear_fraction<-NA
#'       d10x$cell_status<-NA
#'       print(d10x)
#'       d10x}
#' })
#' 
#' d10x.list
#' 
#' 
#' 
#' ## cell counts
#' plt_count_raw<-lapply(1:length(d10x.list), function(x) {
#'   df<-data.frame(raw_cell_count=nrow(d10x.list[[x]]@meta.data),individual=unique(d10x.list[[x]]@meta.data$individual))
#'   df})
#' plt_count_raw<-do.call(rbind, plt_count_raw)
#' print(plt_count_raw)
#' 
#' 
#' 
#' #'## QC
#' #'The percentage of reads that map to the mitochondrial genome
#' #'Low-quality / dying cells often exhibit extensive mitochondrial contamination
#' #'We calculate mitochondrial QC metrics with the PercentageFeatureSet function, which calculates the percentage of counts originating from a set of features
#' #'We use the set of all genes starting with MT- as a set of mitochondrial genes
#' 
#' # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#' invisible(lapply(1:length(d10x.list), function(x){
#'   d10x.list[[x]][["percent.mt"]] <<- PercentageFeatureSet(d10x.list[[x]], pattern = "^MT-")}))
#' 
#' 
#' 
#' #'Low-quality cells or empty droplets will often have very few genes
#' #'Cell doublets or multiplets may exhibit an aberrantly high gene count
#' 
#' 
#' # Visualize QC metrics
#' #nFeature number of unique genes
#' #nCount number of total molecules
#' plt_QC_data<-do.call(rbind, lapply(1:length(d10x.list), function(x) d10x.list[[x]]@meta.data))
#' save(plt_QC_data, file=here("data","IFALD_QC_metrics_PBMC.Rdata"))
#' 
#' #load(here("data","IFALD_QC_metrics.Rdata"))
#' 
#' 
#' 
#' qc_plts<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
#'   geom_point() +
#'   scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
#'   geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
#'   geom_hline(yintercept = 6000) +theme_bw()+th
#' save_plts(qc_plts, "IFALD_intital_QC_plts_PBMC", w=6,h=4)
#' 
#' qc_plts_chem<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
#'   geom_point() + facet_wrap(~chemistry)+
#'   scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
#'   geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
#'   geom_hline(yintercept = 6000) +theme_bw()+th
#' save_plts(qc_plts_chem, "IFALD_intital_QC_plts_chemistry_PBMC", w=12,h=4)
#' 
#' qc_plts_individual<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
#'   geom_point() + facet_wrap(~individual)+
#'   scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
#'   geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
#'   geom_hline(yintercept = 6000) +theme_bw()+th
#' save_plts(qc_plts_chem, "IFALD_intital_QC_plts_individual_PBMC", w=12,h=4)
#' 
#' MT_plt<-ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
#'   geom_vline(xintercept = 25)+ theme_bw()+xlab("Percent Mitochondrial")+th
#' save_plts(MT_plt, "IFALD_percentMT_plt_PBMC", w=6,h=4)
#' 
#' MT_plt_individual<-ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
#'   facet_wrap(~individual, scales="free_y")+
#'   geom_vline(xintercept = 25)+ theme_bw()+xlab("Percent Mitochondrial")+th
#' save_plts(MT_plt_individual, "IFALD_percentMT_plt_individual_PBMC", w=8,h=4)
#' 
#' 
#' 
#' nf_plts<-ggplot(plt_QC_data, aes(nuclear_fraction,log10(nCount_RNA),colour=cell_status)) +
#'   geom_point() +  ylab("Number of Total Molecules\n(log 10 nCount) ")+xlab("Nuclear Fraction")+theme_bw()+th+
#'   facet_wrap(~individual)
#' nf_plts
#' save_plts(nf_plts, "IFALD_nuclear_fraction_PBMC", w=6,h=4)
#' 
#' ggplot(plt_QC_data, aes(nuclear_fraction,percent.mt,colour=cell_status)) +
#'   geom_point() +  ylab("Number of Total Molecules\n(log 10 nCount) ")+xlab("Nuclear Fraction")+theme_bw()+th+
#'   facet_wrap(~individual)+
#'   scale_color_manual(values=c("grey","red","cornflowerblue"),name="Nuclear\nFraction")
#' 
#' 
#' #'We filter cells that have unique feature counts over 6,000 or less than 500
#' #'We filter cells that have >10% mitochondrial counts
#' #'we will also filter doublets as called by scrublet
#' d10x.list.raw<-d10x.list
#' 
#' invisible(lapply(1:length(d10x.list), function(x){
#'   d10x.list[[x]] <<- subset(d10x.list[[x]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 25)
#' }))
#' 
#' d10x.list
#' 
#' ## cell counts after QC
#' plt_count_QC<-lapply(1:length(d10x.list), function(x) {
#'   df<-data.frame(qc_cell_count=nrow(d10x.list[[x]]@meta.data),individual=unique(d10x.list[[x]]@meta.data$individual))
#'   df})
#' plt_count_QC<-do.call(rbind, plt_count_QC)
#' print(plt_count_QC)
#' 
#' counts<-merge(plt_count_raw, plt_count_QC, by="individual")
#' meta<-merge(meta,counts,by.x="Sample_ID", by.y="individual")
#' 
#' cell_count<-grid.arrange(ggplot(meta, aes(AgeGroup, raw_cell_count,fill=AgeGroup))+
#'                            geom_boxplot()+geom_point()+
#'                            theme_bw()+geom_text(aes(label=Sample_ID), hjust=-0.25, size=3)+xlab("Age Group")+
#'                            ylab("Total Cell Number")+th+fillscale_age+ylim(0,60000)+
#'                            theme(legend.position = "none")+ggtitle("Before Quality Control"),
#'                          ggplot(meta, aes(AgeGroup, qc_cell_count,fill=AgeGroup))+
#'                            geom_boxplot()+geom_point()+
#'                            theme_bw()+geom_text(aes(label=Sample_ID), hjust=-0.25, size=3)+xlab("Age Group")+
#'                            ylab("Total Cell Number")+th+fillscale_age+ylim(0,60000)+
#'                            theme(legend.position = "none")+ggtitle("After Quality Control"), ncol=2)
#' 
#' save_plts(cell_count, "IFALD_QC_cellcount_age_PBMC", w=8,h=4)
#' 
#' 
#' d10x <- d10x.list[[1]]
#' 
#' d10x
#' 
#' #saveRDS(d10x, file = here("../../../projects/macparland/RE/PediatricAdult/processed_data","IFALD_d10x_adult_ped_raw_PBMC.rds"))
#' saveRDS(d10x, file = here("/media/redgar/Seagate Portable Drive/processed_data","IFALD_d10x_adult_ped_raw_PBMC.rds"))
#' 
# d10x<-readRDS(file = here("/media/redgar/Seagate Portable Drive/processed_data","IFALD_d10x_adult_ped_raw_PBMC.rds"))
# 
# 
# ################
# ## Normalize scale and UMAP
# ################
# d10x <- NormalizeData(d10x)
# d10x <- FindVariableFeatures(d10x, selection.method = "vst", nfeatures = 2000)
# d10x <- ScaleData(d10x) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))
# 
# # dimension reduction
# d10x <- RunPCA(d10x, ndims.print = 1:10, nfeatures.print = 10)
# d10x <- RunUMAP(d10x, dims = 1:30)
# d10x <- RunTSNE(d10x, dims = 1:30)
# 
# 
# 
# ######################
# ## cell cycle gene expression
# ######################
# # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# # segregate this list into markers of G2/M phase and markers of S phase
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# 
# d10x <- CellCycleScoring(d10x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# 
# pca_cellcycle<-DimPlot(d10x, reduction="pca",  group.by = "Phase")
# save_plts(pca_cellcycle, "IFALD_pca_cellcycle_PBMC", w=6,h=4)
# 
# pca_nfeature<-FeaturePlot(d10x, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
# save_plts(pca_nfeature, "IFALD_pca_nfeature_PBMC", w=6,h=4)
# 
# 
# 
# ## regress out cell cycle and other covariates
# #Transformed data will be available in the SCT assay, which is set as the default after running sctransform
# #By default, sctransform accounts for cellular sequencing depth, or nUMIs.
# d10x <- SCTransform(d10x, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), verbose = FALSE)
# 
# # dimension reduction
# d10x <- RunPCA(d10x, verbose = FALSE)
# d10x <- RunUMAP(d10x, dims = 1:30)
# d10x <- RunTSNE(d10x, dims = 1:30)
# 
# # cluster
# d10x <- FindNeighbors(d10x, reduction = "pca", dims = 1:20)
# d10x <- FindClusters(d10x, resolution = 0.5)
# 
# 
# 
# ###############
# ## visualize
# ###############
# SCT_cluster_umap<-DimPlot(d10x, reduction = "umap", pt.size=0.25, label=T)
# save_plts(SCT_cluster_umap, "IFALD_SCT_cluster_umap_PBMC", w=6,h=4)
# 
# SCT_cluster_tsne<-DimPlot(d10x, reduction = "tsne", pt.size=0.25, label=T)
# save_plts(SCT_cluster_tsne, "IFALD_SCT_cluster_tsne_PBMC", w=6,h=4)
# 
# cell_pca_SCT<-DimPlot(d10x, reduction="pca", group.by="Phase")
# save_plts(cell_pca_SCT, "IFALD_cell_PCA_afterSCT_PBMC", w=6,h=4)
# 
# drop_pca_SCT<-DimPlot(d10x, reduction="pca", group.by="cell_status")
# save_plts(drop_pca_SCT, "IFALD_cell_PCA_afterSCT_DropletQC_PBMC", w=6,h=4)
# 
# nFeature_UMAP_SCT<-FeaturePlot(d10x, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
# save_plts(nFeature_UMAP_SCT, "IFALD_nfeature_UMAP_afterSCT_PBMC", w=6,h=4)
# 
# 
# 
# ## Low quality cluster?
# MT_umap_SCT<-FeaturePlot(d10x, features = "percent.mt", min.cutoff = "q9", pt.size=1)
# save_plts(MT_umap_SCT, "IFALD_MT_umap_PBMC", w=6,h=4)
# 
# ncount_umap_SCT<-FeaturePlot(d10x, features = "nCount_RNA", min.cutoff = "q9", pt.size=1)
# save_plts(ncount_umap_SCT, "IFALD_ncount_umap_SCT_PBMC", w=6,h=4)
# 
# nfeature_umap_SCT<-FeaturePlot(d10x, features = "nFeature_RNA", min.cutoff = "q9", pt.size=1)
# save_plts(nfeature_umap_SCT, "IFALD_nfeature_umap_SCT_PBMC", w=6,h=4)
# 
# ###############
# ## Cell Cluster Labeling
# ###############
# DefaultAssay(d10x)<-"RNA"
# d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")
# 
# signatures<-read.csv(here("data/PBMC_markers.csv"))
# 
# d10x_exp <- GetAssayData(d10x)
# results = SCINA(d10x_exp, signatures, max_iter = 100, convergence_n = 10,
#                 convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
# d10x$SCINA_broad<-results$cell_labels
# gc()
# 
# DimPlot(d10x, reduction = "umap", pt.size=0.25, label=T, group.by = "SCINA_broad")
# DimPlot(d10x, reduction = "umap", pt.size=0.25, label=T)
# 
# signatur_gene<-unique(c(as.character(signatures[2,]),as.character(signatures[1,])))
# DotPlot(object = d10x, features = signatur_gene)
# 
# FeaturePlot(d10x, features=c("CD8A","NKG7","FCGR3B","TRAC","TRDC","CD3D"))
# DotPlot(object = d10x, features = c("FCGR3B","TRAC","TRDC","CD3D","CD27"))
# 
# DimPlot(d10x, reduction = "umap", pt.size=0.25, label=T, group.by = "")
# ggplot(d10x@meta.data, aes(seurat_clusters, nuclear_fraction))+geom_violin()
# ggplot(d10x@meta.data, aes(seurat_clusters, nFeature_RNA))+geom_violin()
# ggplot(d10x@meta.data, aes(seurat_clusters, nCount_RNA))+geom_violin()
# 
# 
# ## Sum SCINA to cluster level
# d10x@meta.data$CellType_refined<-as.character(d10x@meta.data$SCINA_broad)
# d10x@meta.data$CellType_refined[which(d10x@meta.data$seurat_clusters%in%c("1","2","6","9","10"))]<-"Mature B-cells"
# d10x@meta.data$CellType_refined[which(d10x@meta.data$seurat_clusters%in%c("16"))]<-"Plasma cells"
# d10x@meta.data$CellType_refined[which(d10x@meta.data$seurat_clusters%in%c("0","13"))]<-"Naive CD4 T-cells"
# d10x@meta.data$CellType_refined[which(d10x@meta.data$seurat_clusters%in%c("5","11"))]<-"Memory CD4 T-cells"
# d10x@meta.data$CellType_refined[which(d10x@meta.data$seurat_clusters%in%c("3","7"))]<-"CD8 T-cells"
# d10x@meta.data$CellType_refined[which(d10x@meta.data$seurat_clusters%in%c("4","8","14"))]<-"NK cells"
# d10x@meta.data$CellType_refined[which(d10x@meta.data$seurat_clusters%in%c("15"))]<-"DC"
# d10x@meta.data$CellType_refined[which(d10x@meta.data$seurat_clusters%in%c("12"))]<-"Neutrophil"
# 
# DimPlot(d10x, reduction = "umap", pt.size=0.25, label=T, group.by = "CellType_refined")


###############
## Integration with liver samples
###############
#https://satijalab.org/seurat/articles/integration_rpca.html
print("RUNNING INTEGRATION")


#d10x_PBMC<-readRDS(file = here("/media/redgar/Seagate Portable Drive/processed_data","IFALD_d10x_adult_ped_raw_PBMC.rds"))

d10x_PBMC<-readRDS(file = here("../../../projects/macparland/RE/PediatricAdult/processed_data","IFALD_d10x_adult_ped_raw_PBMC.rds"))
d10x_liver<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))

d10x.list<- list(PBMC = d10x_PBMC, liver = d10x_liver)

## run integration across donor and hopefully that will also smooth out differences with chemistry?
## data is already split by donor

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
save_plts(SCT_cluster_umap, "IFALD_rPCA_cluster_umap_PBMC", w=6,h=4)

SCT_cluster_tsne<-DimPlot(d10x.combined, reduction = "tsne", pt.size=0.25, label=T)
save_plts(SCT_cluster_tsne, "IFALD_rPCA_cluster_tsne_PBMC", w=6,h=4)

chem_umap_sct<-DimPlot(d10x.combined, reduction = "umap", group.by = "chemistry", pt.size=0.25)
save_plts(chem_umap_sct, "IFALD_chem_rPCA_umap_PBMC", w=6,h=4)

age_umap_sct<-DimPlot(d10x.combined, reduction = "umap", group.by = "AgeGroup", pt.size=0.25)+fillscale_age
save_plts(age_umap_sct, "IFALD_age_rPCA_umap_PBMC", w=6,h=4)

MT_umap_sct<-FeaturePlot(d10x.combined, reduction = "umap", features = "percent.mt", pt.size=0.25)
save_plts(MT_umap_sct, "IFALD_MT_rPCA_umap_PBMC", w=5,h=4)

individual_umap_sct<-DimPlot(d10x.combined, reduction = "umap", group.by = "individual", pt.size=0.25)
save_plts(individual_umap_sct, "IFALD_individual_rPCA_UMAP_PBMC", w=6,h=4)

individual_split<-DimPlot(d10x.combined, reduction = "umap", group.by = "individual", split.by="individual",pt.size=0.25)
save_plts(individual_split, "IFALD_individual_facet_rPCA_UMAP_PBMC", w=20,h=4)

age_split<-DimPlot(d10x.combined, reduction = "umap", group.by = "individual", split.by="AgeGroup",pt.size=0.25)
save_plts(age_split, "IFALD_age_facet_rPCA_UMAP_PBMC", w=10,h=5)


save(d10x.combined, file=paste(here("data/"),"IFALD_adult_ped_PBMC_integrated.rds", sep=""))


print(sessionInfo())




