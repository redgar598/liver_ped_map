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

#' #dataset_loc <- here("/media/redgar/Seagate Portable Drive/ped_liver_map_raw")
#' dataset_loc <- here("../../../projects/macparland/RE/PediatricAdult/ped_liver_map_raw")
#' 
#' samples<-list.files(dataset_loc)
#' samples<-samples[-grep("data_transfer",samples)]
#' print(samples)
#' 
#' #meta<-read.table(here("data/data_transfer_updated_mar20_2023_IFALD.csv"), header=T, sep=",")
#' meta<-read.table(here(dataset_loc,"data_transfer_updated_mar20_2023_IFALD.csv"), header=T, sep=",")
#' 
#' samples<-samples[which(samples%in%meta$file)]
#' 
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
#' # Show QC metrics for the first 5 cells
#' print(head(d10x.list[[2]]@meta.data, 5))
#' 
#' #'Low-quality cells or empty droplets will often have very few genes
#' #'Cell doublets or multiplets may exhibit an aberrantly high gene count
#' 
#' 
#' # Visualize QC metrics
#' #nFeature number of unique genes
#' #nCount number of total molecules
#' plt_QC_data<-do.call(rbind, lapply(1:length(d10x.list), function(x) d10x.list[[x]]@meta.data))
#' save(plt_QC_data, file=here("data","IFALD_QC_metrics.Rdata"))
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
#' save_plts(qc_plts, "IFALD_intital_QC_plts", w=6,h=4)
#' 
#' qc_plts_chem<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
#'   geom_point() + facet_wrap(~chemistry)+
#'   scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
#'   geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
#'   geom_hline(yintercept = 6000) +theme_bw()+th
#' save_plts(qc_plts_chem, "IFALD_intital_QC_plts_chemistry", w=12,h=4)
#' 
#' qc_plts_individual<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
#'   geom_point() + facet_wrap(~individual)+
#'   scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
#'   geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
#'   geom_hline(yintercept = 6000) +theme_bw()+th
#' save_plts(qc_plts_chem, "IFALD_intital_QC_plts_individual", w=12,h=4)
#' 
#' MT_plt<-ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
#'   geom_vline(xintercept = 25)+ theme_bw()+xlab("Percent Mitochondrial")+th
#' save_plts(MT_plt, "IFALD_percentMT_plt", w=6,h=4)
#' 
#' MT_plt_individual<-ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
#'   facet_wrap(~individual, scales="free_y")+
#'   geom_vline(xintercept = 25)+ theme_bw()+xlab("Percent Mitochondrial")+th
#' save_plts(MT_plt_individual, "IFALD_percentMT_plt_individual", w=8,h=4)
#' 
#' 
#' 
#' nf_plts<-ggplot(plt_QC_data, aes(nuclear_fraction,log10(nCount_RNA),colour=cell_status)) +
#'   geom_point() +  ylab("Number of Total Molecules\n(log 10 nCount) ")+xlab("Nuclear Fraction")+theme_bw()+th+
#'   facet_wrap(~individual)
#' nf_plts
#' save_plts(nf_plts, "IFALD_nuclear_fraction", w=6,h=4)
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
#' save_plts(cell_count, "IFALD_QC_cellcount_age", w=8,h=4)
#' 
#' 
#' d10x <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "adult_ped_map")#add.cell.ids = alldata_names2,
#' 
#' d10x
#' 
#' saveRDS(d10x, file = here("data","IFALD_d10x_adult_ped_raw.rds"))
#' 
#' 
#' ################
#' ## Normalize scale and UMAP
#' ################
#' d10x <- NormalizeData(d10x)
#' d10x <- FindVariableFeatures(d10x, selection.method = "vst", nfeatures = 2000)
#' d10x <- ScaleData(d10x) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))
#' 
#' # dimension reduction
#' d10x <- RunPCA(d10x, ndims.print = 1:10, nfeatures.print = 10)
#' d10x <- RunUMAP(d10x, dims = 1:30)
#' d10x <- RunTSNE(d10x, dims = 1:30)
#' 
#' 
#' 
#' ######################
#' ## cell cycle gene expression
#' ######################
#' # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
#' # segregate this list into markers of G2/M phase and markers of S phase
#' s.genes <- cc.genes$s.genes
#' g2m.genes <- cc.genes$g2m.genes
#' 
#' d10x <- CellCycleScoring(d10x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#' 
#' pca_cellcycle<-DimPlot(d10x, reduction="pca",  group.by = "Phase")
#' save_plts(pca_cellcycle, "IFALD_pca_cellcycle", w=6,h=4)
#' 
#' pca_nfeature<-FeaturePlot(d10x, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
#' save_plts(pca_nfeature, "IFALD_pca_nfeature", w=6,h=4)
#' 
#' 
#' 
#' ## regress out cell cycle and other covariates
#' #Transformed data will be available in the SCT assay, which is set as the default after running sctransform
#' #By default, sctransform accounts for cellular sequencing depth, or nUMIs.
#' d10x <- SCTransform(d10x, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), verbose = FALSE)
#' 
#' # dimension reduction
#' d10x <- RunPCA(d10x, verbose = FALSE)
#' d10x <- RunUMAP(d10x, dims = 1:30)
#' d10x <- RunTSNE(d10x, dims = 1:30)
#' 
#' # cluster
#' d10x <- FindNeighbors(d10x, reduction = "pca", dims = 1:20)
#' d10x <- FindClusters(d10x, resolution = 0.5)
#' 
#' 
#' 
#' ###############
#' ## visualize
#' ###############
#' SCT_cluster_umap<-DimPlot(d10x, reduction = "umap", pt.size=0.25, label=T)
#' save_plts(SCT_cluster_umap, "IFALD_SCT_cluster_umap", w=6,h=4)
#' 
#' SCT_cluster_tsne<-DimPlot(d10x, reduction = "tsne", pt.size=0.25, label=T)
#' save_plts(SCT_cluster_tsne, "IFALD_SCT_cluster_tsne", w=6,h=4)
#' 
#' cell_pca_SCT<-DimPlot(d10x, reduction="pca", group.by="Phase")
#' save_plts(cell_pca_SCT, "IFALD_cell_PCA_afterSCT", w=6,h=4)
#' 
#' drop_pca_SCT<-DimPlot(d10x, reduction="pca", group.by="cell_status")
#' save_plts(drop_pca_SCT, "IFALD_cell_PCA_afterSCT_DropletQC", w=6,h=4)
#' 
#' nFeature_UMAP_SCT<-FeaturePlot(d10x, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
#' save_plts(nFeature_UMAP_SCT, "IFALD_nfeature_UMAP_afterSCT", w=6,h=4)
#' 
#' 
#' 
#' chem_umap_sct<-DimPlot(d10x, reduction = "umap", group.by = "chemistry", pt.size=0.25)
#' save_plts(chem_umap_sct, "IFALD_chem_SCT_umap", w=6,h=4)
#' 
#' age_umap_sct<-DimPlot(d10x, reduction = "umap", group.by = "AgeGroup", pt.size=0.25)+fillscale_age
#' save_plts(age_umap_sct, "IFALD_age_SCT_umap", w=6,h=4)
#' 
#' individual_umap_sct<-DimPlot(d10x, reduction = "umap", group.by = "individual", pt.size=1)
#' save_plts(individual_umap_sct, "IFALD_individual_SCT_UMAP", w=6,h=4)
#' 
#' 
#' ## Low quality cluster?
#' MT_umap_SCT<-FeaturePlot(d10x, features = "percent.mt", min.cutoff = "q9", pt.size=1)
#' save_plts(pca_nfeature, "IFALD_pca_nfeature", w=6,h=4)
#' 
#' ncount_umap_SCT<-FeaturePlot(d10x, features = "nCount_RNA", min.cutoff = "q9", pt.size=1)
#' save_plts(ncount_umap_SCT, "IFALD_ncount_umap_SCT", w=6,h=4)
#' 
#' nfeature_umap_SCT<-FeaturePlot(d10x, features = "nFeature_RNA", min.cutoff = "q9", pt.size=1)
#' save_plts(nfeature_umap_SCT, "IFALD_nfeature_umap_SCT", w=6,h=4)
#' 
#' 
#' 
#' 
#'                                 #
#'                                 # ###############
#'                                 # ## Integration
#'                                 # ###############
#'                                 # #https://satijalab.org/seurat/articles/integration_rpca.html
#'                                 # print("RUNNING INTEGRATION")
#'                                 # ## Start back with raw data
#'                                 # d10x <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "adult_ped_map")#add.cell.ids = alldata_names2,
#'                                 # d10x
#'                                 #
#'                                 # # split the dataset into a list of two seurat objects (3' and 5')
#'                                 # d10x.list.chem <- SplitObject(d10x, split.by = "chemistry")
#'                                 #
#'                                 #
#'                                 # # normalize, identify variable features and score cell cycle for each dataset independently
#'                                 # s.genes <- cc.genes$s.genes
#'                                 # g2m.genes <- cc.genes$g2m.genes
#'                                 #
#'                                 # d10x.list.chem <- lapply(X = d10x.list.chem, FUN = function(x) {
#'                                 #   x <- NormalizeData(x)
#'                                 #   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
#'                                 #   x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#'                                 # })
#'                                 #
#'                                 # # select features that are repeatedly variable across datasets for integration run PCA on each
#'                                 # # dataset using these features
#'                                 # features <- SelectIntegrationFeatures(object.list = d10x.list.chem)
#'                                 # d10x.list.chem <- lapply(X = d10x.list.chem, FUN = function(x) {
#'                                 #   #x <- ScaleData(x, features = features, verbose = FALSE)
#'                                 #   x <- ScaleData(x, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), features = features, verbose = FALSE)
#'                                 #   x <- RunPCA(x, features = features, verbose = FALSE)
#'                                 # })
#'                                 #
#'                                 #
#'                                 #
#'                                 # ## Identify anchors
#'                                 # chem.anchors <- FindIntegrationAnchors(object.list = d10x.list.chem, anchor.features = features, reduction = "rpca")
#'                                 # d10x.combined <- IntegrateData(anchorset = chem.anchors)
#'                                 #
#'                                 # DefaultAssay(d10x.combined) <- "integrated"
#'                                 #
#'                                 # print("INTEGRATED")
#' 
#' 
#' ###############
#' ## Integration
#' ###############
#' #https://satijalab.org/seurat/articles/integration_rpca.html
#' print("RUNNING INTEGRATION")
#' 
#' ## run integration across donor and hopefully that will also smooth out differences with chemistry?
#' ## data is already split by donor
#' 
#' # normalize, identify variable features and score cell cycle for each dataset independently
#' s.genes <- cc.genes$s.genes
#' g2m.genes <- cc.genes$g2m.genes
#' 
#' d10x.list <- lapply(X = d10x.list, FUN = function(x) {
#'   x <- NormalizeData(x)
#'   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
#'   x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#' })
#' 
#' # select features that are repeatedly variable across datasets for integration run PCA on each
#' # dataset using these features
#' features <- SelectIntegrationFeatures(object.list = d10x.list)
#' d10x.list <- lapply(X = d10x.list, FUN = function(x) {
#'   #x <- ScaleData(x, features = features, verbose = FALSE)
#'   x <- ScaleData(x, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), features = features, verbose = FALSE)
#'   x <- RunPCA(x, features = features, verbose = FALSE)
#' })
#' 
#' 
#' 
#' ## Identify anchors
#' chem.anchors <- FindIntegrationAnchors(object.list = d10x.list, anchor.features = features, reduction = "rpca")
#' d10x.combined <- IntegrateData(anchorset = chem.anchors)
#' 
#' DefaultAssay(d10x.combined) <- "integrated"
#' 
#' print("INTEGRATED")
#' 
#' 
#' # Run the standard workflow for visualization and clustering
#' d10x.combined <- ScaleData(d10x.combined, verbose = FALSE)
#' d10x.combined <- RunPCA(d10x.combined, npcs = 30, verbose = FALSE)
#' d10x.combined <- RunUMAP(d10x.combined, reduction = "pca", dims = 1:30)
#' d10x.combined <- RunTSNE(d10x.combined, dims = 1:30)
#' 
#' d10x.combined <- FindNeighbors(d10x.combined, reduction = "pca", dims = 1:30)
#' d10x.combined <- FindClusters(d10x.combined, resolution = 0.5)
#' 
#' d10x.combined
#' 
#' 
#' ###########
#' ## Visualize integration
#' ###########
#' SCT_cluster_umap<-DimPlot(d10x.combined, reduction = "umap", pt.size=0.25, label=T)
#' save_plts(SCT_cluster_umap, "IFALD_rPCA_cluster_umap", w=6,h=4)
#' 
#' SCT_cluster_tsne<-DimPlot(d10x.combined, reduction = "tsne", pt.size=0.25, label=T)
#' save_plts(SCT_cluster_tsne, "IFALD_rPCA_cluster_tsne", w=6,h=4)
#' 
#' chem_umap_sct<-DimPlot(d10x.combined, reduction = "umap", group.by = "chemistry", pt.size=0.25)
#' save_plts(chem_umap_sct, "IFALD_chem_rPCA_umap", w=6,h=4)
#' 
#' age_umap_sct<-DimPlot(d10x.combined, reduction = "umap", group.by = "AgeGroup", pt.size=0.25)+fillscale_age
#' save_plts(age_umap_sct, "IFALD_age_rPCA_umap", w=6,h=4)
#' 
#' MT_umap_sct<-FeaturePlot(d10x.combined, reduction = "umap", features = "percent.mt", pt.size=0.25)
#' save_plts(MT_umap_sct, "IFALD_MT_rPCA_umap", w=5,h=4)
#' 
#' individual_umap_sct<-DimPlot(d10x.combined, reduction = "umap", group.by = "individual", pt.size=0.25)
#' save_plts(individual_umap_sct, "IFALD_individual_rPCA_UMAP", w=6,h=4)
#' 
#' individual_split<-DimPlot(d10x.combined, reduction = "umap", group.by = "individual", split.by="individual",pt.size=0.25)
#' save_plts(individual_split, "IFALD_individual_facet_rPCA_UMAP", w=20,h=4)
#' 
#' age_split<-DimPlot(d10x.combined, reduction = "umap", group.by = "individual", split.by="AgeGroup",pt.size=0.25)
#' save_plts(age_split, "IFALD_age_facet_rPCA_UMAP", w=10,h=5)
#' 
# 
# ############################
# ##### Example markers to plot for the liver [From Diana]
# ############################
# 
# ######
# #Immune cells
# ######
# 
# ## Macrophages
# DotPlot(d10x.combined, features = c( "PTPRC", "CD68", "MARCO","CD5L","VSIG4", "MAF", "LYZ", "CSTA", "S100A8", "S10049", "CD14", "CD74", "GPBAR1", "ID3"), cols=c("blue", "red")) + RotatedAxis()
# 
# ## NK/T/B cells
# DotPlot(d10x.combined, features = c("PTPRC", "CD2", "CD3E", "IL7R", "KLRB1",
#                                     "NKG7", "GZMA", "GZMB", "GZMK" , "PRF1","CD4", "CD8A","CD247", "TRAC","TRDC", "TRGC1", "TRGC2", "TRBC1",
#                                     "TRBC2", "S1PR1", "CD28", "CD27", "SELL", "CCR7", "CXCR4","CCR4","FAS",
#                                     "FOXP3", "CTLA4", "LAG3", "TNFRSF4","TNFRSF18", "ICOS" ,"CD69", "CD79A", "CD79B", "IGHG1", "MS4A1",
#                                     "LTB", "CD52", "IGHD", "CD19", "ID3"), cols=c("blue", "red")) + RotatedAxis()
# 
# # LEC and LSEC
# DotPlot(d10x.combined, features = c("CALCRL", "VWF", "RAMP2", "STAB2", "LYVE1", "PECAM1", "ENG", "FCGR2B", "F8", "SPARCL1", "ID1", "SOX18", "CD32B", "ID3"), cols=c("blue", "red")) + RotatedAxis()
# 
# ######
# #Hepatocytes
# ######
# DotPlot(d10x.combined, features=c("ALB", "HAMP", "ARG1", "PCK1", "AFP", "BCHE", "HAL", "SCD", "CPS1", "CYP3A4",
#                                   "ELF3", "CRP", "GSTA2", "AKR1C1", "MGST1", "CYP3A5", "ALDH1A1", "ADH1A", "CYP2E1",
#                                   "GLS2", "SDS", "GLUL", "AKR1D1", "HPR",
#                                   "HMGCS1", "IGSF23", "ACSS2", "G6PC", "ID3"),
#         cols=c("blue", "red")) + RotatedAxis()
# 
# 
# 
# ######
# #Cholangiocytes
# ######
# DotPlot(d10x.combined,features = c( "EPCAM", "SOX9", "KRT1", "KRT7", "ANXA4", "KRT18", "ID3"), cols=c("blue", "red")) + RotatedAxis()
# 
# ######
# #HSCs
# ######
# DotPlot(d10x.combined,features = c( "RBP1", "LRAT", "PDE3B", "ACTA2", "AOX1", "PDE3D", "PDE4D", "SPARC", "TAGLN", "COL1A1", "COL1A2", "COL3A1", "TIMP1", "DCN", "MYL9", "TPM2", "MEG3", "BGN", "IGFBP7", "IGFBP3", "CYR61", "IGFBP6", "CCL2", "COLEC11", "CTGF", "HGF", "ID3"), cols=c("blue", "red")) + RotatedAxis()
# FeaturePlot(d10x.combined, reduction = "umap", features = c("PTPRC", "CD3D", "CD68", "CD79A","TRDC", "NKG7", "KRT7", "CALCRL", "ACTA2", "MS4A1", "CYP3A4", "SCD", "FCN2", "CD4", "CD8A", "FCER1A", "MARCO", "LYZ", "VSIG4", "FOLR2", "ID3"), ncol = 4)
# key_markers<-FeaturePlot(d10x.combined, reduction = "umap", features = c("EPCAM", "SOX9", "SELL", "PTPRC",
#                                                                          "TRDC", "NKG7", "CALCRL", "VWF",
#                                                                          "MARCO", "LYZ","COL1A1","IGFBP3"), ncol = 4)
# save_plts(key_markers, "IFALD_markers_rPCA_UMAP", w=25,h=20)
# 
#################
## Rough annotation
#################
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

# genes<-unique(c(Macrophage_genes, NK_T_genes, B_genes,LEC_genes,Hepatocyte_genes,Cholangiocytes_genes,HSCs_genes))
# 
# d10x.exp<-as.data.frame(d10x.combined[["RNA"]]@data)
# d10x.exp.GOI<-d10x.exp[genes,]
# d10x.exp.GOI$gene<-rownames(d10x.exp.GOI)
# d10x.exp.GOI<-melt(d10x.exp.GOI)#
# 
# meta<-d10x.combined@meta.data
# meta$cell<-rownames(meta)
# 
# plt<-merge(d10x.exp.GOI, meta,by.x="variable", by.y="cell")
# 
# plt$variable<-as.character(plt$variable)
# 
# ## possible cells types
# cluster_marker_mean<-function(gene_list, type){
#   plt_epi<-plt[which(plt$gene%in%gene_list),]
#   mean_type<-as.data.frame(tapply(plt_epi$value, plt_epi$seurat_clusters, mean))
#   colnames(mean_type)<-type
#   mean_type
# }
# 
# cell_rough<-cbind(cluster_marker_mean(Macrophage_genes, "Myeloid"),
#                   cluster_marker_mean(NK_T_genes, "NK_T"),
#                   cluster_marker_mean(LEC_genes, "LEC"),
#                   cluster_marker_mean(Hepatocyte_genes, "Hepatocyte"),
#                   cluster_marker_mean(Cholangiocytes_genes, "Cholangiocytes"),
#                   cluster_marker_mean(HSCs_genes, "HSC"),
#                   cluster_marker_mean(B_genes, "B_cell"))
# 
# 
# cell_rough$CellType_rough<-sapply(1:nrow(cell_rough), function(x) {
#   compart<-colnames(cell_rough)[which(cell_rough[x,] == max(cell_rough[x,]))]
#   if(length(compart)==1){compart}else{"Unclear"}
# })
# 
# cell_rough$seurat_clusters<-rownames(cell_rough)
# 
# ## to explore max values
# save(cell_rough, file=here("data/IFALD_cell_rough_maxmean.RData"))
# 
# ## second option to plot (to maybe manually relable hepatocytes)
# not_hep_cell<-colnames(cell_rough)[which(!(colnames(cell_rough)%in%c("Hepatocyte","CellType_rough" ,"seurat_clusters")))]
# 
# cell_rough$second_best_cell<-sapply(1:nrow(cell_rough), function(x){
#   if(cell_rough$CellType_rough[x]!="Hepatocyte"){cell_rough$CellType_rough[x]}else{
#     not_hep_mean_max<-max(cell_rough[x, not_hep_cell])
#     paste("Hep_",not_hep_cell[which(cell_rough[x, not_hep_cell]==not_hep_mean_max)], sep="")}
# })
# 
# meta<-d10x.combined@meta.data
# meta$cell<-rownames(meta)
# plt_summary<-merge(meta, cell_rough[,c("seurat_clusters","CellType_rough","second_best_cell")], by="seurat_clusters")
# plt_summary<-plt_summary[match(rownames(d10x.combined@meta.data),plt_summary$cell),]
# identical(plt_summary$cell, rownames(d10x.combined@meta.data))
# 
# rownames(plt_summary)<-plt_summary$cell
# 
# d10x.combined<- AddMetaData(d10x.combined, plt_summary)
# d10x.combined
# table(d10x.combined$CellType_rough)
# 
# table(d10x.combined$second_best_cell, d10x.combined$CellType_rough)
# 
# 
# 
# ##############
# ## Save integrated to look local
# ##############
# save(d10x.combined, file=paste(here("data/"),"IFALD_adult_ped_integrated.rds", sep=""))
# cell_label<-d10x.combined@meta.data
# save(cell_label, file=paste(here("data/"),"IFALD_adult_ped_cellRough.rds", sep=""))
# 
# 
# load(here("data","IFALD_adult_ped_integrated.rds"))
# 
# 
# #####
# ## plot cell types
# #####
# d10x.combined@meta.data$CellType_rough<-as.factor(d10x.combined@meta.data$CellType_rough)
# levels(d10x.combined@meta.data$CellType_rough)<-c("B-cells","Cholangiocytes",
#                                                   "Hepatocytes",
#                                                   "HSC","LSEC","Myeloid cells","NK and T cells")
# 
# roughcell_cluster_umap<-DimPlot(d10x.combined, reduction = "umap",group.by="CellType_rough", pt.size=0.15, label=T)+colscale_cellType+ggtitle("")+
#   annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(d10x.combined))))
# roughcell_cluster_umap
# save_plts(roughcell_cluster_umap, "IFALD_rPCA_roughcellType_cluster_umap", w=6,h=4)
# roughcell_cluster_umap<-DimPlot(d10x.combined, reduction = "umap",group.by="CellType_rough", pt.size=0.15)+colscale_cellType+ggtitle("")+
#   annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(d10x.combined))))
# roughcell_cluster_umap
# save_plts(roughcell_cluster_umap, "IFALD_rPCA_roughcellType_cluster_umap_nolab", w=6,h=4)
# 
# roughcell_cluster_tsne<-DimPlot(d10x.combined, reduction = "tsne",group.by="CellType_rough", pt.size=0.15)+colscale_cellType+ggtitle("")+
#   annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(d10x.combined))))
# roughcell_cluster_tsne
# save_plts(roughcell_cluster_tsne, "IFALD_rPCA_roughcellType_cluster_tsne_nolab", w=6,h=4)
# 
# ########
# ## DropletQC
# ########
# nuclearfraction_cluster_umap<-FeaturePlot(d10x.combined, reduction = "umap",feature="nuclear_fraction", pt.size=0.15, label=F)+ggtitle("")+
#   annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(d10x.combined))))
# nuclearfraction_cluster_umap
# save_plts(nuclearfraction_cluster_umap, "IFALD_rPCA_nuclearfraction_cluster_umap", w=6,h=4)
# damaged_empty_cluster_umap<-DimPlot(d10x.combined, reduction = "umap",group.by="cell_status", pt.size=0.15)+ggtitle("")+
#   annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(d10x.combined))))+
#   scale_color_manual(values=c("grey","red","cornflowerblue"),name="Nuclear\nFraction")
# damaged_empty_cluster_umap
# save_plts(damaged_empty_cluster_umap, "IFALD_rPCA_damaged_empty_cluster_umap", w=6,h=4)
# 
# #other QC
# nFeature_cluster_umap<-FeaturePlot(d10x.combined, reduction = "umap",feature="nFeature_RNA", pt.size=0.15, label=F)+ggtitle("")+
#   annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(d10x.combined))))
# nFeature_cluster_umap
# save_plts(nFeature_cluster_umap, "IFALD_rPCA_nFeature_cluster_umap", w=6,h=4)
# 



#################
## SCINA
#################
d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

signatures<-read.csv(here("data/Liver_Markers - Human_for_SCINA.csv"))

d10x_exp <- GetAssayData(d10x)
results = SCINA(d10x_exp, signatures, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
d10x$SCINA_broad<-results$cell_labels

### Cell subsets
RBCsignatures<-read.csv(here("data/Liver_Markers - Erythrocytes.csv"))
Neutrosignatures<-read.csv(here("data/Liver_Markers - Neurtophil.csv"))
Tsignatures<-read.csv(here("data/Liver_Markers - Tcell.csv"))
Bsignatures<-read.csv(here("data/Liver_Markers - Bcell.csv"))
myeloidsignatures<-read.csv(here("data/Liver_Markers - Myeloid.csv"))

d10x.combined_RBC<-subset(d10x, subset = SCINA_broad == "Erythrocytes")
d10x_exp <- GetAssayData(d10x.combined_RBC)
results = SCINA(d10x_exp, RBCsignatures, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
d10x.combined_RBC$SCINA_refined<-results$cell_labels

d10x.combined_myeloid<-subset(d10x, subset = SCINA_broad == "Myeloid")
d10x_exp <- GetAssayData(d10x.combined_myeloid)
results = SCINA(d10x_exp, myeloidsignatures, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
d10x.combined_myeloid$SCINA_refined<-results$cell_labels

d10x.combined_bcell<-subset(d10x, subset = SCINA_broad == "B_cell")
d10x_exp <- GetAssayData(d10x.combined_bcell)
results = SCINA(d10x_exp, Bsignatures, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
d10x.combined_bcell$SCINA_refined<-results$cell_labels

d10x.combined_tcell<-subset(d10x, subset = SCINA_broad == "T_cell")
d10x_exp <- GetAssayData(d10x.combined_tcell)
results = SCINA(d10x_exp, Tsignatures, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
d10x.combined_tcell$SCINA_refined<-results$cell_labels

d10x.combined_neutro<-subset(d10x, subset = SCINA_broad == "Neutrophil")
d10x_exp <- GetAssayData(d10x.combined_neutro)
results = SCINA(d10x_exp, Neutrosignatures, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
d10x.combined_neutro$SCINA_refined<-results$cell_labels

SCINA_cell_labels<-rbind(d10x.combined_RBC@meta.data[,c("cell","SCINA_broad","SCINA_refined")],
                         d10x.combined_myeloid@meta.data[,c("cell","SCINA_broad","SCINA_refined")],
                         d10x.combined_bcell@meta.data[,c("cell","SCINA_broad","SCINA_refined")],
                         d10x.combined_tcell@meta.data[,c("cell","SCINA_broad","SCINA_refined")],
                         d10x.combined_neutro@meta.data[,c("cell","SCINA_broad","SCINA_refined")])
save(SCINA_cell_labels, file=paste(here("data/"),"IFALD_adult_ped_SCINA_cell_labels.rds", sep=""))


# 
# ########
# ## refining cell labels
# ########
# d10x.combined@meta.data$second_best_cell<-as.factor(d10x.combined@meta.data$second_best_cell)
# levels(d10x.combined@meta.data$second_best_cell)<-c("B-cells","Cholangiocytes","Cholangiocytes\n(Hepatocyte Like)",
#                                                    "Myeloid cells\n(Hepatocyte Like)",# "LSEC\n(Hepatocyte Like)", "HSC\n(Hepatocyte Like)",
#                                                     "HSC","LSEC","Myeloid cells","NKT cells")#"NKT cells\n(Hepatocyte Like)",
# 
# 
# secondbest_cluster_umap<-DimPlot(d10x.combined, reduction = "umap",group.by="second_best_cell", pt.size=0.15, label=T)+colscale_cellType+ggtitle("")+
#   annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(d10x.combined))))
# secondbest_cluster_umap
# save_plts(secondbest_cluster_umap, "IFALD_rPCA_secondbest_cellType_cluster_umap", w=6,h=4)
# secondbest_cluster_umap<-DimPlot(d10x.combined, reduction = "umap",group.by="second_best_cell", pt.size=0.15)+colscale_cellType+ggtitle("")+
#   annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(d10x.combined))))
# secondbest_cluster_umap
# save_plts(secondbest_cluster_umap, "IFALD_rPCA_secondbest__cluster_umap_nolab", w=6,h=4)
# 
# 
# DefaultAssay(d10x.combined) <- "RNA"
# 
# pdf(file = here("figures/IFALD_dot_plots.pdf"), w=10, h=10)
# DotPlot(object = d10x.combined, features = B_genes)+xlab("B Cell Marker")
# DotPlot(object = d10x.combined, features = T_genes)+xlab("T Cell Marker")
# DotPlot(object = d10x.combined, features = NK_genes)+xlab("NK Cell Marker")
# DotPlot(object = d10x.combined, features = LEC_genes)+xlab("LSEC Marker")
# DotPlot(object = d10x.combined, features = Hepatocyte_genes[1:15])+xlab("Hepatocyte Marker")
# DotPlot(object = d10x.combined, features = Hepatocyte_genes[16:29])+xlab("Hepatocyte Marker")
# DotPlot(object = d10x.combined, features = Cholangiocytes_genes)+xlab("Cholangiocyte Marker")
# DotPlot(object = d10x.combined, features = HSCs_genes[1:14])+xlab("HSC Marker")
# DotPlot(object = d10x.combined, features = HSCs_genes[15:27])+xlab("HSC Marker")
# DotPlot(object = d10x.combined, features = Macrophage_genes)+xlab("Macrophage Marker")
# dev.off()
# 
# 
# # ## donor integrated map, clusters needing checks: 12, 16, 19, 26, 31
# mt_percent_box<-ggplot(d10x.combined@meta.data, aes(seurat_clusters, percent.mt))+geom_violin( fill='lightgrey')+theme_bw()
# save_plts(mt_percent_box, "IFALD_mt_percent_box", w=12,h=2)
# 
#             ggplot(d10x.combined@meta.data, aes(seurat_clusters, nuclear_fraction))+geom_violin()
#             table(d10x.combined@meta.data$seurat_clusters, d10x.combined@meta.data$cell_status)
#             plt_entropy_individual<-entropy_d10(d10x.combined, "individual")
#             entropy_individual<-entropy_plt(plt_entropy_individual, "individual", d10x.combined)
#             entropy_individual
#             plt_entropy_individual<-entropy_d10(d10x.combined, "cell_status")
#             entropy_individual<-entropy_plt(plt_entropy_individual, "cell_status", d10x.combined)
#             entropy_individual
# 
# # 26 Highly ALB contaminated myeloid
# hep_QC<-FeaturePlot(d10x.combined, reduction = "umap",features="ALB",pt.size=0.15)
# save_plts(hep_QC, "IFALD_ALB_rPCA_UMAP", w=6,h=5)
# 
# # 19 seems to be red blood cells
# RBC<-FeaturePlot(d10x.combined, features = c("HBB","HBA2","HBA1","FCGR3A"), min.cutoff = "q9", pt.size=0.15)
# save_plts(RBC, "IFALD_RBC_rPCA_UMAP", w=8,h=7)
# 
# DimPlot(d10x.combined, reduction = "umap")+scale_color_manual(values=c(rep("grey",19),"red",rep("grey", 12)))
# DimPlot(d10x.combined, reduction = "tsne")+scale_color_manual(values=c(rep("grey",19),"red",rep("grey", 12)))
# 
# DimPlot(d10x.combined, reduction = "umap")+scale_color_manual(values=c(rep("grey",16),"red",rep("grey", 15)))
# DimPlot(d10x.combined, reduction = "umap")+scale_color_manual(values=c(rep("grey",26),"red",rep("grey", 5)))
# DimPlot(d10x.combined, reduction = "umap")+scale_color_manual(values=c(rep("grey",29),"red",rep("grey", 2)))
# DimPlot(d10x.combined, reduction = "umap")+scale_color_manual(values=c(rep("grey",24),"red",rep("grey", 7)))
# 
# DimPlot(d10x.combined, reduction = "umap")+scale_color_manual(values=c(rep("grey",18),"red",rep("grey", 13)))
# DimPlot(d10x.combined, reduction = "tsne")+scale_color_manual(values=c(rep("grey",18),"red",rep("grey", 13)))
# DimPlot(d10x.combined, reduction = "umap")+scale_color_manual(values=c(rep("grey",31),"red"))
# 
# 
# #############
# # #relabel some clusters
# #############
# DefaultAssay(d10x.combined) <- "integrated"
# d10x.combined@meta.data$CellType_rough<-as.character(d10x.combined@meta.data$CellType_rough)
# d10x.combined@meta.data$CellType_rough[which(d10x.combined@meta.data$seurat_clusters%in%c("26"))]<-"Low Quality"
# d10x.combined@meta.data$CellType_rough[which(d10x.combined@meta.data$seurat_clusters=="19")]<-"Erythrocytes"
# 
# 
# DimPlot(d10x.combined, reduction = "umap",group.by="CellType_rough", pt.size=0.15, label=T)+colscale_cellType+ggtitle("")+
#   annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(d10x.combined))))
# 
# 
# # some are TPSAB1+, AREG+ resident Mast cells
# MAST_markers<-FeaturePlot(d10x.combined, features = c("TPSAB1", "AREG"))
# save_plts(MAST_markers, "IFALD_MAST_markers_rPCA_UMAP", w=8,h=4)
# 
# MAST_markers<-FeaturePlot(d10x.combined, features = c("TPSAB1", "AREG"))
# save_plts(MAST_markers, "IFALD_MAST_markers_rPCA_UMAP_fullmap", w=8,h=4)
# 
# 
# ##########
# ## Individual plots
# ##########
# d10x.combined$age_id<-paste(d10x.combined$individual, d10x.combined$AgeGroup)
# #d10x.combined$age_id[which(d10x.combined$age_id=="C85_caud5pr Ped")]<-"C85_caud5pr Ped (Frozen)"
# d10x.combined$age_id<-as.factor(d10x.combined$age_id)
# #order by % cell passing QC
# d10x.combined$age_id<-factor(d10x.combined$age_id, c("C88_caud3pr Adult", "C97_caud3pr Adult","C70_caud5pr Adult",
#                                                      "C68_caud3pr Adult","C82_caud3pr Adult","C54_caud3pr Adult","C39_caud3prNPC Adult", "C39_caud3prTLH Adult",
#                                                      "C85_caud3pr Ped","C96_caud3pr Ped","C93_caud3pr Ped", "C64_caud5pr Ped","C104_bx5pr Ped",
#                                                      "IFALD006 Ped", "IFALD073 Ped","IFALD030 Ped"))
# levels(d10x.combined$age_id)<-gsub(" ","\n", levels(d10x.combined$age_id))
# 
# individual_split<-DimPlot(d10x.combined, reduction = "umap", group.by = "CellType_rough", split.by="age_id",pt.size=0.25, ncol=5)+colscale_cellType+ggtitle("")
# save_plts(individual_split, "IFALD_individual_roughCell_facet_rPCA_UMAP", w=22,h=16)
# 
# cell_num_all<-as.data.frame(table(d10x.combined@meta.data$AgeGroup))
# colnames(cell_num_all)<-c("AgeGroup","CellCount")
# age_split<-DimPlot(d10x.combined, reduction = "umap", group.by = "CellType_rough", split.by="AgeGroup",pt.size=0.25)+colscale_cellType+ggtitle("")+
#   geom_text(aes(x=-9, y=-14, label=paste0("n = ",comma(CellCount))), cell_num_all)
# save_plts(age_split, "IFALD_age_roughCell_facet_rPCA_UMAP", w=10,h=5)
# 
# 
# 
# ########
# ## mean genes per cell
# ########
# d10x.combined@meta.data %>%
#   group_by(individual) %>%
#   summarise(max = mean(nFeature_RNA, na.rm=TRUE))
# 
# 
# ##############
# ## Really no hepatocytes in ped samples
# ##############
# ## diff in ped (none in ped?)
# hep_age<-FeaturePlot(d10x.combined, reduction = "umap", features = c("ALB", "CPS1", "CYP3A4",
#                                                                      "MGST1", "CYP2E1"),
#                      ncol = 2, split.by='AgeGroup')
# save_plts(hep_age, "IFALD_hep_markers_age_rPCA_UMAP", w=6,h=14)
# 
# 
# ##############
# ## RBC
# ##############
# d10x.combined_RBC<-subset(d10x.combined, subset = seurat_clusters %in% c(19))
# d10x.combined_RBC <- RunPCA(d10x.combined_RBC, npcs = 30, verbose = FALSE)
# d10x.combined_RBC <- RunUMAP(d10x.combined_RBC, reduction = "pca", dims = 1:30)
# d10x.combined_RBC <- FindNeighbors(d10x.combined_RBC, reduction = "pca", dims = 1:30)
# d10x.combined_RBC <- FindClusters(d10x.combined_RBC, resolution = 0.2)
# 
# RBC_overlap<-FeaturePlot(d10x.combined_RBC, features = c("HBB", "FCGR3A"), blend = TRUE)
# save_plts(RBC_overlap, "IFALD_RBC_receptor_overlap_rPCA_UMAP", w=14,h=4)
# 
# DimPlot(d10x.combined_RBC, reduction = "umap",group.by="individual", pt.size=1)
# FeaturePlot(d10x.combined_RBC, reduction = "umap",features="nFeature_RNA", pt.size=1, label=T)
# 
# 
# d10x.combined_RBC@meta.data$CellType_rough<-as.character(d10x.combined_RBC@meta.data$CellType_rough)
# d10x.combined_RBC@meta.data$CellType_rough[which(d10x.combined_RBC@meta.data$seurat_clusters%in%c("0","1","2"))]<-"Erythrocytes"
# d10x.combined_RBC@meta.data$CellType_rough[which(d10x.combined_RBC@meta.data$seurat_clusters%in%c("3"))]<-"Myeloid Erythrocytes (phagocytosis)"
# 
# 
# DimPlot(d10x.combined_RBC, reduction = "umap",group.by="CellType_rough", pt.size=1)+colscale_cellType+ggtitle("")+
#   annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(d10x.combined))))
# 
# 
# 
# ##############
# ## Myeloid clustering
# ##############
# d10x.combined_myeloid<-subset(d10x.combined, subset = CellType_rough %in% c("Myeloid cells"))
# d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)
# d10x.combined_myeloid <- RunUMAP(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
# d10x.combined_myeloid <- FindNeighbors(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
# d10x.combined_myeloid <- FindClusters(d10x.combined_myeloid, resolution = 0.2)
# 
# myeloid_cluster_umap<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=T)
# myeloid_cluster_umap
# save_plts(myeloid_cluster_umap, "IFALD_myeloid_cluster_umap", w=5,h=4)
# 
# myeloid_cluster_umap_individual<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=T, split.by = "age_id", ncol=5)
# myeloid_cluster_umap_individual
# save_plts(myeloid_cluster_umap_individual, "IFALD_myeloid_cluster_umap_individual", w=10,h=6)
# 
# 
# cell_cluster_count<-d10x.combined_myeloid@meta.data %>%  group_by(age_id, seurat_clusters) %>%
#   summarise(n = n()) %>%
#   mutate(freq = n / sum(n))
# 
# cell_cluster_count<-as.data.frame(cell_cluster_count)
# cell_cluster_count<-merge(cell_cluster_count, d10x.combined_myeloid@meta.data[,c("age_id","AgeGroup")], by="age_id")
# cell_cluster_count<-cell_cluster_count[!duplicated(cell_cluster_count),]
# 
# bar_individual<-ggplot(cell_cluster_count, aes(fill=seurat_clusters, y=n, x=age_id)) +
#   geom_bar(position="fill", stat="identity", color="black")+theme_bw()+th+
#   facet_wrap(~AgeGroup, scale="free_x")
# save_plts(bar_individual, "IFALD_bar_individual_myeloid", w=14,h=6)
# 
# ## low quality clusters?
# myeloid_cluster_umap_nf<-FeaturePlot(d10x.combined_myeloid, reduction = "umap", feature="nuclear_fraction", pt.size=0.25, label=T)
# myeloid_cluster_umap_nf
# save_plts(myeloid_cluster_umap_nf, "IFALD_myeloid_cluster_umap_nuclearfraction", w=5,h=4)
# 
# myeloid_cluster_umap_cellstatus<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, group.by="cell_status",label=T)+scale_color_manual(values=c("grey","red","cornflowerblue"),name="Nuclear\nFraction")
# myeloid_cluster_umap_cellstatus
# save_plts(myeloid_cluster_umap_cellstatus, "IFALD_myeloid_cluster_umap_cellstatus", w=6,h=4)
# 
# umap_MTmyeloid<-FeaturePlot(d10x.combined_myeloid, features = "percent.mt", min.cutoff = "q9", pt.size=1)
# umap_MTmyeloid
# save_plts(umap_MTmyeloid, "IFALD_umap_MTmyeloid", w=5,h=4)
# 
# cluster_MT<-d10x.combined_myeloid@meta.data %>%  group_by(seurat_clusters) %>%
#   dplyr::summarize(Mean = mean(percent.mt, na.rm=TRUE))
# cell_cluster_count<-as.data.frame(cluster_MT)
# 
# box_MT<-ggplot(d10x.combined_myeloid@meta.data, aes(seurat_clusters, percent.mt)) +
#   geom_violin(fill="lightgrey", color="lightgrey")+xlab("Myeloid Cluster")+ylab("Percent MT")+
#   geom_boxplot(width=0.1, outlier.size = 0.05)+theme_bw()+th
# box_MT
# save_plts(box_MT, "IFALD_box_MT", w=8,h=3)
# 
# umap_nfeaturemyeloid<-FeaturePlot(d10x.combined_myeloid, features = "nFeature_RNA", min.cutoff = "q9", pt.size=1)
# save_plts(umap_nfeaturemyeloid, "IFALD_umap_nfeaturemyeloid", w=5,h=4)
# 
# cluster_nFeature<-d10x.combined_myeloid@meta.data %>%  group_by(seurat_clusters) %>%
#   dplyr::summarize(Mean = mean(nFeature_RNA, na.rm=TRUE))
# cell_cluster_count<-as.data.frame(cluster_nFeature)
# 
# box_nfeature<-ggplot(d10x.combined_myeloid@meta.data, aes(seurat_clusters, nFeature_RNA)) +
#   geom_violin(fill="lightgrey", color="lightgrey")+xlab("Myeloid Cluster")+ylab("Number of Feature \n(nFeature_RNA)")+
#   geom_boxplot(width=0.1, outlier.size = 0.05)+theme_bw()+th
# box_nfeature
# save_plts(box_nfeature, "IFALD_box_nfeature", w=8,h=3)
# 
# 
# ############
# ## Myeloid labelling
# ############
# recent_recruit_myeloid<-c("S100A8","S100A9","CD68","LYZ")
# kuffer_signature<-c("VSIG4","MARCO","CD5L","HMOX1")
# neutro_gene<-c("CSF3R","FCGR3B","NAMPT","CXCR2","DEFA3","DEFA4")
# MHCII<-c("HLA-DRA","CD74","HLA-DPB1","HLA-DQB1")
# LSEC<-c("CALCRL","STAB2","FCN2","FCN3")
# 
# 
# FeaturePlot(d10x.combined_myeloid, features = recent_recruit_myeloid, min.cutoff = "q9", pt.size=0.25)
# FeaturePlot(d10x.combined_myeloid, features = kuffer_signature, min.cutoff = "q9", pt.size=0.25)
# FeaturePlot(d10x.combined_myeloid, features = neutro_gene, min.cutoff = "q9", pt.size=0.25)
# FeaturePlot(d10x.combined_myeloid, features = MHCII, min.cutoff = "q9", pt.size=0.25)
# FeaturePlot(d10x.combined_myeloid, features = LSEC, min.cutoff = "q9", pt.size=0.25)
# 
# 
# FeaturePlot(d10x.combined_myeloid, features = c("HBB","HBA2","FCGR3A","MARCO"), min.cutoff = "q9", pt.size=1)
# #https://www.frontiersin.org/articles/10.3389/fimmu.2019.02035/full
# FeaturePlot(d10x.combined_myeloid, features = c("CD14","CD5L","FCGR3A","MARCO"), min.cutoff = "q9", pt.size=0.25)
# 
# FeaturePlot(d10x.combined_myeloid, features = c("LYZ", "S100A8", "CD14", "S100A10", "HLA-DRA", "CD74", "IFI30", "HLA-DPB1", "SECISBP2L"), min.cutoff = "q9", pt.size=0.25)# SLAN: SECISBP2L
# 
# ## Dot plot of all markers
# DefaultAssay(d10x.combined_myeloid) <- "RNA"
# pdf(file = here("figures/IFALD_dot_plots_myeloid.pdf"), w=10, h=10)
# DotPlot(object = d10x.combined_myeloid, features = B_genes)+xlab("B Cell Marker")
# DotPlot(object = d10x.combined_myeloid, features = T_genes)+xlab("T Cell Marker")
# DotPlot(object = d10x.combined_myeloid, features = NK_genes)+xlab("NK Cell Marker")
# DotPlot(object = d10x.combined_myeloid, features = LEC_genes)+xlab("LSEC Marker")
# DotPlot(object = d10x.combined_myeloid, features = Hepatocyte_genes[1:15])+xlab("Hepatocyte Marker")
# DotPlot(object = d10x.combined_myeloid, features = Hepatocyte_genes[16:29])+xlab("Hepatocyte Marker")
# DotPlot(object = d10x.combined_myeloid, features = Cholangiocytes_genes)+xlab("Cholangiocyte Marker")
# DotPlot(object = d10x.combined_myeloid, features = HSCs_genes[1:14])+xlab("HSC Marker")
# DotPlot(object = d10x.combined_myeloid, features = HSCs_genes[15:27])+xlab("HSC Marker")
# DotPlot(object = d10x.combined_myeloid, features = Macrophage_genes)+xlab("Macrophage Marker")
# dev.off()
# DefaultAssay(d10x.combined_myeloid) <- "integrated"
# 
# ## one unclear cluster
# cluster4.markers <-  FindMarkers(d10x.combined_myeloid, ident.1 = 6,ident.2 = 2, min.pct = 0.25)
# head(cluster4.markers, n = 10)
# head(cluster4.markers[which(cluster4.markers$avg_log2FC>0),], n = 20)
# FeaturePlot(d10x.combined_myeloid, features = c("IDO1","HLA-DOB","C1orf54","CLEC9A"), min.cutoff = "q9", pt.size=1)
# FeaturePlot(d10x.combined_myeloid, features = c("IDO1","CD83","FOXP1","CLEC9A"), min.cutoff = "q9", pt.size=1)
# # https://www.frontiersin.org/articles/10.3389/fimmu.2022.1006501/full
# #https://pesquisa.bvsalud.org/global-literature-on-novel-coronavirus-2019-ncov/resource/fr/covidwho-992926
# 
# d10x.combined_myeloid@meta.data$CellType_rough<-as.character(d10x.combined_myeloid@meta.data$CellType_rough)
# d10x.combined_myeloid@meta.data$CellType_rough[which(d10x.combined_myeloid@meta.data$seurat_clusters%in%c("0"))]<-"RR Myeloid"
# d10x.combined_myeloid@meta.data$CellType_rough[which(d10x.combined_myeloid@meta.data$seurat_clusters%in%c("1","3"))]<-"KC Like"
# d10x.combined_myeloid@meta.data$CellType_rough[which(d10x.combined_myeloid@meta.data$seurat_clusters%in%c("4"))]<-"Neutrophil"
# d10x.combined_myeloid@meta.data$CellType_rough[which(d10x.combined_myeloid@meta.data$seurat_clusters%in%c("5"))]<-"Doublet"
# d10x.combined_myeloid@meta.data$CellType_rough[which(d10x.combined_myeloid@meta.data$seurat_clusters%in%c("2"))]<-"Macrophage (MHCII high)"
# d10x.combined_myeloid@meta.data$CellType_rough[which(d10x.combined_myeloid@meta.data$seurat_clusters%in%c("6"))]<-"Macrophage (CLEC9A high)"
# 
# DefaultAssay(d10x.combined_myeloid) <- "RNA"
# marker_mye<-DotPlot(object = d10x.combined_myeloid, features = c(recent_recruit_myeloid, kuffer_signature, neutro_gene, MHCII, LSEC, "CD3D", "IL32"))
# save_plts(marker_mye, "IFALD_myeloid_dotplot_doublets", w=20,h=5)
# DefaultAssay(d10x.combined_myeloid) <- "integrated"
# 
# myeloid_cluster_umap<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=T, group.by = "CellType_rough")+colscale_cellType
# myeloid_cluster_umap
# save_plts(myeloid_cluster_umap, "IFALD_myeloid_cluster_umap_labelled", w=7,h=5)
# 
# plt_entropy_individual<-entropy_d10(d10x.combined_myeloid, "individual")
# entropy_individual<-entropy_plt(plt_entropy_individual, "individual", d10x.combined_myeloid)
# save_plts(entropy_individual, "IFALD_entropy_individual_myeloid", w=15,h=10)
# 
# plt_entropy_age<-entropy_d10(d10x.combined_myeloid, "AgeGroup")
# entropy_age<-entropy_plt(plt_entropy_age, "AgeGroup", d10x.combined_myeloid)
# save_plts(entropy_age, "IFALD_entropy_age_myeloid", w=15,h=10)
# 
# 
# 
# 
# 
# ############
# ## B cell labelling
# ############
# d10x.combined_bcell<-subset(d10x.combined, subset = CellType_rough %in% c("B-cells"))
# d10x.combined_bcell <- RunPCA(d10x.combined_bcell, npcs = 30, verbose = FALSE)
# d10x.combined_bcell <- RunUMAP(d10x.combined_bcell, reduction = "pca", dims = 1:30)
# d10x.combined_bcell <- FindNeighbors(d10x.combined_bcell, reduction = "pca", dims = 1:30)
# d10x.combined_bcell <- FindClusters(d10x.combined_bcell, resolution = 0.3)
# 
# bcell_cluster_umap<-DimPlot(d10x.combined_bcell, reduction = "umap", pt.size=0.25, label=T)
# bcell_cluster_umap
# save_plts(bcell_cluster_umap, "IFALD_bcell_cluster_umap", w=5,h=4)
# 
# b_genes_noIG<-c("POU2F2", "FCER2", "MS4A1", "LTB","CD37","CD79B", "CD19")
# immunoglobins<-c("IGKC","IGHG1","IGLC2", "IGHA1")
# 
# FeaturePlot(d10x.combined_bcell, features = b_genes_noIG, min.cutoff = "q9", pt.size=1)
# FeaturePlot(d10x.combined_bcell, features = immunoglobins, min.cutoff = "q9", pt.size=1)
# 
# 
# bcell_cluster_umap_cellstatus<-DimPlot(d10x.combined_bcell, reduction = "umap", pt.size=0.25, group.by="cell_status",label=T)+scale_color_manual(values=c("grey","red","cornflowerblue"),name="Nuclear\nFraction")
# bcell_cluster_umap_cellstatus
# save_plts(bcell_cluster_umap_cellstatus, "IFALD_bcell_cluster_umap_cellstatus", w=6,h=4)
# umap_MTbcell<-FeaturePlot(d10x.combined_bcell, features = "percent.mt", min.cutoff = "q9", pt.size=1)
# save_plts(umap_MTbcell, "IFALD_umap_MTbcell", w=5,h=4)
# 
# phase_bcell<-DimPlot(d10x.combined_bcell, reduction = "umap", group.by = "Phase" )+scale_color_manual(values=c("#e6ab02","#386cb0","#1b9e77"))
# save_plts(phase_bcell, "IFALD_umap_phase_bcell", w=5,h=4)
# 
# phase_marker<-FeaturePlot(d10x.combined_bcell, reduction = "umap", features = c("MKI67","TOP2A"), ncol = 2)
# save_plts(phase_marker, "IFALD_umap_phaseMarker_bcell", w=10,h=4)
# 
# ## one unclear cluster
# # platlets
# #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7049341/
# cluster4.markers <-  FindMarkers(d10x.combined_bcell, ident.1 = 9, min.pct = 0.25)
# head(cluster4.markers, n = 10)
# platelets_markers<-FeaturePlot(d10x.combined_bcell, reduction = "umap", features = c("PPBP","NRGN"), ncol = 2)
# save_plts(platelets_markers, "IFALD_platelets_markers_cluster_umap", w=8,h=4)
# 
# FeaturePlot(d10x.combined_bcell, reduction = "umap", features = c("PPBP","NRGN","ITGA2B","CSF3R","CD61"), ncol = 2)
# 
# 
# 
# ## Dot plot of all markers
# DefaultAssay(d10x.combined_bcell) <- "RNA"
# pdf(file = here("figures/IFALD_dot_plots_bcell.pdf"), w=10, h=10)
# DotPlot(object = d10x.combined_bcell, features = B_genes)+xlab("B Cell Marker")
# DotPlot(object = d10x.combined_bcell, features = T_genes)+xlab("T Cell Marker")
# DotPlot(object = d10x.combined_bcell, features = NK_genes)+xlab("NK Cell Marker")
# DotPlot(object = d10x.combined_bcell, features = LEC_genes)+xlab("LSEC Marker")
# DotPlot(object = d10x.combined_bcell, features = Hepatocyte_genes[1:15])+xlab("Hepatocyte Marker")
# DotPlot(object = d10x.combined_bcell, features = Hepatocyte_genes[16:29])+xlab("Hepatocyte Marker")
# DotPlot(object = d10x.combined_bcell, features = Cholangiocytes_genes)+xlab("Cholangiocyte Marker")
# DotPlot(object = d10x.combined_bcell, features = HSCs_genes[1:14])+xlab("HSC Marker")
# DotPlot(object = d10x.combined_bcell, features = HSCs_genes[15:27])+xlab("HSC Marker")
# DotPlot(object = d10x.combined_bcell, features = Macrophage_genes)+xlab("Macrophage Marker")
# dev.off()
# DefaultAssay(d10x.combined_bcell) <- "integrated"
# 
# 
# d10x.combined_bcell@meta.data$CellType_rough<-as.character(d10x.combined_bcell@meta.data$CellType_rough)
# d10x.combined_bcell@meta.data$CellType_rough[which(d10x.combined_bcell@meta.data$seurat_clusters%in%c("0","7","8"))]<-"Mature B-cells"
# d10x.combined_bcell@meta.data$CellType_rough[which(d10x.combined_bcell@meta.data$seurat_clusters%in%c("1","2","5"))]<-"Plasma cells"
# d10x.combined_bcell@meta.data$CellType_rough[which(d10x.combined_bcell@meta.data$seurat_clusters%in%c("6","4"))]<-"Cycling T-cells"
# d10x.combined_bcell@meta.data$CellType_rough[which(d10x.combined_bcell@meta.data$seurat_clusters%in%c("3"))]<-"Cycling Myeloid"
# d10x.combined_bcell@meta.data$CellType_rough[which(d10x.combined_bcell@meta.data$seurat_clusters%in%c("10"))]<-"Doublet"
# d10x.combined_bcell@meta.data$CellType_rough[which(d10x.combined_bcell@meta.data$seurat_clusters%in%c("9"))]<-"Platelets"
# 
# bcell_cluster_umap<-DimPlot(d10x.combined_bcell, reduction = "umap", pt.size=0.25, label=T, group.by = "CellType_rough")+colscale_cellType
# bcell_cluster_umap
# save_plts(bcell_cluster_umap, "IFALD_BCell_cluster_umap_labelled", w=6,h=4)
# 
# plt_entropy_individual<-entropy_d10(d10x.combined_bcell, "individual")
# entropy_individual<-entropy_plt(plt_entropy_individual, "individual", d10x.combined_bcell)
# save_plts(entropy_individual, "IFALD_entropy_individual_bcell", w=15,h=10)
# 
# plt_entropy_Treatment<-entropy_d10(d10x.combined_bcell, "Treatment")
# entropy_Treatment<-entropy_plt(plt_entropy_Treatment, "Treatment", d10x.combined_bcell)
# save_plts(entropy_Treatment, "IFALD_entropy_treatment_bcell", w=15,h=10)
# 
# 
# 
# ###########
# ## Sub cluster NK T
# ###########
# d10x.combined_NK_T<-subset(d10x.combined, subset = CellType_rough %in% c("NK and T cells"))
# d10x.combined_NK_T <- RunPCA(d10x.combined_NK_T, npcs = 30, verbose = FALSE)
# d10x.combined_NK_T <- RunUMAP(d10x.combined_NK_T, reduction = "pca", dims = 1:30)
# d10x.combined_NK_T <- FindNeighbors(d10x.combined_NK_T, reduction = "pca", dims = 1:30)
# d10x.combined_NK_T <- FindClusters(d10x.combined_NK_T, resolution = 0.2)
# 
# NK_T_umap<-DimPlot(d10x.combined_NK_T, label=T)
# NK_T_umap
# save_plts(NK_T_umap, "IFALD_NK_T_umap", w=8,h=6)
# 
# 
# tcell_cluster_umap_cellstatus<-DimPlot(d10x.combined_NK_T, reduction = "umap", pt.size=0.25, group.by="cell_status",label=T)+scale_color_manual(values=c("grey","red","cornflowerblue"),name="Nuclear\nFraction")
# tcell_cluster_umap_cellstatus
# save_plts(tcell_cluster_umap_cellstatus, "IFALD_tcell_cluster_umap_cellstatus", w=6,h=4)
# umap_MTtcell<-FeaturePlot(d10x.combined_NK_T, features = "percent.mt", min.cutoff = "q9", pt.size=1)
# save_plts(umap_MTtcell, "IFALD_umap_MTtcell", w=5,h=4)
# 
# FeaturePlot(d10x.combined_NK_T, reduction = "umap", features = NK_genes, ncol = 2)
# FeaturePlot(d10x.combined_NK_T, reduction = "umap", features = gd_genes)
# FeaturePlot(d10x.combined_NK_T, reduction = "umap", features = T_genes, ncol = 2)
# FeaturePlot(d10x.combined_NK_T, reduction = "umap", features = neutro_gene, ncol = 2)
# 
# ## Dot plot of all markers
# DefaultAssay(d10x.combined_NK_T) <- "RNA"
# pdf(file = here("figures/IFALD_dot_plots_tcell.pdf"), w=10, h=10)
# DotPlot(object = d10x.combined_NK_T, features = B_genes)+xlab("B Cell Marker")
# DotPlot(object = d10x.combined_NK_T, features = T_genes)+xlab("T Cell Marker")
# DotPlot(object = d10x.combined_NK_T, features = NK_genes)+xlab("NK Cell Marker")
# DotPlot(object = d10x.combined_NK_T, features = LEC_genes)+xlab("LSEC Marker")
# DotPlot(object = d10x.combined_NK_T, features = Hepatocyte_genes[1:15])+xlab("Hepatocyte Marker")
# DotPlot(object = d10x.combined_NK_T, features = Hepatocyte_genes[16:29])+xlab("Hepatocyte Marker")
# DotPlot(object = d10x.combined_NK_T, features = Cholangiocytes_genes)+xlab("Cholangiocyte Marker")
# DotPlot(object = d10x.combined_NK_T, features = HSCs_genes[1:14])+xlab("HSC Marker")
# DotPlot(object = d10x.combined_NK_T, features = HSCs_genes[15:27])+xlab("HSC Marker")
# DotPlot(object = d10x.combined_NK_T, features = Macrophage_genes)+xlab("Macrophage Marker")
# dev.off()
# DefaultAssay(d10x.combined_NK_T) <- "integrated"
# 
# 
# ## one unclear cluster
# cluster4.markers <-  FindMarkers(d10x.combined_NK_T, ident.1 = 5, min.pct = 0.25)
# head(cluster4.markers, n = 20)
# FeaturePlot(d10x.combined_NK_T, reduction = "umap", features = c("CD7","CLNK","IL18"), ncol = 2)
# cluster4.markers <-  FindMarkers(d10x.combined_NK_T, ident.1 = 7, min.pct = 0.25)
# head(cluster4.markers, n = 20)
# FeaturePlot(d10x.combined_NK_T, reduction = "umap", features = c("IGHA1","IGHG4"), ncol = 2)
# 
# 
# d10x.combined_NK_T@meta.data$CellType_rough<-as.character(d10x.combined_NK_T@meta.data$CellType_rough)
# d10x.combined_NK_T@meta.data$CellType_rough[which(d10x.combined_NK_T@meta.data$seurat_clusters%in%c("0"))]<-"NKT cells"
# d10x.combined_NK_T@meta.data$CellType_rough[which(d10x.combined_NK_T@meta.data$seurat_clusters%in%c("1","2","4"))]<-"CD3+ T-cells"
# d10x.combined_NK_T@meta.data$CellType_rough[which(d10x.combined_NK_T@meta.data$seurat_clusters%in%c("3"))]<-"gd T-cells"
# d10x.combined_NK_T@meta.data$CellType_rough[which(d10x.combined_NK_T@meta.data$seurat_clusters%in%c("6","7"))]<-"Doublet"
# d10x.combined_NK_T@meta.data$CellType_rough[which(d10x.combined_NK_T@meta.data$seurat_clusters%in%c("4"))]<-"Cycling T-cells"
# d10x.combined_NK_T@meta.data$CellType_rough[which(d10x.combined_NK_T@meta.data$seurat_clusters%in%c("5"))]<-"CLNK T-cells"
# 
# 
# 
# NKT_cluster_umap<-DimPlot(d10x.combined_NK_T, reduction = "umap", pt.size=0.25, label=T, group.by = "CellType_rough")+colscale_cellType
# NKT_cluster_umap
# save_plts(NKT_cluster_umap, "IFALD_NKT_cluster_umap_labelled", w=6,h=4)
# 
# 
# plt_entropy_individual<-entropy_d10(d10x.combined_NK_T, "individual")
# entropy_individual<-entropy_plt(plt_entropy_individual, "individual", d10x.combined_NK_T)
# save_plts(entropy_individual, "IFALD_entropy_individual_tcell", w=15,h=10)
# 
# plt_entropy_age<-entropy_d10(d10x.combined_NK_T, "AgeGroup")
# entropy_age<-entropy_plt(plt_entropy_age, "AgeGroup", d10x.combined_NK_T)
# save_plts(entropy_age, "IFALD_entropy_age_tcell", w=15,h=10)
# 
# 
# 
# 
# ##############
# ## Relabel subtypes
# ##############
# 
# d10x.combined@meta.data$CellType_refined<-sapply(1:nrow(d10x.combined@meta.data), function(x){
#   if(rownames(d10x.combined@meta.data)[x]%in%rownames(d10x.combined_myeloid@meta.data)){
#     d10x.combined_myeloid@meta.data$CellType_rough[which(rownames(d10x.combined_myeloid@meta.data)==rownames(d10x.combined@meta.data)[x])]
#   }else{
#     if(rownames(d10x.combined@meta.data)[x]%in%rownames(d10x.combined_bcell@meta.data)){
#       d10x.combined_bcell@meta.data$CellType_rough[which(rownames(d10x.combined_bcell@meta.data)==rownames(d10x.combined@meta.data)[x])]
#     }else{
#       if(rownames(d10x.combined@meta.data)[x]%in%rownames(d10x.combined_NK_T@meta.data)){
#         d10x.combined_NK_T@meta.data$CellType_rough[which(rownames(d10x.combined_NK_T@meta.data)==rownames(d10x.combined@meta.data)[x])]
#       }else{
#         if(rownames(d10x.combined@meta.data)[x]%in%rownames(d10x.combined_RBC@meta.data)){
#           d10x.combined_RBC@meta.data$CellType_rough[which(rownames(d10x.combined_RBC@meta.data)==rownames(d10x.combined@meta.data)[x])]
#           }else{d10x.combined@meta.data$CellType_rough[x]}
#       }}}})
# 
# d10x.combined@meta.data$CellType_refined<-as.factor(d10x.combined@meta.data$CellType_refined)
# levels(d10x.combined@meta.data$CellType_refined)<-c("CD3+ T-cells","Cholangiocytes","CLNK T-cells","Cycling Myeloid", "Cycling T-cells",
#                                                     "Doublet","Erythrocytes","gd T-cells",
#                                                     "Hepatocytes","HSC","KC Like",
#                                                     "Low Quality","LSEC","Macrophage\n(CLEC9A high)",
#                                                     "Macrophage\n(MHCII high)", "Mature B-cells",#"Mast cell",
#                                                     "Myeloid Erythrocytes\n(phagocytosis)", "Neutrophil","NK-like cells",
#                                                     "Plasma cells","Platelets","RR Myeloid")
# 
# d10x.combined@meta.data$CellType_refined<-factor(d10x.combined@meta.data$CellType_refined, levels = c("Mature B-cells","Plasma cells",
#                                                                                                       "CD3+ T-cells","gd T-cells","NK-like cells","CLNK T-cells","Cycling T-cells",
#                                                                                                       "Cholangiocytes","LSEC",
#                                                                                                       "RR Myeloid","Macrophage\n(MHCII high)","KC Like","Macrophage\n(CLEC9A high)","Cycling Myeloid",
#                                                                                                       "Neutrophil",#"Mast cell",
#                                                                                                       "Myeloid Erythrocytes\n(phagocytosis)","Erythrocytes","Platelets",
#                                                                                                       "HSC","Hepatocytes","Doublet","Low Quality"))
# 
# all_refined_cluster_umap<-DimPlot(d10x.combined, reduction = "umap", pt.size=0.25, label=T,label.size = 3, group.by = "CellType_refined")+
#   colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
#   annotate("text",x=-13, y=-13, label=paste0("n = ",comma(ncol(d10x.combined))))
# all_refined_cluster_umap
# save_plts(all_refined_cluster_umap, "IFALD_refined_cellType_map", w=12,h=8)
# 
# 
# all_refined_cluster_umap_nolab<-DimPlot(d10x.combined, reduction = "umap", pt.size=0.25, group.by = "CellType_refined")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
#   annotate("text",x=-13, y=-13, label=paste0("n = ",comma(ncol(d10x.combined))))
# all_refined_cluster_umap_nolab
# save_plts(all_refined_cluster_umap_nolab, "IFALD_refined_cellType_map_nolabel", w=12,h=8)
# 
# individual_split<-DimPlot(d10x.combined, reduction = "umap", group.by = "CellType_refined", split.by="age_id",pt.size=0.25, ncol=4)+colscale_cellType+ggtitle("")
# save_plts(individual_split, "IFALD_individual_roughCell_facet_rPCA_UMAP_refined", w=22,h=14)
# 
# cell_num_all<-as.data.frame(table(d10x.combined@meta.data$AgeGroup))
# colnames(cell_num_all)<-c("AgeGroup","CellCount")
# age_split<-DimPlot(d10x.combined, reduction = "umap", group.by = "CellType_refined", split.by="AgeGroup",pt.size=0.25)+colscale_cellType+ggtitle("")+
#   geom_text(aes(x=-13, y=-13, label=paste0("n = ",comma(CellCount))), cell_num_all)
# save_plts(age_split, "IFALD_age_roughCell_facet_rPCA_UMAP_refined", w=15,h=7)
# 
# 
# d10x.combined$age_condition<-paste(d10x.combined$AgeGroup, d10x.combined$Treatment, sep=" ")
# cell_num_all<-as.data.frame(table(d10x.combined@meta.data$age_condition))
# colnames(cell_num_all)<-c("age_condition","CellCount")
# 
# d10x.combined$age_condition<-factor(d10x.combined$age_condition, levels=c( "Ped Healthy","Ped IFALD", "Adult Healthy"))
# 
# age_condition_split<-DimPlot(d10x.combined, reduction = "umap", group.by = "CellType_refined", split.by="age_condition", pt.size=0.05, ncol=2)+colscale_cellType+ggtitle("")+
#   geom_text(aes(x=-13, y=-13, label=paste0("n = ",comma(CellCount))), cell_num_all)+xlab("UMAP 1")+ylab("UMAP 2")
# save_plts(age_condition_split, "IFALD_age_condition_roughCell_facet_rPCA_UMAP_refined", w=15,h=10)
# 
# # 
# # d10x.combined_health<-subset(d10x.combined, subset = age_condition %in% c("Ped Healthy","Adult Healthy"))
# # d10x.combined_health$age_condition<-as.character(d10x.combined_health$age_condition)
# # cell_num_all_health<-cell_num_all[which(cell_num_all$age_condition%in% c("Ped Healthy","Adult Healthy")),]
# # age_condition_split<-DimPlot(d10x.combined_health, reduction = "umap", group.by = "CellType_refined", split.by="age_condition", pt.size=0.05, ncol=2)+colscale_cellType+ggtitle("")+
# #   geom_text(aes(x=-13, y=-13, label=paste0("n = ",comma(CellCount))), cell_num_all_health)+xlab("UMAP 1")+ylab("UMAP 2")
# # save_plts(age_condition_split, "IFALD_age_condition_roughCell_facet_rPCA_UMAP_refined_healthyonly", w=15,h=5.5)
# 
# 
# 
# 
# ##############
# ### entropy in clusters
# ##############
# plt_entropy_individual<-entropy_d10(d10x.combined, "individual")
# plt_entropy_age<-entropy_d10(d10x.combined, "AgeGroup")
# plt_entropy_chem<-entropy_d10(d10x.combined, "chemistry")
# plt_entropy_treatment<-entropy_d10(d10x.combined, "age_condition")
# 
# entropy_individual<-entropy_plt(plt_entropy_individual, "individual", d10x.combined)
# save_plts(entropy_individual, "IFALD_entropy_individual_allclusters", w=15,h=10)
# 
# entropy_age<-entropy_plt(plt_entropy_age, "AgeGroup", d10x.combined)
# save_plts(entropy_age, "IFALD_entropy_age_allclusters", w=15,h=10)
# 
# entropy_chem<-entropy_plt(plt_entropy_chem, "chemistry", d10x.combined)
# save_plts(entropy_chem, "IFALD_entropy_chemistry_allclusters", w=15,h=10)
# 
# entropy_Treatment<-entropy_plt(plt_entropy_treatment, "age_condition", d10x.combined)
# save_plts(entropy_Treatment, "IFALD_entropy_Treatment_allclusters", w=15,h=10)
# 
# 
# 
# ##############
# ## Save integrated with refined cluster labels
# ##############
# save(d10x.combined, file=paste(here("data/"),"IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds", sep=""))
# cell_label<-d10x.combined@meta.data
# save(cell_label, file=paste(here("data/"),"IFALD_adult_ped_cellRefined_withDropletQC.rds", sep=""))
# 





print(sessionInfo())




