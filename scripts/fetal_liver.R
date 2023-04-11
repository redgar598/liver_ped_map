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

source("scripts/00_pretty_plots.R")
source("scripts/00_entropy_d10x.R")


#############
#https://developmental.cellatlas.io/fetal-liver
#E-MTAB-7407
# downloaded the cellranger outputs
#############

                      # fetal_liver_h5ad <- read_h5ad(here("/media/redgar/Seagate Portable Drive/fetal_liver/download.h5ad"))
                      # # write_h5ad(fetal_liver_h5ad, here("/media/redgar/Seagate Portable Drive/fetal_liver/fetal_liver.h5ad"))
                      # write_csvs(fetal_liver_h5ad, here("/media/redgar/Seagate Portable Drive/fetal_liver/fetal_liver.csv"))



## Meta
#meta<-read.csv(here("/media/redgar/Seagate Portable Drive/fetal_liver/E-MTAB-7407.sdrf.txt"), sep="\t")
meta<-read.csv(here("/cluster/projects/macparland/RE/PediatricAdult/fetal_liver/E-MTAB-7407.sdrf.txt"), sep="\t")

meta_liver<-meta[which(meta$Factor.Value.organism.part.=="liver"),]

#dataset_loc <- here("/media/redgar/Seagate Portable Drive/fetal_liver")
dataset_loc <- here("/cluster/projects/macparland/RE/PediatricAdult/fetal_liver")

samples<-list.files(dataset_loc)
samples<-samples[-grep("MTAB|h5ad|liver",samples)]
samples<-samples[which(samples%in%meta_liver$Extract.Name)]
print(samples)


d10x.list <- sapply(1:length(samples), function(y){
  print(file.path(dataset_loc,paste(samples[y],"/GRCh38", sep="")))
  d10x <- Read10X(file.path(dataset_loc,paste(samples[y],"/GRCh38", sep="")))
        #colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),caud,sep="-")
        # print(dim(d10x))
  #' Initialize the Seurat object with the raw (non-normalized data).
  d10x<-CreateSeuratObject(counts = d10x, project = "fetal_liver", min.cells = 0, min.features = 0)

  # if(dir.exists(file.path(dataset_loc,paste(samples[y],"/GRCh38/raw_feature_bc_matrix", sep="")))) {
  #   ## SoupX needs clusters so quickly make clusters for each sample
  #   d10x    <- SCTransform(d10x, verbose = F)
  #   d10x    <- RunPCA(d10x, verbose = F)
  #   d10x    <- RunUMAP(d10x, dims = 1:30, verbose = F)
  #   d10x    <- FindNeighbors(d10x, dims = 1:30, verbose = F)
  #   d10x    <- FindClusters(d10x, verbose = T)
  #   meta_clusters    <- d10x@meta.data
  # 
  #   sc = load10X(file.path(dataset_loc,paste(samples[y],"/outs", sep="")))
  #   sc = setClusters(sc, setNames(meta_clusters$seurat_clusters, rownames(meta_clusters)))
  # 
  #   ######
  #   ## Load data and estimate soup profile
  #   ######
  #   # Estimate rho
  #   sc = autoEstCont(sc, forceAccept=TRUE)
  #   #Genes with highest expression in background. These are often enriched for ribosomal proteins.
  #   print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20))
  #   print(unique(sc$metaData$rho))
  #   # Clean the data
  #   out = adjustCounts(sc)
  # 
  #   ## Save a metric of soupness (ALB change after soupX)
  #   DR = sc$metaData[,sc$DR]
  #   df = DR
  #   old = colSums(sc$toc["ALB",rownames(df),drop=FALSE])
  #   new = colSums(out["ALB",rownames(df),drop=FALSE])
  #   relChange = (old-new)/old
  #   df$old = old
  #   df$new = new
  #   df$relALBChange=relChange
  #   df$cell<-rownames(df)
  # 
  #   ## make seurat object of adjusted counts
  #   d10x = CreateSeuratObject(out)
  # 
  #   #add meta data to each seurat object
  #   meta_cell<-data.frame(cell=colnames(d10x), Extract.Name=caud)
  #   meta_cell_add<-merge(meta_cell, meta, by.x="Extract.Name", by.y="Sample_ID")
  #   meta_cell_add<-merge(meta_cell_add, df[,c("cell","relALBChange")], by="cell")
  #   meta_cell_add<-meta_cell_add[match(colnames(d10x), meta_cell_add$cell),]
  #   print(identical(meta_cell_add$cell, colnames(d10x)))
  #   rownames(meta_cell_add)<-meta_cell_add$cell
  #   d10x<- AddMetaData(d10x, meta_cell_add)
  # }else{
    #add meta data to each seurat object
    meta_cell<-data.frame(cell=colnames(d10x), Extract.Name=samples[y])
    meta_cell_add<-merge(meta_cell, meta_liver[,c("Characteristics.age.","Characteristics.sex.","Characteristics.individual.","Characteristics.clinical.information.","Characteristics.facs.sorting.","Extract.Name")],
                                               by="Extract.Name")
    meta_cell_add$Barcodes<-sapply(1:nrow(meta_cell_add), function(x) strsplit(meta_cell_add$cell[x],"-")[[1]][1])
    
    
    ## Cell annotation 
    cell_labels<-read.csv(file.path(dataset_loc,paste(samples[y],"/",samples[y],".csv", sep="")))
    
    meta_cell_add<-merge(meta_cell_add, cell_labels, by="Barcodes")
    meta_cell_add<-meta_cell_add[match(colnames(d10x), meta_cell_add$cell),]
    print(identical(meta_cell_add$cell, colnames(d10x)))
    rownames(meta_cell_add)<-meta_cell_add$cell
    d10x<- AddMetaData(d10x, meta_cell_add)
  #}

  # ######
  # ## dropletQC
  # ######
  # if(file.exists(file.path(dataset_loc,paste(samples[y],"/outs", sep=""),"possorted_genome_bam.bam"))){
  #   nf1 <- nuclear_fraction_tags(
  #     outs = file.path(dataset_loc,paste(samples[y],"/outs", sep="")),
  #     tiles = 1, cores = 1, verbose = FALSE)
  #   head(nf1)
  # 
  #   print(identical(rownames(nf1), colnames(d10x)))
  #   d10x<- AddMetaData(d10x, nf1)
  #   d10x
  # 
  #   nf.umi <- data.frame(nf=d10x$nuclear_fraction,
  #                        umi=d10x$nCount_RNA)
  # 
  #   # Run identify_empty_drops
  #   empty_drop <- identify_empty_drops(nf_umi=nf.umi)
  #   empty_drop$Extract.Name<-d10x$Extract.Name
  #   empty_drop_damagedcell <- identify_damaged_cells(empty_drop, verbose = FALSE, output_plots = F)
  # 
  #   head(empty_drop_damagedcell[[1]])
  #   table(empty_drop_damagedcell[[1]]$cell_status)
  # 
  #   print(identical(rownames(empty_drop_damagedcell[[1]]), colnames(d10x)))
  #   d10x<- AddMetaData(d10x, empty_drop_damagedcell[[1]])
  #   d10x$nf<-NULL
  #   d10x$umi<-NULL
  #   print(d10x)
  #   d10x}else{
  #     d10x$nuclear_fraction<-NA
  #     d10x$cell_status<-NA
  #     print(d10x)
  #     d10x}
})

d10x.list



## cell counts
plt_count_raw<-lapply(1:length(d10x.list), function(x) {
  df<-data.frame(raw_cell_count=nrow(d10x.list[[x]]@meta.data),Extract.Name=unique(d10x.list[[x]]@meta.data$Extract.Name))
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
save_fetal_plts(qc_plts, "intital_QC_plts", w=6,h=4)

qc_plts_Extract.Name<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
  geom_point() + facet_wrap(~Extract.Name)+
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th
save_fetal_plts(qc_plts_Extract.Name, "intital_QC_plts_Extract.Name", w=12,h=4)

MT_plt<-ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
  geom_vline(xintercept = 25)+ theme_bw()+xlab("Percent Mitochondrial")+th
save_fetal_plts(MT_plt, "percentMT_plt", w=6,h=4)

MT_plt_Extract.Name<-ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
  facet_wrap(~Extract.Name, scales="free_y")+
  geom_vline(xintercept = 25)+ theme_bw()+xlab("Percent Mitochondrial")+th
save_fetal_plts(MT_plt_Extract.Name, "percentMT_plt_Extract.Name", w=8,h=4)



# nf_plts<-ggplot(plt_QC_data, aes(nuclear_fraction,log10(nCount_RNA),colour=cell_status)) +
#   geom_point() +  ylab("Number of Total Molecules\n(log 10 nCount) ")+xlab("Nuclear Fraction")+theme_bw()+th+
#   facet_wrap(~Extract.Name)
# nf_plts
# save_fetal_plts(nf_plts, "nuclear_fraction", w=6,h=4)

# ggplot(plt_QC_data, aes(nuclear_fraction,percent.mt,colour=cell_status)) +
#   geom_point() +  ylab("Number of Total Molecules\n(log 10 nCount) ")+xlab("Nuclear Fraction")+theme_bw()+th+
#   facet_wrap(~Extract.Name)+
#   scale_color_manual(values=c("grey","red","cornflowerblue"),name="Nuclear\nFraction")


#'We filter cells that have unique feature counts over 6,000 or less than 500
#'We filter cells that have >10% mitochondrial counts
#'we will also filter doublets as called by scrublet
invisible(lapply(1:length(d10x.list), function(x){
  d10x.list[[x]] <<- subset(d10x.list[[x]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 25)
}))

d10x.list

## cell counts after QC
plt_count_QC<-lapply(1:length(d10x.list), function(x) {
  df<-data.frame(qc_cell_count=nrow(d10x.list[[x]]@meta.data),Extract.Name=unique(d10x.list[[x]]@meta.data$Extract.Name))
  df})
plt_count_QC<-do.call(rbind, plt_count_QC)
print(plt_count_QC)

counts<-merge(plt_count_raw, plt_count_QC, by="Extract.Name")
meta<-merge(meta_liver,counts,by="Extract.Name")

cell_count<-grid.arrange(ggplot(meta, aes(Characteristics.age., raw_cell_count,fill=Characteristics.age.))+
                           geom_boxplot()+geom_point()+
                           theme_bw()+geom_text(aes(label=Extract.Name), hjust=-0.25, size=3)+xlab("Age Group")+
                           ylab("Total Cell Number")+th+ylim(0,5000)+
                           theme(legend.position = "none")+ggtitle("Before Quality Control"),
                         ggplot(meta, aes(Characteristics.age., qc_cell_count,fill=Characteristics.age.))+
                           geom_boxplot()+geom_point()+
                           theme_bw()+geom_text(aes(label=Extract.Name), hjust=-0.25, size=3)+xlab("Age Group")+
                           ylab("Total Cell Number")+th+ylim(0,5000)+
                           theme(legend.position = "none")+ggtitle("After Quality Control"), ncol=2)

save_fetal_plts(cell_count, "QC_cellcount_age", w=8,h=4)


d10x <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "adult_ped_map")#add.cell.ids = alldata_names2,

d10x

saveRDS(d10x, file = here("data","d10x_fetal_raw.rds"))

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
          # save_fetal_plts(pca_cellcycle, "pca_cellcycle", w=6,h=4)
          # 
          # pca_nfeature<-FeaturePlot(d10x, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
          # save_fetal_plts(pca_nfeature, "pca_nfeature", w=6,h=4)
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
          # save_fetal_plts(SCT_cluster_umap, "SCT_cluster_umap", w=6,h=4)
          # 
          # SCT_cluster_tsne<-DimPlot(d10x, reduction = "tsne", pt.size=0.25, label=T)
          # save_fetal_plts(SCT_cluster_tsne, "SCT_cluster_tsne", w=6,h=4)
          # 
          # cell_pca_SCT<-DimPlot(d10x, reduction="pca", group.by="Phase")
          # save_fetal_plts(cell_pca_SCT, "cell_PCA_afterSCT", w=6,h=4)
          # 
          # nFeature_UMAP_SCT<-FeaturePlot(d10x, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
          # save_fetal_plts(nFeature_UMAP_SCT, "nfeature_UMAP_afterSCT", w=6,h=4)
          # 
          # age_umap_sct<-DimPlot(d10x, reduction = "umap", group.by = "Characteristics.age.", pt.size=0.25)+fillscale_age
          # save_fetal_plts(age_umap_sct, "age_SCT_umap", w=6,h=4)
          # 
          # Extract.Name_umap_sct<-DimPlot(d10x, reduction = "umap", group.by = "Extract.Name", pt.size=1)
          # save_fetal_plts(Extract.Name_umap_sct, "Extract.Name_SCT_UMAP", w=6,h=4)
          # 
          # Cell.Labels_umap_sct<-DimPlot(d10x, reduction = "umap", group.by = "Cell.Labels", pt.size=0.5)
          # save_fetal_plts(Cell.Labels_umap_sct, "celltype_SCT_UMAP", w=10,h=6)
          # 
          # 
          # ## Low quality cluster?
          # MT_umap_SCT<-FeaturePlot(d10x, features = "percent.mt", min.cutoff = "q9", pt.size=1)
          # save_fetal_plts(pca_nfeature, "pca_nfeature", w=6,h=4)
          # 
          # ncount_umap_SCT<-FeaturePlot(d10x, features = "nCount_RNA", min.cutoff = "q9", pt.size=1)
          # save_fetal_plts(ncount_umap_SCT, "ncount_umap_SCT", w=6,h=4)
          # 
          # nfeature_umap_SCT<-FeaturePlot(d10x, features = "nFeature_RNA", min.cutoff = "q9", pt.size=1)
          # save_fetal_plts(nfeature_umap_SCT, "nfeature_umap_SCT", w=6,h=4)
          # 



                                #
                                # ###############
                                # ## Integration
                                # ###############
                                # #https://satijalab.org/seurat/articles/integration_rpca.html
                                # print("RUNNING INTEGRATION")
                                # ## Start back with raw data
                                # d10x <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "adult_ped_map")#add.cell.ids = alldata_names2,
                                # d10x
                                #
                                # # split the dataset into a list of two seurat objects (3' and 5')
                                # d10x.list.chem <- SplitObject(d10x, split.by = "chemistry")
                                #
                                #
                                # # normalize, identify variable features and score cell cycle for each dataset independently
                                # s.genes <- cc.genes$s.genes
                                # g2m.genes <- cc.genes$g2m.genes
                                #
                                # d10x.list.chem <- lapply(X = d10x.list.chem, FUN = function(x) {
                                #   x <- NormalizeData(x)
                                #   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
                                #   x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
                                # })
                                #
                                # # select features that are repeatedly variable across datasets for integration run PCA on each
                                # # dataset using these features
                                # features <- SelectIntegrationFeatures(object.list = d10x.list.chem)
                                # d10x.list.chem <- lapply(X = d10x.list.chem, FUN = function(x) {
                                #   #x <- ScaleData(x, features = features, verbose = FALSE)
                                #   x <- ScaleData(x, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), features = features, verbose = FALSE)
                                #   x <- RunPCA(x, features = features, verbose = FALSE)
                                # })
                                #
                                #
                                #
                                # ## Identify anchors
                                # chem.anchors <- FindIntegrationAnchors(object.list = d10x.list.chem, anchor.features = features, reduction = "rpca")
                                # d10x.combined <- IntegrateData(anchorset = chem.anchors)
                                #
                                # DefaultAssay(d10x.combined) <- "integrated"
                                #
                                # print("INTEGRATED")


###############
## Integration
###############
#https://satijalab.org/seurat/articles/integration_rpca.html
print("RUNNING INTEGRATION")

## run integration across age
d10x <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "adult_ped_map")#add.cell.ids = alldata_names2,
d10x

# split the dataset into a list of two seurat objects (3' and 5')
d10x.list.age <- SplitObject(d10x, split.by = "Characteristics.age.")

# normalize, identify variable features and score cell cycle for each dataset independently
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

d10x.list.age <- lapply(X = d10x.list.age, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = d10x.list.age)
d10x.list.age <- lapply(X = d10x.list.age, FUN = function(x) {
  #x <- ScaleData(x, features = features, verbose = FALSE)
  x <- ScaleData(x, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})



## Identify anchors
age.anchors <- FindIntegrationAnchors(object.list = d10x.list.age, anchor.features = features, reduction = "rpca")
d10x.combined <- IntegrateData(anchorset = age.anchors)

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
save_fetal_plts(SCT_cluster_umap, "rPCA_cluster_umap", w=6,h=4)

age_umap_sct<-DimPlot(d10x.combined, reduction = "umap", group.by = "Characteristics.age.", pt.size=0.25)+fillscale_age
save_fetal_plts(age_umap_sct, "age_rPCA_umap", w=6,h=4)

MT_umap_sct<-FeaturePlot(d10x.combined, reduction = "umap", features = "percent.mt", pt.size=0.25)
save_fetal_plts(MT_umap_sct, "MT_rPCA_umap", w=5,h=4)

Extract.Name_umap_sct<-DimPlot(d10x.combined, reduction = "umap", group.by = "Extract.Name", pt.size=1)
save_fetal_plts(Extract.Name_umap_sct, "Extract.Name_rPCA_UMAP", w=6,h=4)

Cell.Labels_umap_sct<-DimPlot(d10x.combined, reduction = "umap", group.by = "Cell.Labels", pt.size=0.5)
save_fetal_plts(Cell.Labels_umap_sct, "celltype_rPCA_UMAP", w=12,h=6)






##############
## Save integrated with refined cluster labels
##############
save(d10x.combined, file=paste(here("data/"),"fetal_integrated.rds", sep=""))
cell_label<-d10x.combined@meta.data
save(cell_label, file=paste(here("data/"),"fetal_celllabels.rds", sep=""))






print(sessionInfo())




