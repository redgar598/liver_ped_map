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


#cellranger_out_loc <- here("/media/redgar/Seagate Portable Drive/ped_liver_map_raw")
#cellranger_out_realigned <- here("/media/redgar/Seagate Portable Drive/realign_samples")
cellranger_out_loc <- here("../../../projects/macparland/RE/PediatricAdult/ped_liver_map_raw")
cellranger_out_realigned <- here("../../../projects/macparland/RE/PediatricAdult/realign_samples")

samples<-c(list.files(cellranger_out_loc),list.files(cellranger_out_realigned))
print(samples)

meta<-read.table(here("data/data_transfer_updated_feb12_2024_IFALD_PBMC.csv"), header=T, sep=",")

samples<-samples[which(samples%in%meta$file)]
samples<-samples[-grep("PBMC",samples)]


d10x.list <- sapply(1:length(samples), function(y){
  if(length(grep("realign",samples[y]))==1){
    dataset_loc<-cellranger_out_realigned
  }else{
    dataset_loc<-cellranger_out_loc
  }

  caud<-meta$Sample_ID[which(meta$file == samples[y])]
  print(caud)
  print(file.path(dataset_loc,paste(samples[y],"/outs", sep=""),"filtered_feature_bc_matrix"))
  d10x <- Read10X(file.path(dataset_loc,paste(samples[y],"/outs", sep=""),"filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),caud,sep="-")
  # print(dim(d10x))
  #' Initialize the Seurat object with the raw (non-normalized data).
  d10x<-CreateSeuratObject(counts = d10x, project = "ped_adult_map", min.cells = 0, min.features = 0)

  #add meta data to each seurat object
  meta_cell<-data.frame(cell=colnames(d10x), individual=caud)
  meta_cell_add<-merge(meta_cell, meta, by.x="individual", by.y="Sample_ID")
  meta_cell_add$relALBChange<-NA
  meta_cell_add<-meta_cell_add[match(colnames(d10x), meta_cell_add$cell),]
  print(identical(meta_cell_add$cell, colnames(d10x)))
  rownames(meta_cell_add)<-meta_cell_add$cell
  d10x<- AddMetaData(d10x, meta_cell_add)

})

d10x.list



#######################################################################################################################################

d10x.list<-lapply(1:length(d10x.list), function(x){
  sobj<-d10x.list[[x]]
  
  # This script filters cells based on MALAT1 expression
  # Code is written assuming the use of a Seurat object, but
  # should be able to be applied to any single-cell object
  
  # Look at histogram of MALAT1 counts
  hist(sobj@assays$RNA@counts["MALAT1",], freq = FALSE, breaks=100)
  
  # If there is a peak (even a small peak) at zero, followed by a dip, then a more
  # normal distribution, run the following code.
  
  # Choose rough x value that surrounds the minimum between 0 and where slope starts to rise again; this is max_counts
  # Counts is vector of MALAT1 counts, min and max are range of count values to look between
  # Lower bw if minimum doesn't look totally accurate
  # Change lwd to change thickness of red line
  # Change breaks to change bins of histogram
  define_malat1_threshold <- function(counts, max_counts, min_counts = 0, bw = 0.1, lwd = 2, breaks = 100) {
    pdf(here("figures",paste(unique(sobj$individual),"MALAT1.pdf", sep="")))
    # Visualise MALAT1 histogram
    hist(counts, breaks = breaks, freq = FALSE, main=unique(sobj$individual))
    # Calculate the density values
    density_data <- density(counts, bw = bw)
    # Isolate x value at minimum
    
    ######## Updated this line as in one sample the first values in counts was negative so threw off the indexing
    threshold <- round(density_data$x[density_data$x >= min_counts & density_data$x <= max_counts][which(density_data$y[density_data$x >= min_counts & density_data$x <= max_counts] == min(density_data$y[density_data$x >= min_counts & density_data$x <= max_counts]))])
   
    # Visualise on histogram
    abline(v = threshold, col = "red", lwd = lwd)
    dev.off() 
    return(threshold)
  }
  
  # Run this function on MALAT1 reads, eg:
  threshold <- define_malat1_threshold(sobj@assays$RNA@counts["MALAT1",], max_counts = 100)
  
  ######## Manually decided which samples had a peak and a dip
  MALAT1_filter<-c("C85_caud3pr","C96_caud3pr","C39_caud3prNPC","IFALD073","IFALD030" ,"C105_caud5pr","C115")
  
  if(unique(sobj$individual)%in%MALAT1_filter){
    # Use this value to subset your cells
    cells_subset <- WhichCells(sobj, expression = MALAT1 > threshold, slot = "counts")
    sobj_subset <- subset(sobj, cells = cells_subset)
    
    print(paste(ncol(sobj), "cell originally", ncol(sobj_subset), "after MALAT1 filter"))
    
    sobj_subset
  }else{
    sobj
  }
  

})


#######################################################################################################################################


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
save(plt_QC_data, file=here("data","MALAT1_QC_metrics.Rdata"))

#load(here("data","IFALD_QC_metrics.Rdata"))




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

save_plts(cell_count, "MALAT1_QC_cellcount_age", w=8,h=4)


d10x <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "adult_ped_map")#add.cell.ids = alldata_names2,

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
## Integration
###############
#https://satijalab.org/seurat/articles/integration_rpca.html
print("RUNNING INTEGRATION")

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
save_plts(key_markers, "MALAT1_markers_rPCA_UMAP", w=25,h=20)

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

genes<-unique(c(Macrophage_genes, NK_T_genes, B_genes,LEC_genes,Hepatocyte_genes,Cholangiocytes_genes,HSCs_genes))

d10x.exp<-as.data.frame(d10x.combined[["RNA"]]@data)
d10x.exp.GOI<-d10x.exp[genes,]
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


## second option to plot (to maybe manually relable hepatocytes)
not_hep_cell<-colnames(cell_rough)[which(!(colnames(cell_rough)%in%c("Hepatocyte","CellType_rough" ,"seurat_clusters")))]

cell_rough$second_best_cell<-sapply(1:nrow(cell_rough), function(x){
  if(cell_rough$CellType_rough[x]!="Hepatocyte"){cell_rough$CellType_rough[x]}else{
    not_hep_mean_max<-max(cell_rough[x, not_hep_cell])
    paste("Hep_",not_hep_cell[which(cell_rough[x, not_hep_cell]==not_hep_mean_max)], sep="")}
})

meta<-d10x.combined@meta.data
meta$cell<-rownames(meta)
plt_summary<-merge(meta, cell_rough[,c("seurat_clusters","CellType_rough","second_best_cell")], by="seurat_clusters")
plt_summary<-plt_summary[match(rownames(d10x.combined@meta.data),plt_summary$cell),]
identical(plt_summary$cell, rownames(d10x.combined@meta.data))

rownames(plt_summary)<-plt_summary$cell

d10x.combined<- AddMetaData(d10x.combined, plt_summary)
d10x.combined
table(d10x.combined$CellType_rough)

table(d10x.combined$second_best_cell, d10x.combined$CellType_rough)



##############
## Save integrated to look local
##############
#save(d10x.combined, file=paste(here("data/"),"IFALD_adult_ped_integrated.rds", sep=""))
saveRDS(d10x.combined, file = here("../../../projects/macparland/RE/PediatricAdult/processed_data","MALAT1_adult_ped_integrated.rds"))

cell_label<-d10x.combined@meta.data
#save(cell_label, file=paste(here("data/"),"IFALD_adult_ped_cellRough.rds", sep=""))
saveRDS(cell_label, file = here("../../../projects/macparland/RE/PediatricAdult/processed_data","MALAT1_adult_ped_cellRough.rds"))






######################
## Compare filtered cells to original
######################

d10x.dropletQC<-readRDS(file = here("/media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024","IFALD_adult_ped_cellRough.rds"))
d10x.dropletQC$cell2<-sapply(1:nrow(d10x.dropletQC), function(x) paste(strsplit(d10x.dropletQC$cell[x], "-")[[1]][1], "-", d10x.dropletQC$individual[x], sep="" ))

load(here("/media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024","IFALD_QC_metrics.Rdata"))
plt_QC_data_noMALAT1<-plt_QC_data
plt_QC_data_noMALAT1$cell2<-sapply(1:nrow(plt_QC_data_noMALAT1), function(x) paste(strsplit(plt_QC_data_noMALAT1$cell[x], "-")[[1]][1], "-", plt_QC_data_noMALAT1$individual[x], sep="" ))

load(here("/media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024","MALAT1_QC_metrics.Rdata"))
plt_QC_data_MALAT1<-plt_QC_data


dim(plt_QC_data_MALAT1)
dim(plt_QC_data_noMALAT1)

length(which(plt_QC_data_noMALAT1$cell2 %in% plt_QC_data_MALAT1$cell))



plt_QC_data_noMALAT1$MALAT1<-"Passed"
plt_QC_data_noMALAT1$MALAT1[which(!(plt_QC_data_noMALAT1$cell2 %in% plt_QC_data_MALAT1$cell))]<-"Filtered"
table(plt_QC_data_noMALAT1$MALAT1)


plt_QC_data_noMALAT1$MT_nFeature_filter<-"Passed"
plt_QC_data_noMALAT1$MT_nFeature_filter[which(!(plt_QC_data_noMALAT1$cell2 %in% d10x.dropletQC$cell2))]<-"Filtered"
table(plt_QC_data_noMALAT1$MT_nFeature_filter)

plt_QC_data_noMALAT1$Both_filtering<-sapply(1:nrow(plt_QC_data_noMALAT1), function(x){
  if(plt_QC_data_noMALAT1$MT_nFeature_filter[x]=="Passed" & plt_QC_data_noMALAT1$MALAT1[x]=="Passed"){"Passed\nin Both"}else{
    if(plt_QC_data_noMALAT1$MT_nFeature_filter[x]=="Filtered" & plt_QC_data_noMALAT1$MALAT1[x]=="Passed"){"Filtered by\nMT and nFeature\nfilter only"}else{
      if(plt_QC_data_noMALAT1$MT_nFeature_filter[x]=="Passed" & plt_QC_data_noMALAT1$MALAT1[x]=="Filtered"){"Filtered by\nMALAT1 only"}else{"Filtered\nby both"}
    }
  }
})

table(plt_QC_data_noMALAT1$Both_filtering)



QC_plot<-lapply(c("nCount_RNA","nFeature_RNA","percent.mt","relALBChange","nuclear_fraction"), function(y){
  ggplot(plt_QC_data_noMALAT1, aes_string("Both_filtering", y)) + geom_violin(fill="lightgrey", color="lightgrey")+geom_boxplot(width=0.1, outlier.shape = NA) +
    theme_bw()+th_present+xlab("")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
})

plot_grid(plotlist=QC_plot, ncol=2, axis="h")
save_plts(plot_grid(plotlist=QC_plot, ncol=2, axis="h"), "MALAT1_filtering_comparison", h=20, w=10)


table(plt_QC_data_noMALAT1$Both_filtering, plt_QC_data_noMALAT1$cell_status)






whodis<-d10x.dropletQC[which(d10x.dropletQC$cell2 %in% plt_QC_data_noMALAT1$cell2[which(plt_QC_data_noMALAT1$Both_filtering=="Filtered by\nMALAT1 only")]),]

round(table(whodis$CellType_rough)/nrow(whodis)*100, 2)
round(table(d10x.dropletQC$CellType_rough)/nrow(d10x.dropletQC)*100, 2)


round((table(whodis$CellType_rough)/table(d10x.dropletQC$CellType_rough))*100, 2)
