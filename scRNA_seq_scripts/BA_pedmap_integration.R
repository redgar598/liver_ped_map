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
library(cowplot)
library(RColorBrewer)



#########################################################
## BA
#########################################################
#' 
#' BA_mat_loc <- here("/cluster/projects/macparland/RE/PediatricAdult/BA_GSE176189")
#' 
#' samples<-c(list.files(BA_mat_loc))
#' print(samples)
#' 
#' meta<-read.table(here("/cluster/projects/macparland/RE/PediatricAdult/BA_GSE176189/GSE176189_series_matrix.csv"), header=T, sep=",")
#' meta<-t(meta)
#' colnames(meta)<-meta[1,]
#' meta<-meta[-1,]
#' meta<-data.frame(meta,  quote=FALSE)
#' colnames(meta)<-gsub("X.","",colnames(meta))
#' meta<-meta[,c("Sample_geo_accession","Sample_characteristics_ch1")]
#' meta$Sample_name<-rownames(meta)
#' 
#' meta<-meta[grep("CELL",meta$Sample_name),]
#' 
#' 
#' d10x.list <- sapply(1:nrow(meta), function(y){
#'   
#'   sampl_file<-samples[grep(meta$Sample_geo_accession[y], samples)]
#'   print(meta$Sample_name[y])
#'   
#'   if(meta$Sample_name[y]=="E_CELL"){
#'     dir_filt<-"YXC18"}else{
#'     dir_filt<-meta$Sample_name[y]
#'     }
#'   
#'   sample_path<-file.path(BA_mat_loc,paste(sampl_file,"/",strsplit(dir_filt,"_")[[1]][1],"_","filtered_feature_bc_matrix", sep=""))
#'   d10x <- Read10X(sample_path)
#'   # print(dim(d10x))
#'   #' Initialize the Seurat object with the raw (non-normalized data).
#'   d10x<-CreateSeuratObject(counts = d10x, project = "BA_GSE176189", min.cells = 0, min.features = 0)
#' 
#'   #add meta data to each seurat object
#'   meta_cell<-data.frame(cell=colnames(d10x), individual=meta$Sample_name[y])
#'   meta_cell_add<-merge(meta_cell, meta, by.x="individual", by.y="Sample_name")
#'   meta_cell_add<-meta_cell_add[match(colnames(d10x), meta_cell_add$cell),]
#'   print(identical(meta_cell_add$cell, colnames(d10x)))
#'   rownames(meta_cell_add)<-meta_cell_add$cell
#'   d10x<- AddMetaData(d10x, meta_cell_add)
#'   }
#' )
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
#' 
#' 
#' 
#' qc_plts<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
#'   geom_point() +
#'   scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
#'   geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
#'   geom_hline(yintercept = 6000) +theme_bw()+th
#' 
#' qc_plts_chem<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
#'   geom_point() + facet_wrap(~Sample_characteristics_ch1)+
#'   scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
#'   geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
#'   geom_hline(yintercept = 6000) +theme_bw()+th
#' 
#' qc_plts_individual<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
#'   geom_point() + facet_wrap(~individual)+
#'   scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
#'   geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
#'   geom_hline(yintercept = 6000) +theme_bw()+th
#' 
#' MT_plt<-ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
#'   geom_vline(xintercept = 25)+ theme_bw()+xlab("Percent Mitochondrial")+th
#' 
#' MT_plt_individual<-ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
#'   facet_wrap(~individual, scales="free_y")+
#'   geom_vline(xintercept = 25)+ theme_bw()+xlab("Percent Mitochondrial")+th
#' 
#' 
#' 
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
#' meta<-merge(meta,counts,by.x="Sample_name", by.y="individual")
#' 
#' cell_count<-grid.arrange(ggplot(meta, aes(Sample_characteristics_ch1, raw_cell_count,fill=Sample_characteristics_ch1))+
#'                            geom_boxplot()+geom_point()+
#'                            theme_bw()+geom_text(aes(label=Sample_name), hjust=-0.25, size=3)+xlab("Diagnosis")+
#'                            ylab("Total Cell Number")+th+fillscale_age+ylim(0,22000)+
#'                            theme(legend.position = "none")+ggtitle("Before Quality Control"),
#'                          ggplot(meta, aes(Sample_characteristics_ch1, qc_cell_count,fill=Sample_characteristics_ch1))+
#'                            geom_boxplot()+geom_point()+
#'                            theme_bw()+geom_text(aes(label=Sample_name), hjust=-0.25, size=3)+xlab("Diagnosis")+
#'                            ylab("Total Cell Number")+th+fillscale_age+ylim(0,22000)+
#'                            theme(legend.position = "none")+ggtitle("After Quality Control"), ncol=2)
#' 
#' save_plts(cell_count, "BA_QC_cellcount", w=8,h=4)
#' 
#' 
#' d10x <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "BA_GSE176189")#add.cell.ids = alldata_names2,
#' 
#' d10x
#' 
#' ##saveRDS(d10x, file = here("data","BA_d10x_adult_ped_raw.rds"))
#' saveRDS(d10x, file = here("../../../projects/macparland/RE/PediatricAdult/BA_GSE176189","BA_d10x_adult_ped_raw.rds"))
#' 
##################################################################
## Taylor
#########################################################
taylor_mat_loc <- here("/cluster/projects/macparland/RE/PediatricAdult/Taylor_GSE163650")

samples<-c(list.files(taylor_mat_loc))
print(samples)

samples<-samples[grep("filtered",samples)]

d10x.list <- sapply(1:length(samples), function(y){
  
  print(samples[y])
  condition<-strsplit(samples[y],"_")[[1]][2]
  sample_id<-strsplit(samples[y],"_")[[1]][1]
  
  sample_path<-file.path(taylor_mat_loc,samples[y])
  d10x <-Read10X_h5(sample_path, use.names = TRUE, unique.features = TRUE)
  
  # print(dim(d10x))
  #' Initialize the Seurat object with the raw (non-normalized data).
  d10x<-CreateSeuratObject(counts = d10x, project = "taylor", min.cells = 0, min.features = 0)
  
  #add meta data to each seurat object
  meta_cell<-data.frame(cell=colnames(d10x), individual=sample_id, condition=condition)
  rownames(meta_cell)<-meta_cell$cell
  d10x<- AddMetaData(d10x, meta_cell)
}
)

d10x.list


# 
# ## cell counts
# plt_count_raw<-lapply(1:length(d10x.list), function(x) {
#   df<-data.frame(raw_cell_count=nrow(d10x.list[[x]]@meta.data),individual=unique(d10x.list[[x]]@meta.data$individual))
#   df})
# plt_count_raw<-do.call(rbind, plt_count_raw)
# print(plt_count_raw)
# 


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

#' 
#' # Visualize QC metrics
#' #nFeature number of unique genes
#' #nCount number of total molecules
#' plt_QC_data<-do.call(rbind, lapply(1:length(d10x.list), function(x) d10x.list[[x]]@meta.data))
#' 
#' 
#' 
#' qc_plts<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
#'   geom_point() +
#'   scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
#'   geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
#'   geom_hline(yintercept = 6000) +theme_bw()+th
#' 
#' qc_plts_chem<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
#'   geom_point() + facet_wrap(~Sample_characteristics_ch1)+
#'   scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
#'   geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
#'   geom_hline(yintercept = 6000) +theme_bw()+th
#' 
#' qc_plts_individual<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
#'   geom_point() + facet_wrap(~individual)+
#'   scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
#'   geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
#'   geom_hline(yintercept = 6000) +theme_bw()+th
#' 
#' MT_plt<-ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
#'   geom_vline(xintercept = 25)+ theme_bw()+xlab("Percent Mitochondrial")+th
#' 
#' MT_plt_individual<-ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
#'   facet_wrap(~individual, scales="free_y")+
#'   geom_vline(xintercept = 25)+ theme_bw()+xlab("Percent Mitochondrial")+th
#' 
#' 
#' #'We filter cells that have unique feature counts over 6,000 or less than 500
#' #'We filter cells that have >10% mitochondrial counts
#' #'we will also filter doublets as called by scrublet
#' d10x.list.raw<-d10x.list

invisible(lapply(1:length(d10x.list), function(x){
  d10x.list[[x]] <<- subset(d10x.list[[x]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 25)
}))

d10x.list
# 
# ## cell counts after QC
# plt_count_QC<-lapply(1:length(d10x.list), function(x) {
#   df<-data.frame(qc_cell_count=nrow(d10x.list[[x]]@meta.data),individual=unique(d10x.list[[x]]@meta.data$individual))
#   df})
# plt_count_QC<-do.call(rbind, plt_count_QC)
# print(plt_count_QC)
# 
# counts<-merge(plt_count_raw, plt_count_QC, by="individual")
# counts<-melt(counts)
# 
# cell_count<-ggplot(counts, aes(variable, value))+
#   geom_boxplot()+geom_point()+
#   theme_bw()+geom_text(aes(label=individual), hjust=-0.25, size=3)+xlab("Individual")+
#   ylab("Total Cell Number")+th+
#   theme(legend.position = "none")+ggtitle("Before Quality Control")
# #save_plts(cell_count, "Taylor_QC_cellcount", w=5,h=4)


d10x <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "BA_GSE176189")#add.cell.ids = alldata_names2,

d10x

##saveRDS(d10x, file = here("data","BA_d10x_adult_ped_raw.rds"))
saveRDS(d10x, file = here("../../../projects/macparland/RE/PediatricAdult/Taylor_GSE163650","Taylor_GSE163650_raw.rds"))

#########################################################
## Integration
#########################################################
d10x<-readRDS(file = here("../../../projects/macparland/RE/PediatricAdult/processed_data/","IFALD_d10x_adult_ped_raw.rds"))
load(here("../../../projects/macparland/RE/PediatricAdult/processed_data/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)
d10x <- AddMetaData(d10x, metadata = cell_label)
d10x$Sex[which(d10x$individual%in%c("C113","C115"))]<-"M"


d10x_BA<-readRDS(file = here("../../../projects/macparland/RE/PediatricAdult/BA_GSE176189","BA_d10x_adult_ped_raw.rds"))
d10x_taylor<-readRDS(file = here("../../../projects/macparland/RE/PediatricAdult/Taylor_GSE163650","Taylor_GSE163650_raw.rds"))

d10x_ped_BA_list <- list(Ped = d10x,BA = d10x_BA,taylor=d10x_taylor)
d10x_ped_BA_list

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

d10x_ped_BA_list <- lapply(X = d10x_ped_BA_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = d10x_ped_BA_list)
d10x_ped_BA_list <- lapply(X = d10x_ped_BA_list, FUN = function(x) {
  #x <- ScaleData(x, features = features, verbose = FALSE)
  x <- ScaleData(x, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})



## Identify anchors
chem.anchors <- FindIntegrationAnchors(object.list = d10x_ped_BA_list, anchor.features = features, reduction = "rpca")
d10x.combined <- IntegrateData(anchorset = chem.anchors)

DefaultAssay(d10x.combined) <- "integrated"

print("INTEGRATED")


# Run the standard workflow for visualization and clustering
d10x.combined <- ScaleData(d10x.combined, verbose = FALSE)
d10x.combined <- RunPCA(d10x.combined, npcs = 30, verbose = FALSE)
d10x.combined <- RunUMAP(d10x.combined, reduction = "pca", dims = 1:30)

d10x.combined <- FindNeighbors(d10x.combined, reduction = "pca", dims = 1:30)
d10x.combined <- FindClusters(d10x.combined, resolution = 0.5)

d10x.combined

saveRDS(d10x.combined, file = here("../../../projects/macparland/RE/PediatricAdult/BA_GSE176189","BA_Taylor_d10x_adult_ped_integrated.rds"))



###########################################################################################
#### Cell labelling
###########################################################################################

d10x_BA_taylor<-readRDS(here("/media/redgar/Seagate Portable Drive/BA_Taylor_d10x_adult_ped_integrated.rds"))

d10x_BA_taylor$condition[which(is.na(d10x_BA_taylor$condition))]<-d10x_BA_taylor$age_condition
d10x_BA_taylor$condition[which(is.na(d10x_BA_taylor$condition))]<-d10x_BA_taylor$Sample_characteristics_ch1

DimPlot(d10x_BA_taylor, group.by = "orig.ident")
DimPlot(d10x_BA_taylor, group.by = "CellType_refined")+colscale_cellType
