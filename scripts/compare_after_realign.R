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
source("scripts/00_plot_gene_exp.R")
source("scripts/00_fanciest_UMAP.R")




/media/redgar/Seagate Portable Drive/realign_samples/realign006_realign
/media/redgar/Seagate Portable Drive/ped_liver_map_raw/MacParland_Sonya__HSC-FI_006


d10x <- Read10X("/media/redgar/Seagate Portable Drive/realign_samples/realign006_realign/outs/filtered_feature_bc_matrix")
colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),"realigned",sep="-")
# print(dim(d10x))
#' Initialize the Seurat object with the raw (non-normalized data).
d10x<-CreateSeuratObject(counts = d10x, project = "realigned", min.cells = 0, min.features = 0)

## SoupX needs clusters so quickly make clusters for each sample
d10x    <- SCTransform(d10x, verbose = F)
d10x    <- RunPCA(d10x, verbose = F)
d10x    <- RunUMAP(d10x, dims = 1:30, verbose = F)
d10x    <- FindNeighbors(d10x, dims = 1:30, verbose = F)
d10x    <- FindClusters(d10x, verbose = T)
meta_clusters    <- d10x@meta.data

sc = load10X("/media/redgar/Seagate Portable Drive/realign_samples/realign006_realign/outs")
sc = setClusters(sc, setNames(meta_clusters$seurat_clusters, rownames(meta_clusters)))

## Load data and estimate soup profile
# Estimate rho
sc = autoEstCont(sc, forceAccept=TRUE)
#Genes with highest expression in background. These are often enriched for ribosomal proteins.
print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20))
print(unique(sc$metaData$rho))
# Clean the data
out = adjustCounts(sc)

## make seurat object of adjusted counts
d10x = CreateSeuratObject(out)

## dropletQC
nf1 <- nuclear_fraction_tags(
  outs = file.path("/media/redgar/Seagate Portable Drive/realign_samples/realign006_realign/outs"),
  tiles = 1, cores = 1, verbose = FALSE)
head(nf1)

print(identical(rownames(nf1), colnames(d10x)))
d10x<- AddMetaData(d10x, nf1)
d10x

nf.umi <- data.frame(nf=d10x$nuclear_fraction,
                     umi=d10x$nCount_RNA)

# Run identify_empty_drops
empty_drop <- identify_empty_drops(nf_umi=nf.umi)
empty_drop$version<-"realigned"
empty_drop_damagedcell <- identify_damaged_cells(empty_drop, verbose = FALSE, output_plots = F)

head(empty_drop_damagedcell[[1]])
table(empty_drop_damagedcell[[1]]$cell_status)

print(identical(rownames(empty_drop_damagedcell[[1]]), colnames(d10x)))
d10x<- AddMetaData(d10x, empty_drop_damagedcell[[1]])
d10x$nf<-NULL
d10x$umi<-NULL
d10x_realigned<-d10x


#####
## original
#####
d10x <- Read10X("/media/redgar/Seagate Portable Drive/ped_liver_map_raw/MacParland_Sonya__HSC-FI_006/outs/filtered_feature_bc_matrix")
colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),"original",sep="-")
# print(dim(d10x))
#' Initialize the Seurat object with the raw (non-normalized data).
d10x<-CreateSeuratObject(counts = d10x, project = "original", min.cells = 0, min.features = 0)

## SoupX needs clusters so quickly make clusters for each sample
d10x    <- SCTransform(d10x, verbose = F)
d10x    <- RunPCA(d10x, verbose = F)
d10x    <- RunUMAP(d10x, dims = 1:30, verbose = F)
d10x    <- FindNeighbors(d10x, dims = 1:30, verbose = F)
d10x    <- FindClusters(d10x, verbose = T)
meta_clusters    <- d10x@meta.data

sc = load10X("/media/redgar/Seagate Portable Drive/ped_liver_map_raw/MacParland_Sonya__HSC-FI_006/outs")
sc = setClusters(sc, setNames(meta_clusters$seurat_clusters, rownames(meta_clusters)))

## Load data and estimate soup profile
# Estimate rho
sc = autoEstCont(sc, forceAccept=TRUE)
#Genes with highest expression in background. These are often enriched for ribosomal proteins.
print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20))
print(unique(sc$metaData$rho))
# Clean the data
out = adjustCounts(sc)

## make seurat object of adjusted counts
d10x = CreateSeuratObject(out)

## dropletQC
nf1 <- nuclear_fraction_tags(
  outs = file.path("/media/redgar/Seagate Portable Drive/ped_liver_map_raw/MacParland_Sonya__HSC-FI_006/outs"),
  tiles = 1, cores = 1, verbose = FALSE)
head(nf1)

print(identical(rownames(nf1), colnames(d10x)))
d10x<- AddMetaData(d10x, nf1)
d10x

nf.umi <- data.frame(nf=d10x$nuclear_fraction,
                     umi=d10x$nCount_RNA)

# Run identify_empty_drops
empty_drop <- identify_empty_drops(nf_umi=nf.umi)
empty_drop$version<-"original"
empty_drop_damagedcell <- identify_damaged_cells(empty_drop, verbose = FALSE, output_plots = F)

head(empty_drop_damagedcell[[1]])
table(empty_drop_damagedcell[[1]]$cell_status)

print(identical(rownames(empty_drop_damagedcell[[1]]), colnames(d10x)))
d10x<- AddMetaData(d10x, empty_drop_damagedcell[[1]])
d10x$nf<-NULL
d10x$umi<-NULL
d10x_original<-d10x


d10x.list<-list(d10x_original, d10x_realigned)



## cell counts
plt_count_raw<-lapply(1:length(d10x.list), function(x) {
  df<-data.frame(raw_cell_count=nrow(d10x.list[[x]]@meta.data),version=unique(d10x.list[[x]]@meta.data$version))
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





qc_plts_version<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
  geom_point() + facet_wrap(~version)+
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th
save_plts(qc_plts_version, "realign_intital_QC_plts_version", w=8,h=4)

MT_plt<-ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
  geom_vline(xintercept = 25)+ theme_bw()+xlab("Percent Mitochondrial")+th+facet_wrap(~version)
save_plts(MT_plt, "realign_percentMT_plt", w=6,h=4)




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
  df<-data.frame(qc_cell_count=nrow(d10x.list[[x]]@meta.data),version=unique(d10x.list[[x]]@meta.data$version))
  df})
plt_count_QC<-do.call(rbind, plt_count_QC)
print(plt_count_QC)

counts<-merge(plt_count_raw, plt_count_QC, by="version")

cell_count<-grid.arrange(ggplot(counts, aes(version, raw_cell_count,fill=version))+
                           geom_boxplot()+geom_point()+
                           theme_bw()+xlab("Version")+
                           ylab("Total Cell Number")+th+ylim(0,60000)+
                           theme(legend.position = "none")+ggtitle("Before Quality Control"),
                         ggplot(counts, aes(version, qc_cell_count,fill=version))+
                           geom_boxplot()+geom_point()+
                           theme_bw()+xlab("Version")+
                           ylab("Total Cell Number")+th+ylim(0,60000)+
                           theme(legend.position = "none")+ggtitle("After Quality Control"), ncol=2)



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

DimPlot(d10x, group.by = "version")


######################
## cell cycle gene expression
######################
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

d10x <- CellCycleScoring(d10x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

pca_cellcycle<-DimPlot(d10x, reduction="pca",  group.by = "Phase", split.by = "version")

######################
### add cell type
######################
load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label<-cell_label[grep("IFALD006",cell_label$individual),]

cell_label1<-cell_label
rownames(cell_label1)<-gsub("_14","_1",rownames(cell_label1))
cell_label2<-cell_label
rownames(cell_label2)<-gsub("_14","_2",rownames(cell_label2))

cell_label<-rbind(cell_label2, cell_label1)

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
missing_in_old<-cell_label[which(is.na(cell_label$index)),]
missing_in_old$index<-colnames(d10x)[which(!(colnames(d10x)%in%cell_label$index))]
cell_label<-rbind(cell_label, missing_in_old)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]

identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

fanciest_UMAP(d10x,NA,F)
save_plts(fanciest_UMAP(d10x,NA,F), "realigned_IFALD006", w=6,h=5)

  umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x, reduction = "umap"))#
  umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
  meta_myeloid<-d10x@meta.data
  meta_myeloid$cell<-rownames(meta_myeloid)
  plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")
  
  len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
  len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
  arr <- list(x = min(plt_myeloid$UMAP_1), y = min(plt_myeloid$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)

  type<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
      geom_point(size = 0.06, colour= "black", stroke = 1)+
      geom_point(aes(color=version),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
      scale_color_manual(values=c("gold","cornflowerblue"))+theme_bw()+
      annotate("segment", 
               x = arr$x, xend = arr$x + c(arr$x_len, 0), 
               y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
               arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
      theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                         axis.title.x = element_text(size=5,hjust = 0.05),
                         axis.title.y = element_text(size=5,hjust = 0.05,angle = 90))+ guides(colour = guide_legend(override.aes = list(size=3)))
  save_plts(type, "realigned_IFALD006_version", w=6,h=5)
  

## lost in realign
  plt_myeloid$index_repeat<-gsub("_1|_2","",plt_myeloid$index)
  plt_myeloid_original<-plt_myeloid[which(plt_myeloid$version=="original"),]
  plt_myeloid_realigned<-plt_myeloid[which(plt_myeloid$version=="realigned"),]
  missing<-plt_myeloid_original[which(!(plt_myeloid_original$index_repeat%in%plt_myeloid_realigned$index_repeat)),]
  
  plt_myeloid$lost_in_realign<-"kept"
  plt_myeloid$lost_in_realign[(plt_myeloid$index%in%missing$index)]<-"lost"

 ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
    geom_point(aes(color=CellType_refined),size=0.05)+
    geom_point(data=plt_myeloid[which(plt_myeloid$lost_in_realign=="lost"),], size = 0.06, colour= "black", stroke = 1)+
    geom_point(aes(color=CellType_refined),data=plt_myeloid[which(plt_myeloid$lost_in_realign=="lost"),], size=0.05)+
    xlab("UMAP 1")+ylab("UMAP 2")+
    colscale_cellType+
    annotate("segment", 
             x = arr$x, xend = arr$x + c(arr$x_len, 0), 
             y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
             arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
    theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                       axis.title.x = element_text(size=5,hjust = 0.05),
                       axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                       legend.position = "none")
 
ggplot(plt_myeloid, aes(lost_in_realign, nCount_RNA))  + geom_violin()
ggplot(plt_myeloid, aes(lost_in_realign, nFeature_RNA))  + geom_violin()
ggplot(plt_myeloid, aes(lost_in_realign, nuclear_fraction))  + geom_violin()
table(plt_myeloid$lost_in_realign, plt_myeloid$cell_status)
ggplot(plt_myeloid, aes(lost_in_realign, relALBChange))  + geom_violin()
