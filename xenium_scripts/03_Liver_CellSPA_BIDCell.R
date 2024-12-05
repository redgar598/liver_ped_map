library(CellSPA)
library(ggplot2)
library(ggthemes)
library(scater)
library(SingleCellExperiment)
library(SpatialExperiment)
library(here)
library(dplyr)
library(RColorBrewer)
library(scales)
library(tiff)
library(sp)
library(rgeos)

library(ggsignif)
library(cowplot)

theme_set(theme_bw())

# #Map tools was archived
# url <- "https://cran.r-project.org/src/contrib/Archive/maptools/maptools_1.1-8.tar.gz"
# pkgFile <- "maptools_1.1-8.tar.gzz"
# download.file(url = url, destfile = pkgFile)
# # Install package
# install.packages(pkgs=pkgFile, type="source", repos=NULL)
# # Delete package tarball
# unlink(pkgFile)
# 
# # CellSPA install
# BiocManager::install("SydneyBioX/CellSPA")

## remotes::install_version("Matrix", version = "1.6-1.1")
#Then for seurat
#remotes::install_version("Matrix", version = "1.6-3")

source(here("scripts/00_Cell_SPA_functions.R"))
source(here("scripts/00_cellshape_metrics.R"))
source("scripts/00_pretty_plots.R")


## the panel comes with cell type labels
load(file=here("data/cell_type_labels_BIDCell.RData"))



count_files<-c("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_2/cell_gene_matrices/2024_03_21_16_22_10/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_3/cell_gene_matrices/2024_03_21_16_38_48/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_4/cell_gene_matrices/2024_03_21_17_03_05/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C105/cell_gene_matrices/2024_05_17_13_32_03/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C85/cell_gene_matrices/2024_05_17_13_29_40/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C95/cell_gene_matrices/2024_07_04_12_35_29/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C101/cell_gene_matrices/2024_06_28_09_23_40/expr_mat.csv")


samples<-c("C94_2", "C94_3","C94_4","C105","C85","C95","C101")
models<-c("2024_03_21_16_22_10","2024_03_21_16_38_48","2024_03_21_17_03_05","2024_05_17_13_32_03","2024_05_17_13_29_40","2024_07_04_12_35_29","2024_06_28_09_23_40")




lapply(1:length(count_files), function(x){
  
  data_dir<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x],"/cell_gene_matrices/", models[x], "/expr_mat.csv", sep="")
  tiff_path<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x], "/model_outputs/",models[x],"/test_output/epoch_1_step_4000_connected.tif", sep="")
  smpl<-samples[x]
  print(smpl)
  
  spe<-readBIDCell_adapted(data_dir, tiff_path = tiff_path)
  spe <- processingSPE(spe,  qc_range = list(total_transciprts = c(20, 2000), total_genes = c(20, Inf)))
          #spe <- CellSPA::subset(spe, 1:500)
  spe
  
  # Baseline metrics
  spe <- generatePolygon(spe)
  spe <- calBaselineAllMetrics(spe, verbose = TRUE)
  head(rowData(spe))
  head(colData(spe))
  
  cellSPA_metrics<-as.data.frame(colData(spe))
  
  
  ##### centroids
  tiff_res <- readTIFF(tiff_path, as.is = TRUE)
  tiff_res <- reshape2::melt(tiff_res)
  tiff_res <- tiff_res[tiff_res$value != 0, ]
  
  colnames(tiff_res) <- c("coord_y", "coord_x", "cell_id")
  tiff_res <- tiff_res[, c("coord_x", "coord_y", "cell_id")]
  tiff_res$cell_id <- as.numeric(tiff_res$cell_id)
  tiff_res$coord_x <- as.numeric(tiff_res$coord_x)
  tiff_res$coord_y <- as.numeric(tiff_res$coord_y)
  tiff_res <- tiff_res[order(tiff_res$cell_id), ]
  rownames(tiff_res) <- NULL
  
  centroids<-as.data.frame(
    tiff_res %>%
      group_by(cell_id) %>%
      summarise(centroid_x = sum(coord_x) / length(coord_x), centroid_y = sum(coord_y) / length(coord_y)))
  
  plt_umap_xenium_sample<-plt_umap_xenium[which(plt_umap_xenium$sample == smpl),]
  plt_umap_xenium_sample$cell<-sapply(1:nrow(plt_umap_xenium_sample), function(x) strsplit(plt_umap_xenium_sample$cell[x],"_")[[1]][1])
  
  cell_centroid<-merge(plt_umap_xenium_sample, centroids, by.x="cell",by.y="cell_id")
  
  metrics_all_samples<-merge(cell_centroid,cellSPA_metrics[,c("cell_id","cell_area", "elongation", "compactness", "eccentricity", "sphericity",  "solidity", "convexity", "circularity", "density")], by.x="cell", by.y="cell_id")
  
  save(metrics_all_samples, file=paste(here("data/"),smpl,"_centroid_cellSPA_metrics.RData",sep=""))
  })


##########
## load metrics
##########


metrics_list <- lapply(1:5, function(x){
  smpl<-samples[x]
  print(smpl)
  load(paste(here("data/"),smpl,"_centroid_cellSPA_metrics.RData",sep=""))
  metrics_all_samples})


metrics_all_samples<-do.call(rbind, metrics_list)

metrics_all_samples<-metrics_all_samples[which(metrics_all_samples$sample!="C94_2"),]
metrics_all_samples$CellType[which(metrics_all_samples$CellType=="LSEC (Periportal)")]<-"VEC"

size<-ggplot(metrics_all_samples, aes(reorder(CellType, cell_area), cell_area))+
  geom_hline(yintercept=median(metrics_all_samples$cell_area), color="grey30")+
  geom_violin(fill="grey",color="grey")+geom_boxplot(aes(fill=CellType),width=0.1, outlier.shape=NA)+
  fillscale_cellType+xlab("")+ylab("Cell Area")+theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
save_plts(size, "cell_area_cell_type", w=10, h=5)

size_sample<-ggplot(metrics_all_samples, aes(reorder(CellType, cell_area), cell_area))+
  geom_hline(yintercept=median(metrics_all_samples$cell_area), color="grey30")+
  geom_boxplot(aes(fill=sample), outlier.shape=NA,position = position_dodge(preserve = "single"))+
  xlab("")+ylab("Cell Area")+scale_fill_manual(values=c("#D64A56","#ac3c27","cornflowerblue","#6db0dd"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
save_plts(size_sample, "cell_area_cell_type_sample", w=10, h=5)


metrics_all_samples$elongation_directionless<-abs(metrics_all_samples$elongation-1)
ggplot(metrics_all_samples, aes(reorder(CellType, elongation_directionless), elongation_directionless))+
  geom_boxplot(aes(fill=sample), outlier.shape=NA,position = position_dodge(preserve = "single"))+
  xlab("")+ylab("Elongation")+scale_fill_manual(values=c("#D64A56","#ac3c27","cornflowerblue","#6db0dd"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))# one is a square

ggplot(metrics_all_samples, aes(reorder(CellType, elongation), elongation))+
  geom_violin(fill="grey",color="grey")+geom_boxplot(aes(fill=CellType),width=0.1, outlier.shape=NA)+
  fillscale_cellType+xlab("")+ylab("Elongation")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_hline(yintercept=1) # one is a square

ggplot(metrics_all_samples, aes(reorder(CellType, eccentricity), eccentricity))+
  geom_boxplot(aes(fill=sample), outlier.shape=NA,position = position_dodge(preserve = "single"))+
  xlab("")+ylab("Eccentricity")+scale_fill_manual(values=c("#D64A56","#ac3c27","cornflowerblue","#6db0dd"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))# one is a square

ggplot(metrics_all_samples, aes(reorder(CellType, eccentricity), eccentricity))+
  geom_violin(fill="grey",color="grey")+geom_boxplot(aes(fill=CellType),width=0.1, outlier.shape=NA)+
  fillscale_cellType+xlab("")+ylab("Cell Area")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_hline(yintercept=median(metrics_all_samples$eccentricity))


ggplot(metrics_all_samples, aes(reorder(CellType, density), density))+
  geom_violin(fill="grey",color="grey")+geom_boxplot(aes(fill=CellType),width=0.1, outlier.shape=NA)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("Density")+
  fillscale_cellType+
  geom_hline(yintercept=median(metrics_all_samples$density))


ggplot(metrics_all_samples, aes(reorder(CellType, density), density))+
  geom_boxplot(aes(fill=sample), outlier.shape=NA,position = position_dodge(preserve = "single"))+
  xlab("")+ylab("Density")+scale_fill_manual(values=c("#D64A56","#ac3c27","cornflowerblue","#6db0dd"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))# one is a square


ggplot(metrics_all_samples, aes(reorder(CellType, circularity), circularity))+
  geom_violin(fill="grey",color="grey")+geom_boxplot(aes(fill=CellType),width=0.1, outlier.shape=NA)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("Circularity")+
  fillscale_cellType

ggplot(metrics_all_samples, aes(reorder(CellType, compactness), compactness))+
  geom_violin(fill="grey",color="grey")+geom_boxplot(aes(fill=CellType),width=0.1, outlier.shape=NA)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("Compactness")+
  fillscale_cellType# one is a circle

ggplot(metrics_all_samples, aes(reorder(CellType, convexity), convexity))+
  geom_violin(fill="grey",color="grey")+geom_boxplot(aes(fill=CellType),width=0.1, outlier.shape=NA)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("convexity")+
  fillscale_cellType# one is a circle

solid<-ggplot(metrics_all_samples, aes(reorder(CellType, solidity), solidity))+
  geom_hline(yintercept=median(metrics_all_samples$solidity), color="grey30")+
  geom_violin(fill="grey",color="grey")+geom_boxplot(aes(fill=CellType),width=0.1, outlier.shape=NA)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("Solidity")+
  theme(legend.position = "none")+  fillscale_cellType+xlab("")
save_plts(solid, "Solidity_cell_type", w=10, h=5)
