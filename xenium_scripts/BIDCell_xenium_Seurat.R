## Load Libraries
library(here)
library(Seurat)

library(SCINA)
library(reshape2)
library(dplyr)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

library(sp)
library(rgeos)
library(scales)
library(viridis)
library(colorspace)


count_files<-c("liver_BIDCell_output/C94_2/cell_gene_matrices/2024_03_21_16_22_10/expr_mat.csv",
               "liver_BIDCell_output/C94_3/cell_gene_matrices/2024_03_21_16_38_48/expr_mat.csv",
               "liver_BIDCell_output/C94_4/cell_gene_matrices/2024_03_21_17_03_05/expr_mat.csv",
               "liver_BIDCell_output/C105/cell_gene_matrices/2024_05_17_13_32_03/expr_mat.csv",
               "liver_BIDCell_output/C85/cell_gene_matrices/2024_05_17_13_29_40/expr_mat.csv",
               "liver_BIDCell_output/C95/cell_gene_matrices/2024_07_04_12_35_29/expr_mat.csv",
               "liver_BIDCell_output/C101/cell_gene_matrices/2024_06_28_09_23_40/expr_mat.csv")



d10x.list <- sapply(count_files, function(file_path){
  counts<-read.csv(file_path)
  
  counts$X<-NULL
  rownames(counts) <- counts$cell_id
  counts$cell_id <- NULL
  
  # Transpose the data and convert to sparse matrix.
  mat <- as(t(as.matrix(counts)), "sparseMatrix")
  
  seu <- CreateSeuratObject(counts=mat)
  seu$sample<-strsplit(file_path,"/")[[1]][6]
  seu
})

d10x.list

xenium.obj <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "xenium_liver")
rm(d10x.list)
gc()


load(file=here("data/cell_type_labels_BIDCell.RData"))

plt_umap_xenium<-plt_umap_xenium[match(colnames(xenium.obj), plt_umap_xenium$cell),]
identical(colnames(xenium.obj), plt_umap_xenium$cell)
rownames(plt_umap_xenium)<-plt_umap_xenium$cell

xenium.obj<-AddMetaData(xenium.obj, plt_umap_xenium)