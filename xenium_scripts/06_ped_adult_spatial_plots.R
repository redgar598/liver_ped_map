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
library(tiff)

library(sp)
library(scales)
library(viridis)

source("scripts/00_pretty_plots.R")
source("scripts/00_long_functions.R")

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

lapply(1:7, function(x){
  data_dir<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x],"/cell_gene_matrices/", models[x], "/expr_mat.csv", sep="")
  tiff_path<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x], "/model_outputs/",models[x],"/test_output/epoch_1_step_4000_connected.tif", sep="")
  smpl<-samples[x]
  print(smpl)
  
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
  
  ## count data
  counts<-read.csv(data_dir)
  counts$X<-NULL
  rownames(counts) <- counts$cell_id
  counts$cell_id <- NULL
  mat <- as(t(as.matrix(counts)), "sparseMatrix")
  
  
  
  
  cell_centroid<-cell_centroid[(order(cell_centroid$CellType)),]
  hep<-ggplot()+
    geom_point(aes(centroid_x,-centroid_y),cell_centroid, color="black", size=1, shape=19)+
    geom_point(aes(centroid_x,-centroid_y),cell_centroid, color="grey95", size=0.1, shape=19)+
    geom_point(aes(centroid_x,-centroid_y, color=CellType),cell_centroid[which(cell_centroid$CellType%in%c("KC Like","Hepatocyte (Periportal)", "Hepatocyte (Pericentral)", "Hepatocyte (Cycling)","Cholangiocytes")),], size=0.5, shape=19)+
    theme_void()+
    guides(colour = guide_legend(override.aes = list(size=5))) +
    colscale_cellType
  hep<-hep+xenium_scale_bar(cell_centroid, 3)
  save_plts(hep,paste(smpl,"_fancy", sep=""),  w=10, h=6)
  
  
  save_plts(get_leg(hep),paste(smpl,"_fancy_legend", sep=""),  w=2, h=3)
  
  load(here("data","zonation_allcells.RData"))
  
  zonation_plt<-ggplot()+
    geom_point(aes(centroid_x,-centroid_y),zonation_scores[[x]], color="black", size=1, shape=19)+
    geom_point(aes(centroid_x,-centroid_y),zonation_scores[[x]], color="grey95", size=0.1, shape=19)+
    geom_point(aes(centroid_x,-centroid_y, color=zonation_zscore),zonation_scores[[x]][which(zonation_scores[[x]]$CellType%in%c("Hepatocyte (Periportal)", "Hepatocyte (Pericentral)", "Hepatocyte (Cycling)")),], size=0.25, shape=19)+
    theme_void()+scale_color_viridis(name="Hepatocyte\nZonation\nScore",limits = range(-7.5, 3.5))
  zonation_plt
  save_plts(zonation_plt+xenium_scale_bar(zonation_scores[[x]], 3),paste(smpl,"_poster_zonation", sep=""),  w=10, h=6)
  
  save_plts(get_leg(hep),paste(smpl,"_fancy_legend", sep=""),  w=2, h=3)
  
  
  plot_grid(hep+theme(legend.position = "none"), get_leg(hep),
            zonation_plt+theme(legend.position = "none"), get_leg(zonation_plt),
            ncol=2, align="h", axis="lr")
  save_plts(plot_grid(hep+xenium_scale_bar(cell_centroid, 3)+theme(legend.position = "none"), get_leg(hep),
                      zonation_plt+xenium_scale_bar(zonation_scores[[x]], 3)+theme(legend.position = "none"), get_leg(zonation_plt),
                      ncol=2, align="h", axis="lr"),
            paste(smpl,"_zonation_grid", sep=""),  w=14, h=12)
})






# 
# count_files<-c("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_2/cell_gene_matrices/2024_03_21_16_22_10/expr_mat.csv",
#                "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_3/cell_gene_matrices/2024_03_21_16_38_48/expr_mat.csv",
#                "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_4/cell_gene_matrices/2024_03_21_17_03_05/expr_mat.csv",
#                "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C105/cell_gene_matrices/2024_05_17_13_32_03/expr_mat.csv",
#                "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C85/cell_gene_matrices/2024_05_17_13_29_40/expr_mat.csv")
# 
# 
# samples<-c("C94_2", "C94_3","C94_4","C105","C85")
# models<-c("2024_03_21_16_22_10","2024_03_21_16_38_48","2024_03_21_17_03_05","2024_05_17_13_32_03","2024_05_17_13_29_40")
# 
# 
# x=4
# 
# data_dir<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x],"/cell_gene_matrices/", models[x], "/expr_mat.csv", sep="")
# tiff_path<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x], "/model_outputs/",models[x],"/test_output/epoch_1_step_4000_connected.tif", sep="")
# smpl<-samples[x]
# print(smpl)
# 
# ##### centroids
# tiff_res <- readTIFF(tiff_path, as.is = TRUE)
# tiff_res <- reshape2::melt(tiff_res)
# tiff_res <- tiff_res[tiff_res$value != 0, ]
# 
# colnames(tiff_res) <- c("coord_y", "coord_x", "cell_id")
# tiff_res <- tiff_res[, c("coord_x", "coord_y", "cell_id")]
# tiff_res$cell_id <- as.numeric(tiff_res$cell_id)
# tiff_res$coord_x <- as.numeric(tiff_res$coord_x)
# tiff_res$coord_y <- as.numeric(tiff_res$coord_y)
# tiff_res <- tiff_res[order(tiff_res$cell_id), ]
# rownames(tiff_res) <- NULL
# 
# centroids<-as.data.frame(
#   tiff_res %>%
#     group_by(cell_id) %>%
#     summarise(centroid_x = sum(coord_x) / length(coord_x), centroid_y = sum(coord_y) / length(coord_y)))
# 
# plt_umap_xenium_sample<-plt_umap_xenium[which(plt_umap_xenium$sample == smpl),]
# plt_umap_xenium_sample$cell<-sapply(1:nrow(plt_umap_xenium_sample), function(x) strsplit(plt_umap_xenium_sample$cell[x],"_")[[1]][1])
# 
# cell_centroid<-merge(plt_umap_xenium_sample, centroids, by.x="cell",by.y="cell_id")
# 
# ## count data
# counts<-read.csv(data_dir)
# counts$X<-NULL
# rownames(counts) <- counts$cell_id
# counts$cell_id <- NULL
# mat <- as(t(as.matrix(counts)), "sparseMatrix")
# 
# 
# 
# 
# cell_centroid<-cell_centroid[(order(cell_centroid$CellType)),]
# hep<-ggplot()+
#   geom_point(aes(centroid_x,-centroid_y),cell_centroid, color="black", size=1, shape=19)+
#   geom_point(aes(centroid_x,-centroid_y),cell_centroid, color="grey95", size=0.1, shape=19)+
#   geom_point(aes(centroid_x,-centroid_y, color=CellType),cell_centroid[which(cell_centroid$CellType%in%c("KC Like","Hepatocyte (Periportal)", "Hepatocyte (Pericentral)", "Hepatocyte (Cycling)","Cholangiocytes")),], size=0.5, shape=19)+
#   theme_void()+
#   guides(colour = guide_legend(override.aes = list(size=5))) +
#   colscale_cellType
# hep
# save_plts(hep,"C105_fancy",  w=10, h=6)
# 
# 
# save_plts(get_leg(hep),"C105_fancy_legend",  w=2, h=3)
# 
# 
# 
# 
# load(here("data","zonation_allcells.RData"))
# 
# zonation_plt<-ggplot()+
#   geom_point(aes(centroid_x,-centroid_y),zonation_scores[[x]], color="black", size=1, shape=19)+
#   geom_point(aes(centroid_x,-centroid_y),zonation_scores[[x]], color="grey95", size=0.1, shape=19)+
#   geom_point(aes(centroid_x,-centroid_y, color=zonation_zscore),zonation_scores[[x]][which(zonation_scores[[x]]$CellType%in%c("Hepatocyte (Periportal)", "Hepatocyte (Pericentral)", "Hepatocyte (Cycling)")),], size=0.25, shape=19)+
#   theme_void()+scale_color_viridis(name="Hepatocyte\nZonation\nScore",limits = range(-7.5, 3.5))
# zonation_plt
# save_plts(zonation_plt,"C105_poster_zonation",  w=10, h=6)
# 
# save_plts(get_leg(hep),"C105_fancy_legend",  w=2, h=3)
# 
# 
# plot_grid(hep+theme(legend.position = "none"), get_leg(hep),
#           zonation_plt+theme(legend.position = "none"), get_leg(zonation_plt),
#           ncol=2, align="h", axis="lr")
# save_plts(plot_grid(hep+theme(legend.position = "none"), get_leg(hep),
#                     zonation_plt+theme(legend.position = "none"), get_leg(zonation_plt),
#                     ncol=2, align="h", axis="lr"),
#           "C105_zonation_grid",  w=14, h=12)
# 
# #######################################
# # other sample
# #######################################
# x=3
# 
# data_dir<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x],"/cell_gene_matrices/", models[x], "/expr_mat.csv", sep="")
# tiff_path<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x], "/model_outputs/",models[x],"/test_output/epoch_1_step_4000_connected.tif", sep="")
# smpl<-samples[x]
# print(smpl)
# 
# ##### centroids
# tiff_res <- readTIFF(tiff_path, as.is = TRUE)
# tiff_res <- reshape2::melt(tiff_res)
# tiff_res <- tiff_res[tiff_res$value != 0, ]
# 
# colnames(tiff_res) <- c("coord_y", "coord_x", "cell_id")
# tiff_res <- tiff_res[, c("coord_x", "coord_y", "cell_id")]
# tiff_res$cell_id <- as.numeric(tiff_res$cell_id)
# tiff_res$coord_x <- as.numeric(tiff_res$coord_x)
# tiff_res$coord_y <- as.numeric(tiff_res$coord_y)
# tiff_res <- tiff_res[order(tiff_res$cell_id), ]
# rownames(tiff_res) <- NULL
# 
# centroids<-as.data.frame(
#   tiff_res %>%
#     group_by(cell_id) %>%
#     summarise(centroid_x = sum(coord_x) / length(coord_x), centroid_y = sum(coord_y) / length(coord_y)))
# 
# plt_umap_xenium_sample<-plt_umap_xenium[which(plt_umap_xenium$sample == smpl),]
# plt_umap_xenium_sample$cell<-sapply(1:nrow(plt_umap_xenium_sample), function(x) strsplit(plt_umap_xenium_sample$cell[x],"_")[[1]][1])
# 
# cell_centroid<-merge(plt_umap_xenium_sample, centroids, by.x="cell",by.y="cell_id")
# 
# ## count data
# counts<-read.csv(data_dir)
# counts$X<-NULL
# rownames(counts) <- counts$cell_id
# counts$cell_id <- NULL
# mat <- as(t(as.matrix(counts)), "sparseMatrix")
# 
# 
# 
# 
# cell_centroid<-cell_centroid[(order(cell_centroid$CellType)),]
# hep<-ggplot()+
#   geom_point(aes(centroid_x,-centroid_y),cell_centroid, color="black", size=1, shape=19)+
#   geom_point(aes(centroid_x,-centroid_y),cell_centroid, color="grey95", size=0.1, shape=19)+
#   geom_point(aes(centroid_x,-centroid_y, color=CellType),cell_centroid[which(cell_centroid$CellType%in%c("KC Like","Hepatocyte (Periportal)", "Hepatocyte (Pericentral)", "Hepatocyte (Cycling)","Cholangiocytes")),], size=0.5, shape=19)+
#   theme_void()+
#   guides(colour = guide_legend(override.aes = list(size=5))) +
#   colscale_cellType
# hep
# save_plts(hep,"C94_4_fancy",  w=10, h=6)
# 
# 
# save_plts(get_leg(hep),"C94_4_fancy_legend",  w=2, h=3)
# 
# 
# 
# 
# load(here("data","zonation_allcells.RData"))
# 
# zonation_plt<-ggplot()+
#   geom_point(aes(centroid_x,-centroid_y),zonation_scores[[x]], color="black", size=1, shape=19)+
#   geom_point(aes(centroid_x,-centroid_y),zonation_scores[[x]], color="grey95", size=0.1, shape=19)+
#   geom_point(aes(centroid_x,-centroid_y, color=zonation_zscore),zonation_scores[[x]][which(zonation_scores[[x]]$CellType%in%c("Hepatocyte (Periportal)", "Hepatocyte (Pericentral)", "Hepatocyte (Cycling)")),], size=0.25, shape=19)+
#   theme_void()+scale_color_viridis(name="Hepatocyte\nZonation\nScore",limits = range(-7.5, 3.5))
# zonation_plt
# save_plts(zonation_plt,"C94_4_poster_zonation",  w=10, h=6)
# 
# save_plts(get_leg(hep),"C94_4_fancy_legend",  w=2, h=3)
# 
# 
# plot_grid(hep+theme(legend.position = "none"), get_leg(hep),
#           zonation_plt+theme(legend.position = "none"), get_leg(zonation_plt),
#           ncol=2, align="h", axis="lr")
# save_plts(plot_grid(hep+theme(legend.position = "none"), get_leg(hep),
#                     zonation_plt+theme(legend.position = "none"), get_leg(zonation_plt),
#                     ncol=2, align="h", axis="lr"),
#           "C94_4_zonation_grid",  w=14, h=12)
