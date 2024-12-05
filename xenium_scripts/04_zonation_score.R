
## Load Libraries
library(here)
library(Seurat)

library(reshape2)
library(dplyr)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

library(scales)
library(viridis)
library(gridExtra)
library(class)
library(ggsignif)

library(tiff)
library(sp)


source("xenium_scripts/00_pretty_plots.R")
source("xenium_scripts/00_long_functions.R")

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


zonation_scores <- lapply(1:7, function(x){
  print("###################################################################################################")
  print(samples[x])
  print("###################################################################################################")
  
  counts<-read.csv(count_files[x])
  
  counts$X<-NULL
  rownames(counts) <- counts$cell_id
  counts$cell_id <- NULL
  
  # Transpose the data and convert to sparse matrix.
  mat <- as(t(as.matrix(counts)), "sparseMatrix")
  
  seu <- CreateSeuratObject(counts=mat)
  seu$sample<-strsplit(count_files[x],"/")[[1]][6]
  seu
  
  
  plt_umap_xenium_sample<-plt_umap_xenium[which(plt_umap_xenium$sample%in%unique(seu$sample)),]
  plt_umap_xenium_sample$cell<-sapply(1:nrow(plt_umap_xenium_sample), function(y){strsplit(plt_umap_xenium_sample$cell[y],"_")[[1]][1]})
  plt_umap_xenium_sample<-plt_umap_xenium_sample[match(colnames(seu), plt_umap_xenium_sample$cell),]
  identical(colnames(seu), plt_umap_xenium_sample$cell)
  rownames(plt_umap_xenium)<-plt_umap_xenium$cell
  
  seu<-AddMetaData(seu, plt_umap_xenium_sample)
  
  seu_hepatocytes<-subset(seu, subset = CellType %in% c("Hepatocyte (Cycling)","Hepatocyte (Periportal)","Hepatocyte (Pericentral)"))
  
  seu_hepatocytes <- NormalizeData(seu_hepatocytes)
  seu_hepatocytes <- FindVariableFeatures(seu_hepatocytes, selection.method = "vst")
  seu_hepatocytes <- ScaleData(seu_hepatocytes) 
  seu_hepatocytes <- RunPCA(seu_hepatocytes, npcs = 30, features = rownames(seu_hepatocytes))
  
  # DimPlot(seu_hepatocytes, group.by = "CellType")+colscale_cellType
  # DimPlot(seu_hepatocytes, group.by = "CellType",dims = c(2, 3))+colscale_cellType
  # DimPlot(seu_hepatocytes, group.by = "CellType",dims = c(3, 4))+colscale_cellType
  # 
  
  ## PC2 is a nice proxy for portal versus central
  embed<-as.data.frame(Embeddings(seu_hepatocytes, reduction = "pca"))
  
  zonation_meta<-seu_hepatocytes@meta.data
  zonation_meta$zonation_PC<-embed[,2]
  zonation_meta$zonation_zscore <- (zonation_meta$zonation_PC - mean(zonation_meta$zonation_PC))/sd(zonation_meta$zonation_PC)
  
  zonation_meta<-zonation_meta[,c("sample",  "cell","zonation_PC","zonation_zscore")]
  
  #hist(zonation_meta$zonation_PC)
  
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
  plt_umap_xenium_sample$cell<-as.character(sapply(1:nrow(plt_umap_xenium_sample), function(x) strsplit(plt_umap_xenium_sample$cell[x],"_")[[1]][1]))
  centroids$cell_id<-as.character(centroids$cell_id)
  if(nrow(centroids)>100000){centroids$cell_id[which(centroids$cell_id=="1e+05")]<-"100000"}
  
  cell_centroid<-merge(plt_umap_xenium_sample, centroids, by.x="cell",by.y="cell_id")
  
  
  zonation_meta_sample<-zonation_meta[which(zonation_meta$sample == smpl),]
  zonation_meta_sample$cell<-sapply(1:nrow(zonation_meta_sample), function(x) strsplit(zonation_meta_sample$cell[x],"_")[[1]][1])
  zonation_meta_sample$sample<-NULL
  
  zonation_centroid<-merge(cell_centroid, zonation_meta_sample, by="cell")
  
  cell_centroid_not_hep<-cell_centroid[which(!(cell_centroid$cell %in% zonation_centroid$cell)),]
  cell_centroid_not_hep$zonation_PC<-NA
  cell_centroid_not_hep$zonation_zscore<-NA
  
  #########################
  ## maxima and minima
  #########################
  print(paste("Maxima and minima from:", nrow(zonation_centroid), "hepatocytes"))
  # Find optimal periportal parameters
  optimal_params <- find_optimal_parameters(zonation_centroid$centroid_x, zonation_centroid$centroid_y, zonation_centroid$zonation_PC, points_per_maxima = 100)
  print(optimal_params)
  maxima <- find_local_maxima(zonation_centroid$centroid_x, zonation_centroid$centroid_y, zonation_centroid$zonation_PC, optimal_params$threshold, optimal_params$k)
  print(nrow(maxima))
  
  # Find optimal pericentral parameters
  optimal_params <- find_optimal_parameters(zonation_centroid$centroid_x, zonation_centroid$centroid_y, -zonation_centroid$zonation_PC, points_per_maxima = 100)
  print(optimal_params)
  minima  <- find_local_maxima(zonation_centroid$centroid_x, zonation_centroid$centroid_y, -zonation_centroid$zonation_PC, optimal_params$threshold, optimal_params$k)
  print(nrow(minima))
  
  
  # ggplot()+
  #   geom_point(aes(centroid_x, -centroid_y,  color=zonation_PC),zonation_centroid,  size=0.5)+
  #   scale_color_viridis(option="mako",name="Non-Hepatocyte Cells\nDistance to\nPeriportal Region", direction=-1)+
  #   geom_point(aes(x, -y),maxima, color="red")
  # 
  # ggplot()+
  #   geom_point(aes(centroid_x, -centroid_y, color=zonation_PC),zonation_scores[[3]][which(!is.na(zonation_scores[[3]]$zonation_PC)),], size=0.5)+
  #   geom_point(aes(centroid_x, -centroid_y),zonation_scores[[3]][which(is.na(zonation_scores[[3]]$zonation_PC)),], color="grey20",size=0.5)+
  #   scale_color_viridis()+
  #   facet_wrap(~CellType)+
  #   geom_point(aes(x, -y),minima, color="red")
  # 
  
  cell_centroid_not_hep$distance_to_periportal <- calculate_distance_to_maxima(cell_centroid_not_hep$centroid_x, cell_centroid_not_hep$centroid_y, maxima$x, maxima$y)
  cell_centroid_not_hep$distance_to_pericentral <- calculate_distance_to_maxima(cell_centroid_not_hep$centroid_x, cell_centroid_not_hep$centroid_y, minima$x, minima$y)
  
  zonation_centroid$distance_to_periportal<-NA
  zonation_centroid$distance_to_pericentral<-NA
  
  rbind(zonation_centroid, cell_centroid_not_hep)
  
})

## zonation score
zonation_scores[[3]]

ggplot(zonation_scores[[3]], aes(CellType,zonation_PC))+geom_boxplot()

###############
### Orient same periportal and pericentral direction based on hep and colangiocyte labels
###############

plot_grid(ggplot(zonation_scores[[1]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[2]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[3]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[4]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[5]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[6]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[7]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip())

plot_grid(ggplot(zonation_scores[[1]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[2]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[3]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[4]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[5]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[6]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[7]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip())


# 
# zonation_scores[[1]]$zonation_zscore<-(-zonation_scores[[1]]$zonation_zscore)
# zonation_scores[[2]]$zonation_zscore<-(-zonation_scores[[2]]$zonation_zscore)
# 
# plot_grid(ggplot(zonation_scores[[1]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip(),
#           ggplot(zonation_scores[[2]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip(),
#           ggplot(zonation_scores[[3]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip(),
#           ggplot(zonation_scores[[4]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip(),
#           ggplot(zonation_scores[[5]], aes(CellType,zonation_zscore))+geom_boxplot()+coord_flip())

## flip max and min
plot_grid(ggplot(zonation_scores[[1]], aes(CellType,distance_to_periportal))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[2]], aes(CellType,distance_to_periportal))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[3]], aes(CellType,distance_to_periportal))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[4]], aes(CellType,distance_to_periportal))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[6]], aes(CellType,distance_to_periportal))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[7]], aes(CellType,distance_to_periportal))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[5]], aes(CellType,distance_to_periportal))+geom_boxplot()+coord_flip())
plot_grid(ggplot(zonation_scores[[1]], aes(CellType,distance_to_pericentral))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[2]], aes(CellType,distance_to_pericentral))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[3]], aes(CellType,distance_to_pericentral))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[4]], aes(CellType,distance_to_pericentral))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[5]], aes(CellType,distance_to_pericentral))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[6]], aes(CellType,distance_to_pericentral))+geom_boxplot()+coord_flip(),
          ggplot(zonation_scores[[7]], aes(CellType,distance_to_pericentral))+geom_boxplot()+coord_flip())


zonation_scores[[1]]$distance_to_central2<-zonation_scores[[1]]$distance_to_pericentral
zonation_scores[[1]]$distance_to_pericentral<-zonation_scores[[1]]$distance_to_periportal
zonation_scores[[1]]$distance_to_periportal<-zonation_scores[[1]]$distance_to_central2
zonation_scores[[1]]$distance_to_central2<-NULL

zonation_scores[[2]]$distance_to_central2<-zonation_scores[[2]]$distance_to_pericentral
zonation_scores[[2]]$distance_to_pericentral<-zonation_scores[[2]]$distance_to_periportal
zonation_scores[[2]]$distance_to_periportal<-zonation_scores[[2]]$distance_to_central2
zonation_scores[[2]]$distance_to_central2<-NULL


save(zonation_scores, file=here("data","zonation_allcells.RData"))


###############################################

load(here("data","zonation_allcells.RData"))

zonation_PC<-ggplot( )+
  geom_point(aes(centroid_x, -centroid_y, color=zonation_zscore),zonation_scores[[1]][which(!is.na(zonation_scores[[1]]$zonation_PC)),], size=0.5)+
  geom_point(aes(centroid_x, -centroid_y),zonation_scores[[1]][which(is.na(zonation_scores[[1]]$zonation_PC)),],color="grey20", size=0.5)+
  theme_presentation()+scale_color_viridis()
zonation_PC
save_plts_black(zonation_PC, "zonation_PC_C94_2", w=15, h=10)

zonation_PC<-ggplot( )+
  geom_point(aes(centroid_x, -centroid_y, color=zonation_zscore),zonation_scores[[2]][which(!is.na(zonation_scores[[2]]$zonation_PC)),], size=0.5)+
  geom_point(aes(centroid_x, -centroid_y),zonation_scores[[2]][which(is.na(zonation_scores[[2]]$zonation_PC)),],color="grey20", size=0.5)+
  theme_presentation()+scale_color_viridis()
zonation_PC
save_plts_black(zonation_PC, "zonation_PC_C94_3", w=15, h=10)

zonation_PC<-ggplot( )+
  geom_point(aes(centroid_x, -centroid_y, color=zonation_zscore),zonation_scores[[3]][which(!is.na(zonation_scores[[3]]$zonation_PC)),], size=0.5)+
  geom_point(aes(centroid_x, -centroid_y),zonation_scores[[3]][which(is.na(zonation_scores[[3]]$zonation_PC)),],color="grey20", size=0.5)+
  theme_presentation()+scale_color_viridis()
zonation_PC
save_plts_black(zonation_PC, "zonation_PC_C94_4", w=15, h=10)

zonation_PC<-ggplot( )+
  geom_point(aes(centroid_x, -centroid_y, color=zonation_zscore),zonation_scores[[4]][which(!is.na(zonation_scores[[4]]$zonation_PC)),], size=0.5)+
  geom_point(aes(centroid_x, -centroid_y),zonation_scores[[4]][which(is.na(zonation_scores[[4]]$zonation_PC)),],color="grey20", size=0.5)+
  theme_presentation()+scale_color_viridis()
zonation_PC
save_plts_black(zonation_PC, "zonation_PC_C105", w=15, h=10)

zonation_PC<-ggplot( )+
  geom_point(aes(centroid_x, -centroid_y, color=zonation_zscore),zonation_scores[[5]][which(!is.na(zonation_scores[[5]]$zonation_PC)),], size=0.5)+
  geom_point(aes(centroid_x, -centroid_y),zonation_scores[[5]][which(is.na(zonation_scores[[5]]$zonation_PC)),],color="grey20", size=0.5)+
  theme_presentation()+scale_color_viridis()
zonation_PC
save_plts_black(zonation_PC, "zonation_PC_C85", w=15, h=10)

zonation_PC<-ggplot( )+
  geom_point(aes(centroid_x, -centroid_y, color=zonation_zscore),zonation_scores[[6]][which(!is.na(zonation_scores[[6]]$zonation_PC)),], size=0.5)+
  geom_point(aes(centroid_x, -centroid_y),zonation_scores[[6]][which(is.na(zonation_scores[[6]]$zonation_PC)),],color="grey20", size=0.5)+
  theme_presentation()+scale_color_viridis()
zonation_PC
save_plts_black(zonation_PC, "zonation_PC_C95", w=15, h=10)

zonation_PC<-ggplot( )+
  geom_point(aes(centroid_x, -centroid_y, color=zonation_zscore),zonation_scores[[7]][which(!is.na(zonation_scores[[7]]$zonation_PC)),], size=0.5)+
  geom_point(aes(centroid_x, -centroid_y),zonation_scores[[7]][which(is.na(zonation_scores[[7]]$zonation_PC)),],color="grey20", size=0.5)+
  theme_presentation()+scale_color_viridis()
zonation_PC
save_plts_black(zonation_PC, "zonation_PC_C101", w=15, h=10)




## non hep plots
periportal<-ggplot()+
  geom_point(aes(centroid_x, -centroid_y),zonation_scores[[3]][which(!is.na(zonation_scores[[3]]$zonation_PC)),], color="grey20", size=0.5)+
  geom_point(aes(centroid_x, -centroid_y, color=distance_to_periportal),zonation_scores[[3]][which(is.na(zonation_scores[[3]]$zonation_PC)),],size=0.5)+
  theme_presentation()+scale_color_viridis(option="mako",name="Non-Hepatocyte Cells\nDistance to\nPeriportal Region", direction=-1)
save_plts_black(periportal, "Periportal_PC_C94_4", w=15, h=10)

pericentral<-ggplot()+
  geom_point(aes(centroid_x, -centroid_y),zonation_scores[[3]][which(!is.na(zonation_scores[[3]]$zonation_PC)),], color="grey20", size=0.5)+
  geom_point(aes(centroid_x, -centroid_y, color=distance_to_pericentral),zonation_scores[[3]][which(is.na(zonation_scores[[3]]$zonation_PC)),],size=0.5)+
  theme_presentation()+scale_color_viridis(option="mako",name="Non-Hepatocyte Cells\nDistance to\nPericentral Region", direction=-1)
save_plts_black(pericentral, "Pericentral_PC_C94_4", w=15, h=10)

ggplot(zonation_scores[[3]], aes(distance_to_periportal, distance_to_pericentral))+geom_point()

## non hep plots
periportal<-ggplot()+
  geom_point(aes(centroid_x, -centroid_y),zonation_scores[[4]][which(!is.na(zonation_scores[[4]]$zonation_PC)),], color="grey20", size=0.5)+
  geom_point(aes(centroid_x, -centroid_y, color=distance_to_periportal),zonation_scores[[4]][which(is.na(zonation_scores[[4]]$zonation_PC)),],size=0.5)+
  theme_presentation()+scale_color_viridis(option="mako",name="Non-Hepatocyte Cells\nDistance to\nPeriportal Region", direction=-1)
save_plts_black(periportal, "Periportal_PC_C105", w=15, h=10)

pericentral<-ggplot()+
  geom_point(aes(centroid_x, -centroid_y),zonation_scores[[4]][which(!is.na(zonation_scores[[4]]$zonation_PC)),], color="grey20", size=0.5)+
  geom_point(aes(centroid_x, -centroid_y, color=distance_to_pericentral),zonation_scores[[4]][which(is.na(zonation_scores[[4]]$zonation_PC)),],size=0.5)+
  theme_presentation()+scale_color_viridis(option="mako",name="Non-Hepatocyte Cells\nDistance to\nPericentral Region", direction=-1)
save_plts_black(pericentral, "Pericentral_PC_C105", w=15, h=10)

ggplot(zonation_scores[[4]], aes(distance_to_periportal, distance_to_pericentral))+geom_point()


## non hep plots
periportal<-ggplot()+
  geom_point(aes(centroid_x, -centroid_y),zonation_scores[[5]][which(!is.na(zonation_scores[[5]]$zonation_PC)),], color="grey20", size=0.5)+
  geom_point(aes(centroid_x, -centroid_y, color=distance_to_periportal),zonation_scores[[5]][which(is.na(zonation_scores[[5]]$zonation_PC)),],size=0.5)+
  theme_presentation()+scale_color_viridis(option="mako",name="Non-Hepatocyte Cells\nDistance to\nPeriportal Region", direction=-1)
save_plts_black(periportal, "Periportal_PC_C85", w=15, h=10)

pericentral<-ggplot()+
  geom_point(aes(centroid_x, -centroid_y),zonation_scores[[5]][which(!is.na(zonation_scores[[5]]$zonation_PC)),], color="grey20", size=0.5)+
  geom_point(aes(centroid_x, -centroid_y, color=distance_to_pericentral),zonation_scores[[5]][which(is.na(zonation_scores[[5]]$zonation_PC)),],size=0.5)+
  theme_presentation()+scale_color_viridis(option="mako",name="Non-Hepatocyte Cells\nDistance to\nPericentral Region", direction=-1)
save_plts_black(pericentral, "Pericentral_PC_C85", w=15, h=10)

ggplot(zonation_scores[[5]], aes(distance_to_periportal, distance_to_pericentral))+geom_point()

##############
## Cell SPA metrics
##############
smpl<-"C94_4"
load(paste(here("data/"),smpl,"_centroid_cellSPA_metrics.RData",sep=""))

metrics_all_samples$CellType<-NULL
metrics_all_samples$centroid_x<-NULL
metrics_all_samples$centroid_y<-NULL
metrics_all_samples$sample<-NULL

zonation_metrics<-merge(zonation_scores[[3]], metrics_all_samples, by="cell")

ggplot(zonation_metrics, aes(distance_to_periportal, cell_area))+
  geom_point()+stat_smooth(method="lm", se=F)+
  facet_wrap(~CellType)+theme_bw()
ggplot(zonation_metrics, aes(distance_to_pericentral, cell_area))+
  geom_point()+stat_smooth(method="lm", se=F)+
  facet_wrap(~CellType)+theme_bw()


ggplot(zonation_metrics, aes(distance_to_periportal, elongation))+
  geom_point()+stat_smooth(method="lm", se=F)+
  facet_wrap(~CellType)+theme_bw()
ggplot(zonation_metrics, aes(distance_to_pericentral, elongation))+
  geom_point()+stat_smooth(method="lm", se=F)+
  facet_wrap(~CellType)+theme_bw()

ggplot(zonation_metrics, aes(distance_to_periportal, density))+
  geom_point()+stat_smooth(method="lm", se=F)+
  facet_wrap(~CellType)+theme_bw()
ggplot(zonation_metrics, aes(distance_to_pericentral, density))+
  geom_point()+stat_smooth(method="lm", se=F)+
  facet_wrap(~CellType)+theme_bw()





## C85
smpl<-"C85"
load(paste(here("data/"),smpl,"_centroid_cellSPA_metrics.RData",sep=""))

metrics_all_samples$CellType<-NULL
metrics_all_samples$centroid_x<-NULL
metrics_all_samples$centroid_y<-NULL
metrics_all_samples$sample<-NULL

zonation_metrics<-merge(zonation_scores[[1]], metrics_all_samples, by="cell")

ggplot(zonation_metrics, aes(distance_to_periportal, cell_area))+
  geom_point()+stat_smooth(method="lm", se=F)+
  facet_wrap(~CellType)+theme_bw()
ggplot(zonation_metrics, aes(distance_to_pericentral, cell_area))+
  geom_point()+stat_smooth(method="lm", se=F)+
  facet_wrap(~CellType)+theme_bw()


ggplot(zonation_metrics, aes(distance_to_periportal, elongation))+
  geom_point()+stat_smooth(method="lm", se=F)+
  facet_wrap(~CellType)+theme_bw()
ggplot(zonation_metrics, aes(distance_to_pericentral, elongation))+
  geom_point()+stat_smooth(method="lm", se=F)+
  facet_wrap(~CellType)+theme_bw()

ggplot(zonation_metrics, aes(distance_to_periportal, density))+
  geom_point()+stat_smooth(method="lm", se=F)+
  facet_wrap(~CellType)+theme_bw()
ggplot(zonation_metrics, aes(distance_to_pericentral, density))+
  geom_point()+stat_smooth(method="lm", se=F)+
  facet_wrap(~CellType)+theme_bw()





######################
## Maxima mini rep plot
######################
lapply(1:7, function(x){
  print("###################################################################################################")
  print(samples[x])
  print("###################################################################################################")
  
  counts<-read.csv(count_files[x])
  
  counts$X<-NULL
  rownames(counts) <- counts$cell_id
  counts$cell_id <- NULL
  
  # Transpose the data and convert to sparse matrix.
  mat <- as(t(as.matrix(counts)), "sparseMatrix")
  
  seu <- CreateSeuratObject(counts=mat)
  seu$sample<-strsplit(count_files[x],"/")[[1]][6]
  seu
  
  
  plt_umap_xenium_sample<-plt_umap_xenium[which(plt_umap_xenium$sample%in%unique(seu$sample)),]
  plt_umap_xenium_sample$cell<-sapply(1:nrow(plt_umap_xenium_sample), function(y){strsplit(plt_umap_xenium_sample$cell[y],"_")[[1]][1]})
  plt_umap_xenium_sample<-plt_umap_xenium_sample[match(colnames(seu), plt_umap_xenium_sample$cell),]
  identical(colnames(seu), plt_umap_xenium_sample$cell)
  rownames(plt_umap_xenium)<-plt_umap_xenium$cell
  
  seu<-AddMetaData(seu, plt_umap_xenium_sample)
  
  seu_hepatocytes<-subset(seu, subset = CellType %in% c("Hepatocyte (Cycling)","Hepatocyte (Periportal)","Hepatocyte (Pericentral)"))
  
  seu_hepatocytes <- NormalizeData(seu_hepatocytes)
  seu_hepatocytes <- FindVariableFeatures(seu_hepatocytes, selection.method = "vst")
  seu_hepatocytes <- ScaleData(seu_hepatocytes) 
  seu_hepatocytes <- RunPCA(seu_hepatocytes, npcs = 30, features = rownames(seu_hepatocytes))
  
  ## PC2 is a nice proxy for portal versus central
  embed<-as.data.frame(Embeddings(seu_hepatocytes, reduction = "pca"))
  
  zonation_meta<-seu_hepatocytes@meta.data
  zonation_meta$zonation_PC<-embed[,2]
  zonation_meta$zonation_zscore <- (zonation_meta$zonation_PC - mean(zonation_meta$zonation_PC))/sd(zonation_meta$zonation_PC)
  
  zonation_meta<-zonation_meta[,c("sample",  "cell","zonation_PC","zonation_zscore")]
  
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
  
  zonation_meta_sample<-zonation_meta[which(zonation_meta$sample == smpl),]
  zonation_meta_sample$cell<-sapply(1:nrow(zonation_meta_sample), function(x) strsplit(zonation_meta_sample$cell[x],"_")[[1]][1])
  zonation_meta_sample$sample<-NULL
  
  zonation_centroid<-merge(cell_centroid, zonation_meta_sample, by="cell")
  
  cell_centroid_not_hep<-cell_centroid[which(!(cell_centroid$cell %in% zonation_centroid$cell)),]
  cell_centroid_not_hep$zonation_PC<-NA
  cell_centroid_not_hep$zonation_zscore<-NA
  
  #########################
  ## maxima and minima
  #########################
  print(paste("Maxima and minima from:", nrow(zonation_centroid), "hepatocytes"))
  # Find optimal periportal parameters
  optimal_params <- find_optimal_parameters(zonation_centroid$centroid_x, zonation_centroid$centroid_y, zonation_centroid$zonation_PC, points_per_maxima = 100)
  print(optimal_params)
  maxima <- find_local_maxima(zonation_centroid$centroid_x, zonation_centroid$centroid_y, zonation_centroid$zonation_PC, optimal_params$threshold, optimal_params$k)
  print(nrow(maxima))
  
  # Find optimal pericentral parameters
  optimal_params <- find_optimal_parameters(zonation_centroid$centroid_x, zonation_centroid$centroid_y, -zonation_centroid$zonation_PC, points_per_maxima = 100)
  print(optimal_params)
  minima  <- find_local_maxima(zonation_centroid$centroid_x, zonation_centroid$centroid_y, -zonation_centroid$zonation_PC, optimal_params$threshold, optimal_params$k)
  print(nrow(minima))

  max_min<-ggplot()+
    geom_point(aes(centroid_x, -centroid_y, color=zonation_PC),zonation_centroid[which(!is.na(zonation_centroid$zonation_PC)),], size=0.5)+
    geom_point(aes(centroid_x, -centroid_y),zonation_centroid[which(is.na(zonation_centroid$zonation_PC)),], color="grey20",size=0.5)+
    scale_color_viridis()+theme_presentation()+
    geom_point(aes(x, -y),minima, color="red")+
    geom_point(aes(x, -y),maxima, color="white")
  save_plts_black(max_min, paste(samples[x], "_zonation_max_min", sep=""), w=15, h=10)
})