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

source("scripts/00_pretty_plots.R")
#source("scripts/00_long_functions_BIDCell.R")

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



lapply(1:length(samples), function(x){
  
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
  
  
  ###############
  ## Cell type plots
  ###############
  cell_type_plot<-ggplot(cell_centroid, aes(centroid_x,-centroid_y))+geom_point(aes(color=CellType),size=0.05, shape=19)+
    theme_presentation()+ggtitle(smpl)+
    guides(colour = guide_legend(override.aes = list(size=5))) +
    colscale_cellType
  cell_type_plot
  
  save_plts_black(plot_grid(cell_type_plot+ guides(color = FALSE, size = FALSE), 
                            get_leg(cell_type_plot), align="vh", rel_widths = c(1,0.75)), 
                  paste(smpl,"_CellType_BIDCell", sep=""),  w=20, h=10)
  
  
  
  blank_plot<-ggplot(cell_centroid, aes(centroid_x,-centroid_y))+
    geom_point(size=0.05, shape=19, alpha=0.25, color="grey")+
    theme_presentation()+ggtitle(smpl)
  
  save_plts_black(blank_plot,paste(smpl,"_centroids_BIDCell", sep=""),  w=8, h=6)
  
  

  kc_plot<-ggplot()+
    geom_point(aes(centroid_x,-centroid_y),cell_centroid, color="grey20", size=0.05, shape=19)+
    geom_point(aes(centroid_x,-centroid_y, color=CellType),cell_centroid[which(cell_centroid$CellType%in%c("Macrophage MHCII High","KC Like","Mono-Mac")),], size=0.5, shape=19)+
    theme_presentation()+ggtitle(smpl)+
    guides(colour = guide_legend(override.aes = list(size=5))) +
    colscale_cellType
  kc_plot
  save_plts_black(plot_grid(kc_plot+ guides(color = FALSE, size = FALSE), 
                            get_leg(kc_plot), align="vh", rel_widths = c(1,0.75)), 
                  paste(smpl,"_KC_BIDCell", sep=""),  w=20, h=10)
  
  
  hep<-ggplot()+
    geom_point(aes(centroid_x,-centroid_y),cell_centroid, color="grey20", size=0.05, shape=19)+
    geom_point(aes(centroid_x,-centroid_y, color=CellType),cell_centroid[which(cell_centroid$CellType%in%c("Hepatocyte (Periportal)", "Hepatocyte (Pericentral)", "Hepatocyte (Cycling)","Cholangiocytes","Cholangiocyte (Biliary)")),], size=0.5, shape=19)+
    theme_presentation()+ggtitle(smpl)+
    guides(colour = guide_legend(override.aes = list(size=5))) +
    colscale_cellType
  hep
  save_plts_black(plot_grid(hep+ guides(color = FALSE, size = FALSE), 
                            get_leg(hep), align="vh", rel_widths = c(1,0.75)), 
                  paste(smpl,"_hepatocytes_BIDCell", sep=""),  w=20, h=10)

  
  mesenchyme<-ggplot()+
    geom_point(aes(centroid_x,-centroid_y),cell_centroid, color="grey20", size=0.05, shape=19)+
    geom_point(aes(centroid_x,-centroid_y, color=CellType),cell_centroid[which(cell_centroid$CellType%in%c( "HSC (Periportal)","HSC (Quiescent)","HSC","HSC (Activated)","Cholangiocytes","Cholangiocyte (Biliary)")),], size=0.5, shape=19)+
    theme_presentation()+ggtitle(smpl)+
    guides(colour = guide_legend(override.aes = list(size=5))) +
    colscale_cellType
  save_plts_black(plot_grid(mesenchyme+ guides(color = FALSE, size = FALSE), 
                            get_leg(mesenchyme), align="vh", rel_widths = c(1,0.75)), 
                  paste(smpl,"_HSC_BIDCell", sep=""),  w=20, h=10)
  
  mesenchyme<-ggplot()+
    geom_point(aes(centroid_x,-centroid_y),cell_centroid, color="grey20", size=0.05, shape=19)+
    geom_point(aes(centroid_x,-centroid_y, color=CellType),cell_centroid[which(cell_centroid$CellType%in%c("Cholangiocytes","Cholangiocyte (Biliary)", "VEC","LSEC I","LSEC II","LSEC (Periportal)")),], size=0.5, shape=19)+
    theme_presentation()+ggtitle(smpl)+
    guides(colour = guide_legend(override.aes = list(size=5))) +
    colscale_cellType
  save_plts_black(plot_grid(mesenchyme+ guides(color = FALSE, size = FALSE), 
                            get_leg(mesenchyme), align="vh", rel_widths = c(1,0.75)), 
                  paste(smpl,"_LSEC_BIDCell", sep=""),  w=20, h=10)
  })
                                   



######################
## Summary Stats
######################
load(file="/media/redgar/Seagate Portable Drive/liver_BIDCell_output/BIDCell_liver_seurat.RData")

table(xenium.obj$sample)

xenium.obj@meta.data


lapply(samples, function(smple){
  df<-xenium.obj@meta.data[which(xenium.obj@meta.data$sample == smple),]
  
  print(smple)
  print(paste("Median transcripts per cell", median(df$nCount_RNA)))
  print(paste("Mean transcripts per cell", mean(df$nCount_RNA)))
  
  print(paste("Median genes per cell", median(df$nFeature_RNA)))
  print(paste("Mean genes per cell", mean(df$nFeature_RNA)))
})


######################
## 10X QC metrics 
######################
samples<-c("C94_2","C94_3","C94_4","C85","C105","C95","C101")


d10x.list <- sapply(samples, function(smple){
  #load(file=paste(here("data/"),smple, "_object_raw.RData",sep=""))
  print(smple)
  load(file=paste("/media/redgar/Seagate Portable Drive/xenium_liver/",smple,"_object_raw.RData",sep=""))
  xenium.obj$sample<-smple
  xenium.obj
})



ImageFeaturePlot(d10x.list[[1]], features = c("nCount_BlankCodeword","nFeature_BlankCodeword"), max.cutoff='q50', size = 0.75, cols = c("grey20", "red"))+
  coord_flip()+theme_presentation()+
  theme(legend.title = element_text(size=30),
        legend.text = element_text(size=15),
        plot.title = element_text(size=40))

ImageFeaturePlot(d10x.list[[1]], features = c("nCount_ControlCodeword","nFeature_ControlCodeword"), max.cutoff='q90', size = 0.75, cols = c("grey20", "red"))+
  coord_flip()+theme_presentation()+
  theme(legend.title = element_text(size=30),
        legend.text = element_text(size=15),
        plot.title = element_text(size=40))

ImageFeaturePlot(d10x.list[[1]], features = c("nCount_ControlProbe","nFeature_ControlProbe"), max.cutoff='q90', size = 0.75, cols = c("grey20", "red"))+
  coord_flip()+theme_presentation()+
  theme(legend.title = element_text(size=30),
        legend.text = element_text(size=15),
        plot.title = element_text(size=40))


n_feature<-ImageFeaturePlot(d10x.list[[1]], features = c("nFeature_Xenium"), max.cutoff='q90', size = 0.5, cols = c("grey20", "#e41a1c"))+
  coord_flip()+ggtitle("C94_2")+
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        plot.title = element_text(size=15))


n_counts<-ImageFeaturePlot(d10x.list[[1]], features = c("nCount_Xenium"), max.cutoff='q90', size = 0.5, cols = c("grey10", "#41b6c4"))+
  coord_flip()+theme_presentation()+ggtitle("C94_2")+
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        plot.title = element_text(size=15))

save_plts_black(n_counts, "C94_2_ncounts", w=8, h=6)
save_plts_black(n_feature, "C94_2_features", w=8, h=6)


n_feature<-ImageFeaturePlot(d10x.list[[2]], features = c("nFeature_Xenium"), max.cutoff='q90', size = 0.5, cols = c("grey20", "#e41a1c"))+
  coord_flip()+ggtitle("C94_3")+
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        plot.title = element_text(size=15))


n_counts<-ImageFeaturePlot(d10x.list[[2]], features = c("nCount_Xenium"), max.cutoff='q90', size = 0.5, cols = c("grey10", "#41b6c4"))+
  coord_flip()+theme_presentation()+ggtitle("C94_3")+
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        plot.title = element_text(size=15))

save_plts_black(n_counts, "C94_3_ncounts", w=8, h=6)
save_plts_black(n_feature, "C94_3_features", w=8, h=6)



n_feature<-ImageFeaturePlot(d10x.list[[3]], features = c("nFeature_Xenium"), max.cutoff='q90', size = 0.5, cols = c("grey20", "#e41a1c"))+
  coord_flip()+ggtitle("C94_4")+
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        plot.title = element_text(size=15))


n_counts<-ImageFeaturePlot(d10x.list[[3]], features = c("nCount_Xenium"), max.cutoff='q90', size = 0.5, cols = c("grey10", "#41b6c4"))+
  coord_flip()+theme_presentation()+ggtitle("C94_4")+
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        plot.title = element_text(size=15))

save_plts_black(n_counts, "C94_4_ncounts", w=8, h=6)
save_plts_black(n_feature, "C94_4_features", w=8, h=6)



n_feature<-ImageFeaturePlot(d10x.list[[4]], features = c("nFeature_Xenium"), max.cutoff='q90', size = 0.5, cols = c("grey20", "#e41a1c"))+
  coord_flip()+ggtitle("C85")+
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        plot.title = element_text(size=15))


n_counts<-ImageFeaturePlot(d10x.list[[4]], features = c("nCount_Xenium"), max.cutoff='q90', size = 0.5, cols = c("grey10", "#41b6c4"))+
  coord_flip()+theme_presentation()+ggtitle("C85")+
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        plot.title = element_text(size=15))

save_plts_black(n_counts, "C85_ncounts", w=8, h=6)
save_plts_black(n_feature, "C85_features", w=8, h=6)


n_feature<-ImageFeaturePlot(d10x.list[[5]], features = c("nFeature_Xenium"), max.cutoff='q90', size = 0.5, cols = c("grey20", "#e41a1c"))+
  coord_flip()+ggtitle("C105")+
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        plot.title = element_text(size=15))


n_counts<-ImageFeaturePlot(d10x.list[[5]], features = c("nCount_Xenium"), max.cutoff='q90', size = 0.5, cols = c("grey10", "#41b6c4"))+
  coord_flip()+theme_presentation()+ggtitle("C105")+
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        plot.title = element_text(size=15))

save_plts_black(n_counts, "C105_ncounts", w=8, h=6)
save_plts_black(n_feature, "C105_features", w=8, h=6)




