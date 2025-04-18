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


source("scRNA_seq_scripts/00_pretty_plots.R")
source("scRNA_seq_scripts/00_fanciest_UMAP.R")
source("scRNA_seq_scripts/00_plot_gene_exp.R")



load("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/IFALD_adult_ped_integrated_healthy_only.rds")



  umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined_healthy, reduction = "umap"))#
  umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
  meta_myeloid<-d10x.combined_healthy@meta.data
  meta_myeloid$cell<-rownames(meta_myeloid)
  plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")
  
  len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
  len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
  arr <- list(x = min(plt_myeloid$UMAP_1), y = min(plt_myeloid$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)
  
    fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
      geom_point(size = 0.06, colour= "black", stroke = 1)+
      geom_point(aes(color=CellType_refined),size=0.05)+xlab("")+ylab("")+
      colscale_cellType+
      theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                         axis.title.x = element_text(size=10,hjust = 0.05),
                         axis.title.y = element_text(size=10,hjust = 0.05,angle = 90),
                         legend.position = "none")
  


  
save_plts(fanciest_UMAP, "Healthy_only_refined_cellType_umpa_fancy_lowres", w=2,h=2)



