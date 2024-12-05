####################
##check markers in xenium data
#####################
## Load Libraries
library(here)
library(Seurat)

library(reshape2)
library(dplyr)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

library(sp)
library(scales)
library(gridExtra)

source("scripts/00_pretty_plots.R")


## the panel comes with cell type labels
load(file=here("data/cell_type_labels_BIDCell.RData"))




count_files<-c("/media/redgar/Seagate Portable Drive/stroke_BIDCell_output/output-XETG00082__0011329__tJWH049__20231003__204217/cell_gene_matrices/2024_02_20_15_44_00/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/stroke_BIDCell_output/output-XETG00082__0011329__tJWH054__20231003__204217/cell_gene_matrices/2024_02_25_19_43_04/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/stroke_BIDCell_output/output-XETG00082__0011329__tJWH050__20231003__204217/cell_gene_matrices/2024_02_25_14_51_00/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/stroke_BIDCell_output/output-XETG00082__0011333__tJWH052__20231003__204217/cell_gene_matrices/2024_02_25_22_53_17/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/stroke_BIDCell_output/output-XETG00082__0011333__tJWH051__20231003__204217/cell_gene_matrices/2024_02_25_19_18_28/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/stroke_BIDCell_output/output-XETG00082__0011333__tJWH055__20231003__204217/cell_gene_matrices/2024_02_26_01_58_17/expr_mat.csv")


d10x.list <- sapply(count_files, function(file_path){
  counts<-read.csv(file_path)
  
  counts$X<-NULL
  rownames(counts) <- counts$cell_id
  counts$cell_id <- NULL
  
  # Transpose the data and convert to sparse matrix.
  mat <- as(t(as.matrix(counts)), "sparseMatrix")
  
  seu <- CreateSeuratObject(counts=mat)
  seu$sample<-strsplit(strsplit(file_path,"/")[[1]][6],"__")[[1]][3]
  seu
  
  
  
})

d10x.list

xenium.obj <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "xenium_stroke")
rm(d10x.list)
gc()





###############
## UMAP stroke and contra
###############
samples<-c("output-XETG00082__0011329__tJWH049__20231003__204217",
           "output-XETG00082__0011329__tJWH054__20231003__204217",
           "output-XETG00082__0011329__tJWH050__20231003__204217",
           "output-XETG00082__0011333__tJWH052__20231003__204217",
           "output-XETG00082__0011333__tJWH051__20231003__204217",
           "output-XETG00082__0011333__tJWH055__20231003__204217")
metrics_list <- lapply(1:6, function(x){
  smpl<-strsplit(samples[x],"__")[[1]][3]
  print(smpl)
  load(paste(here("data/"),smpl,"_centroid_cellSPA_metrics.RData",sep=""))
  plt_umap_xenium_sample_metrics_centroids})


metrics_all_samples<-do.call(rbind, metrics_list)
metrics_all_samples$condition<-as.factor(metrics_all_samples$sample)

levels(metrics_all_samples$condition)<-c("Control", "Reprogrammed", "Reprogrammed", "Control", "Reprogrammed", "Reprogrammed")



metrics_all_samples$sample_num<-as.factor(metrics_all_samples$sample)
levels(metrics_all_samples$sample_num)<-c(1,3,5,4,2,6)

metrics_all_samples$cell_id<-sapply(1:nrow(metrics_all_samples), function(x){
  paste(metrics_all_samples$cell[x],"_",metrics_all_samples$sample_num[x], sep="")
})

xenium.obj$cell<-colnames(xenium.obj)

xenium.obj.stroke<-subset(xenium.obj, subset = cell %in% metrics_all_samples$cell_id )

metrics_all_samples<-metrics_all_samples[match(colnames(xenium.obj.stroke), metrics_all_samples$cell_id),]
identical(colnames(xenium.obj.stroke), metrics_all_samples$cell_id)
rownames(metrics_all_samples)<-metrics_all_samples$cell_id

xenium.obj.stroke<-AddMetaData(xenium.obj.stroke, metrics_all_samples)




###############
## all samples
###############

xenium.obj.stroke<-subset(xenium.obj.stroke, subset = region %in% c("contralateral_stroke","stroke_site"))

xenium.obj.stroke <- NormalizeData(xenium.obj.stroke)
xenium.obj.stroke <- FindVariableFeatures(xenium.obj.stroke, selection.method = "vst", nfeatures = 2000)
xenium.obj.stroke <- ScaleData(xenium.obj.stroke)
xenium.obj.stroke <- RunPCA(xenium.obj.stroke, npcs = 30, features = rownames(xenium.obj))
xenium.obj.stroke <- RunUMAP(xenium.obj.stroke, dims = 1:30)



umap_mat_myeloid<-as.data.frame(Embeddings(object = xenium.obj.stroke, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-xenium.obj.stroke@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_umap<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

len_x_bar<-((range(plt_umap$UMAP_1))[2]-(range(plt_umap$UMAP_1))[1])/10
len_y_bar<-((range(plt_umap$UMAP_2))[2]-(range(plt_umap$UMAP_2))[1])/10
arr <- list(x = min(plt_umap$UMAP_1)-1, y = min(plt_umap$UMAP_2)-1, x_len = len_x_bar, y_len = len_y_bar)


forlegned_plot<-ggplot(plt_umap, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=CellType),size=3, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  fillscale_cellType+theme_bw()+
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=12))
nice_legend<-get_leg(forlegned_plot)

adjust_lab_by<-diff(range(plt_umap$UMAP_2))*0.05


fanciest_UMAP<-ggplot(plt_umap, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.5, colour= "black", stroke = 1)+
  geom_point(aes(color=CellType),size=0.45)+      
  xlab("")+ylab("")+
  annotate("text",  x = arr$x+(arr$x_len/2), y = arr$y-adjust_lab_by, label="UMAP 1", size=3)+
  annotate("text",  x = arr$x-adjust_lab_by, y = arr$y+(arr$y_len/2), label="UMAP 2", size=3, angle = 90)+
  colscale_cellType+
  annotate("segment",
           x = arr$x, xend = arr$x + c(arr$x_len, 0),
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=12,hjust = 0.05),
                     axis.title.y = element_text(size=12,hjust = 0.05,angle = 90),
                     legend.position = "none")

## Split
cell_num_all<-as.data.frame(table(xenium.obj.stroke@meta.data$region))
colnames(cell_num_all)<-c("region","CellCount")
fanciest_UMAP <- fanciest_UMAP + facet_wrap(~region, ncol=2)+  theme(strip.text = element_text(size=14))

save_plts(plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2)), "UMAP_allsamples_stroke_contralateral", w=14, h=5)






############
## controls only 
############

xenium.obj.stroke<-subset(xenium.obj.stroke, subset = region %in% c("contralateral_stroke","stroke_site"))
xenium.obj.stroke<-subset(xenium.obj.stroke, subset = condition %in% c("Control"))

xenium.obj.stroke <- NormalizeData(xenium.obj.stroke)
xenium.obj.stroke <- FindVariableFeatures(xenium.obj.stroke, selection.method = "vst", nfeatures = 2000)
xenium.obj.stroke <- ScaleData(xenium.obj.stroke)
xenium.obj.stroke <- RunPCA(xenium.obj.stroke, npcs = 30, features = rownames(xenium.obj))
xenium.obj.stroke <- RunUMAP(xenium.obj.stroke, dims = 1:30)



umap_mat_myeloid<-as.data.frame(Embeddings(object = xenium.obj.stroke, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-xenium.obj.stroke@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_umap<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

len_x_bar<-((range(plt_umap$umap_1.y))[2]-(range(plt_umap$umap_1.y))[1])/10
len_y_bar<-((range(plt_umap$umap_2.y))[2]-(range(plt_umap$umap_2.y))[1])/10
arr <- list(x = min(plt_umap$umap_1.y)-1, y = min(plt_umap$umap_2.y)-1, x_len = len_x_bar, y_len = len_y_bar)


forlegned_plot<-ggplot(plt_umap, aes(umap_1.y,umap_2.y))+
  geom_point(aes(fill=CellType),size=3, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  fillscale_cellType+theme_bw()+
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=12))
nice_legend<-get_leg(forlegned_plot)

fanciest_UMAP<-ggplot(plt_umap, aes(umap_1.y,umap_2.y))+
  geom_point(size = 0.25, colour= "black", stroke = 1)+
  geom_point(aes(color=CellType),size=0.2)+xlab("UMAP 1")+ylab("UMAP 2")+
  colscale_cellType+
  annotate("segment",
           x = arr$x, xend = arr$x + c(arr$x_len, 0),
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=12,hjust = 0.05),
                     axis.title.y = element_text(size=12,hjust = 0.05,angle = 90),
                     legend.position = "none")

## Split
cell_num_all<-as.data.frame(table(xenium.obj.stroke@meta.data$region))
colnames(cell_num_all)<-c("region","CellCount")
fanciest_UMAP <- fanciest_UMAP + facet_wrap(~region, ncol=2)+  theme(strip.text = element_text(size=14))
geom_text(aes(x = (min(plt_umap$umap_1.y)+(0.95*len_x_bar))+0.5, y = (min(plt_umap$umap_2.y)+(0.5*len_y_bar))-1.5, label=paste0("n = ",comma(CellCount))), cell_num_all, size=7)

save_plts(plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2)), "UMAP_control_stroke_contralateral", w=14, h=6)




