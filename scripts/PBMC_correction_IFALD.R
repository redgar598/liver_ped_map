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


d10x_PBMC_liver<-readRDS(file = here("/media/redgar/Seagate Portable Drive/processed_data","IFALD_adult_ped_PBMC_integrated.rds"))

#some IFLAD didn't get labelled biopsy
d10x_PBMC_liver@meta.data$Tissue[which(d10x_PBMC_liver@meta.data$individual%in%c("IFALD030", "IFALD006"))]<-"Biopsy"

###########
## Visualize integration
###########
SCT_cluster_umap<-DimPlot(d10x_PBMC_liver, reduction = "umap", pt.size=0.25, label=T)
save_plts(SCT_cluster_umap, "IFALD_rPCA_cluster_umap_PBMC", w=6,h=4)

SCT_cluster_tsne<-DimPlot(d10x_PBMC_liver, reduction = "tsne", pt.size=0.25, label=T)
save_plts(SCT_cluster_tsne, "IFALD_rPCA_cluster_tsne_PBMC", w=6,h=4)


individual_umap_sct<-DimPlot(d10x_PBMC_liver, reduction = "umap", group.by = "Tissue", pt.size=0.5)+scale_color_manual(values=c("cornflowerblue","grey","red"))
save_plts(individual_umap_sct, "IFALD_individual_rPCA_UMAP_PBMC", w=6,h=5)


###########
## Add cell type
###########
load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
cell_label$index<-rownames(cell_label)

load(here("data","IFALD_B_cell_labels.rds"))
cell_label_bcell$index<-rownames(cell_label_bcell)
cell_label_notB<-cell_label[which(!(cell_label$index%in%cell_label_bcell$index)),]
cell_label_bcell<-cell_label_bcell[c("index","CellType_refined")]
cell_label_notB<-cell_label_notB[c("index","CellType_refined")]
cell_label<-rbind(cell_label_notB, cell_label_bcell)

load(here("data","IFALD_PBMC_cell_labels.rds"))
cell_label_PBMC$index<-rownames(cell_label_PBMC)
cell_label_PBMC<-cell_label_PBMC[c("index","CellType_refined")]

cell_label_integrated<-rbind(cell_label, cell_label_PBMC)
cell_label_integrated<-cell_label_integrated[match(colnames(d10x_PBMC_liver), cell_label_integrated$index),]
identical(colnames(d10x_PBMC_liver), cell_label_integrated$index)

d10x_PBMC_liver <- AddMetaData(d10x_PBMC_liver, metadata = cell_label_integrated)

DimPlot(d10x_PBMC_liver, reduction = "umap", group.by = "CellType_refined", pt.size=0.25)+colscale_cellType

## save 073 integration
saveRDS(d10x.combined, file = here("/media/redgar/Seagate Portable Drive/processed_data","IFALD_adult_ped_PBMC_integrated.rds"))


fancy_073<-fanciest_UMAP(d10x.combined, NA,F)
fancy_073

fanciest_UMAP(d10x.combined, "pDC",F)
fanciest_UMAP(d10x.combined, "pre B-cell",F)



umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1), y = min(plt_myeloid$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)


forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=CellType_refined),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  fillscale_cellType+theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)
#guides(colour = guide_legend(override.aes = list(size=0.5),byrow = TRUE))


plt_myeloid$highlight<-"0"
plt_myeloid$highlight[which(plt_myeloid$individual=="IFALD073_PBMC")]<-"1"

fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=CellType_refined),size=0.05)+
  geom_point(data=plt_myeloid[which(plt_myeloid$highlight==1),], size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=CellType_refined),data=plt_myeloid[which(plt_myeloid$highlight==1),], size=0.05)+
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

fanciest_UMAP <- fanciest_UMAP + annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(ncol(d10x.combined))), size=2)


fancy_type<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))    
save_plts(fancy_type, "IFALD_073_liver_PBMC", w=6,h=4)





umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1), y = min(plt_myeloid$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)


forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=CellType_refined),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  fillscale_cellType+theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)
#guides(colour = guide_legend(override.aes = list(size=0.5),byrow = TRUE))

fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=CellType_refined),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  colscale_cellType+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")

## cell count
cell_num_all<-as.data.frame(table(d10x.combined@meta.data$individual))
colnames(cell_num_all)<-c("individual","CellCount")
fanciest_UMAP <- fanciest_UMAP + facet_wrap(~individual, ncol=2)+  geom_text(aes(x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(CellCount))), cell_num_all, size=2)


fancy_type_split<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
save_plts(fancy_type_split, "IFALD_073_liver_PBMC_split", w=12,h=4)


