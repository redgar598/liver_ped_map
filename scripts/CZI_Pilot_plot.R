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




source("scripts/00_pretty_plots.R")
source("scripts/00_entropy_d10x.R")
source("scripts/00_fanciest_UMAP.R")
source("scripts/00_plot_gene_exp.R")


load(here("data","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

#some IFLAD didn't get labelled biopsy
d10x.combined@meta.data$Tissue[which(d10x.combined@meta.data$individual%in%c("IFALD030", "IFALD006"))]<-"Biopsy"


d10x.combined_ped<-subset(d10x.combined, subset = age_condition %in% c("Ped IFALD","Ped Healthy"))
rm(d10x.combined)
gc()





fanciest_UMAP(d10x.combined_ped, NA,F)



##################
## sample type
##################

umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined_ped, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined_ped@meta.data
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
cell_num_all<-as.data.frame(table(d10x.combined_ped@meta.data$Tissue))
colnames(cell_num_all)<-c("Tissue","CellCount")
fanciest_UMAP <- fanciest_UMAP + facet_wrap(~Tissue, ncol=4)+  geom_text(aes(x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(CellCount))), cell_num_all, size=2)

fanciest_UMAP <- fanciest_UMAP+annotate("rect",xmin = -8, xmax = 1.5,ymin = -10, ymax = 2,  color = "black",size=0.25, fill=NA)

fancy_type_split<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
fancy_type_split
save_plts(fancy_type_split, "IFALD_liver_PBMC_split_pedonly_sampletype", w=10,h=4)




###############
## Tissue
###############
d10x.combined$FreshorFrozen<-"Fresh"
d10x.combined$FreshorFrozen[which(d10x.combined$individual%in%c("C85_caud3pr", "C93_caud3pr","C64_caud5pr","C96_caud3pr"))]<-"Frozen"
d10x.combined$tissue_sampletype<-paste(d10x.combined$Tissue, d10x.combined$FreshorFrozen)

tapply(d10x.combined@meta.data$individual, d10x.combined@meta.data$tissue_sampletype, function(x) length(unique(x)))
d10x.combined@meta.data$tissue_sampletype<-as.factor(d10x.combined@meta.data$tissue_sampletype)
levels(d10x.combined@meta.data$tissue_sampletype)<-c("Fresh Unperfused Biopsies\n(4 individuals)" ,  "Fresh Perfused Caudates\n(8 individuals)" , "Frozen Perfused Caudates\n(4 individuals)")


umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1)-3, y = min(plt_myeloid$UMAP_2)-3, x_len = len_x_bar, y_len = len_y_bar)


forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=CellType_refined),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  fillscale_cellType+theme_bw()+ guides(fill = guide_legend(nrow = 2))+
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=8))
nice_legend<-get_leg(forlegned_plot)
#guides(colour = guide_legend(override.aes = list(size=0.5),byrow = TRUE))

fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.02, colour= "black", stroke = 0.8)+
  geom_point(aes(color=CellType_refined),size=0.01)+xlab("UMAP 1")+ylab("UMAP 2")+
  colscale_cellType+
  annotate("segment",
           x = arr$x, xend = arr$x + c(arr$x_len, 0),
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.01,angle = 90),
                     legend.position = "none")

## cell count
cell_num_all<-as.data.frame(table(d10x.combined@meta.data$tissue_sampletype))
colnames(cell_num_all)<-c("tissue_sampletype","CellCount")
fanciest_UMAP <- fanciest_UMAP + facet_wrap(~tissue_sampletype, ncol=1)+  geom_text(aes(x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)-1, label=paste0("n = ",comma(CellCount))), cell_num_all, size=2)

fanciest_UMAP <- fanciest_UMAP+annotate("rect",xmin = -8, xmax = 1.5,ymin = -10, ymax = 2,  color = "black",size=0.25, fill=NA)

fancy_type_split<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
fancy_type_split
save_plts(fanciest_UMAP, "IFALD_liver_PBMC_split_pedonly_tissueFrozen", w=2,h=6)
save_plts(nice_legend, "IFALD_liver_PBMC_split_pedonly_tissueFrozen_legend", w=12,h=1)

