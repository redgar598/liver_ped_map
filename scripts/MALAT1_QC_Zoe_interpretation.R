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


######################
## Compare filtered cells to original
######################

d10x.dropletQC<-readRDS(file = here("/media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024","IFALD_adult_ped_cellRough.rds"))
d10x.dropletQC$cell2<-sapply(1:nrow(d10x.dropletQC), function(x) paste(strsplit(d10x.dropletQC$cell[x], "-")[[1]][1], "-", d10x.dropletQC$individual[x], sep="" ))

load(here("/media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024","IFALD_QC_metrics.Rdata"))
plt_QC_data_noMALAT1<-plt_QC_data
plt_QC_data_noMALAT1$cell2<-sapply(1:nrow(plt_QC_data_noMALAT1), function(x) paste(strsplit(plt_QC_data_noMALAT1$cell[x], "-")[[1]][1], "-", plt_QC_data_noMALAT1$individual[x], sep="" ))

load(here("/media/redgar/Seagate\ Portable\ Drive/ped_map_update_feb2024","MALAT1_QC_metrics.Rdata"))
plt_QC_data_MALAT1<-plt_QC_data


dim(plt_QC_data_MALAT1)
dim(plt_QC_data_noMALAT1)

length(which(plt_QC_data_noMALAT1$cell2 %in% plt_QC_data_MALAT1$cell))



plt_QC_data_noMALAT1$MALAT1<-"Passed"
plt_QC_data_noMALAT1$MALAT1[which(!(plt_QC_data_noMALAT1$cell2 %in% plt_QC_data_MALAT1$cell))]<-"Filtered"
table(plt_QC_data_noMALAT1$MALAT1)


plt_QC_data_noMALAT1$MT_nFeature_filter<-"Passed"
plt_QC_data_noMALAT1$MT_nFeature_filter[which(!(plt_QC_data_noMALAT1$cell2 %in% d10x.dropletQC$cell2))]<-"Filtered"
table(plt_QC_data_noMALAT1$MT_nFeature_filter)

plt_QC_data_noMALAT1$Both_filtering<-sapply(1:nrow(plt_QC_data_noMALAT1), function(x){
  if(plt_QC_data_noMALAT1$MT_nFeature_filter[x]=="Passed" & plt_QC_data_noMALAT1$MALAT1[x]=="Passed"){"Passed\nin Both"}else{
    if(plt_QC_data_noMALAT1$MT_nFeature_filter[x]=="Filtered" & plt_QC_data_noMALAT1$MALAT1[x]=="Passed"){"Filtered by\nMT and nFeature\nfilter only"}else{
      if(plt_QC_data_noMALAT1$MT_nFeature_filter[x]=="Passed" & plt_QC_data_noMALAT1$MALAT1[x]=="Filtered"){"Filtered by\nMALAT1 only"}else{"Filtered\nby both"}
    }
  }
})

table(plt_QC_data_noMALAT1$Both_filtering)



QC_plot<-lapply(c("nCount_RNA","nFeature_RNA","percent.mt","relALBChange","nuclear_fraction"), function(y){
  ggplot(plt_QC_data_noMALAT1, aes_string("Both_filtering", y)) + geom_violin(fill="lightgrey", color="lightgrey")+geom_boxplot(width=0.1, outlier.shape = NA) +
    theme_bw()+th_present+xlab("")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
})

plot_grid(plotlist=QC_plot, ncol=2, axis="h")
save_plts(plot_grid(plotlist=QC_plot, ncol=2, axis="h"), "MALAT1_filtering_comparison", h=20, w=10)


table(plt_QC_data_noMALAT1$Both_filtering, plt_QC_data_noMALAT1$cell_status)






whodis<-d10x.dropletQC[which(d10x.dropletQC$cell2 %in% plt_QC_data_noMALAT1$cell2[which(plt_QC_data_noMALAT1$Both_filtering=="Filtered by\nMALAT1 only")]),]

round(table(whodis$CellType_rough)/nrow(whodis)*100, 2)
round(table(d10x.dropletQC$CellType_rough)/nrow(d10x.dropletQC)*100, 2)


round((table(whodis$CellType_rough)/table(d10x.dropletQC$CellType_rough))*100, 2)







d10x.combined<-readRDS(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024","IFALD_adult_ped_integrated.rds"))
d10x.combined@meta.data$CellType_rough<-as.factor(d10x.combined@meta.data$CellType_rough)
levels(d10x.combined@meta.data$CellType_rough)<-c("B-cells","Cholangiocytes",
                                                  "Hepatocytes",
                                                  "HSC","LSEC","Myeloid cells","NK and T cells")
load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/IFALD_adult_ped_SCINA_cell_labels.RData"))

length(which(SCINA_cell_labels$cell%in%colnames(d10x.combined)))

SCINA_cell_labels<-SCINA_cell_labels[match(colnames(d10x.combined), SCINA_cell_labels$cell),]
identical(colnames(d10x.combined), SCINA_cell_labels$cell)
rownames(SCINA_cell_labels)<-SCINA_cell_labels$cell
d10x.combined <- AddMetaData(d10x.combined, metadata = SCINA_cell_labels)

d10x.combined@meta.data$cell2<-sapply(1:nrow(d10x.combined@meta.data), function(x) paste(strsplit(d10x.combined@meta.data$cell[x], "-")[[1]][1], "-", d10x.combined@meta.data$individual[x], sep="" ))

length(which(d10x.combined@meta.data$cell2 %in% plt_QC_data_noMALAT1$cell2))

plt_QC_data_noMALAT1<-plt_QC_data_noMALAT1[which(plt_QC_data_noMALAT1$cell2 %in% d10x.combined@meta.data$cell2),]
plt_QC_data_noMALAT1<-plt_QC_data_noMALAT1[match(d10x.combined@meta.data$cell2, plt_QC_data_noMALAT1$cell2),]
identical(d10x.combined@meta.data$cell2, plt_QC_data_noMALAT1$cell2)
rownames(plt_QC_data_noMALAT1)<-plt_QC_data_noMALAT1$cell2

umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

plt_myeloid<-merge(plt_myeloid, plt_QC_data_noMALAT1, by="cell2")

cell_type_umap<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=CellType_rough),size=0.25)+xlab("UMAP 1")+ylab("UMAP 2")+
  colscale_cellType+theme_bw()+theme(legend.text = element_text(size=14),
                                     legend.title = element_text(size=16))+ 
  guides(colour = guide_legend(override.aes = list(size=3)))

plt_myeloid<-plt_myeloid[rev(order(plt_myeloid$Both_filtering)),]
MALAT1_filtered<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=Both_filtering),size=0.25)+xlab("UMAP 1")+ylab("UMAP 2")+
  theme_bw()+scale_color_manual(values=c("blue","grey"))+
  theme(legend.text = element_text(size=14),
        legend.title = element_text(size=16))+ 
  guides(colour = guide_legend(override.aes = list(size=3)))


save_plts(plot_grid(cell_type, MALAT1_filtered), "MALAT1_UMAP", w=20,h=8)


ggplot(plt_myeloid[which(plt_myeloid$Both_filtering=="Filtered by\nMALAT1 only"),], aes(UMAP_1,UMAP_2))+
  xlab("UMAP 1")+ylab("UMAP 2")+
  theme_bw()+
  theme(legend.text = element_text(size=14),
        legend.title = element_text(size=16))+ 
  geom_density_2d_filled(bins = 9)+scale_fill_brewer()+geom_point(size=0.25, alpha=0.1)
  
  
