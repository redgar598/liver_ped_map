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


source("scripts/00_pretty_plots.R")
source("scripts/00_fanciest_UMAP.R")
source("scripts/00_plot_gene_exp.R")

load(here("data","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

umap_mat<-as.data.frame(Embeddings(object = d10x.combined, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)
meta<-d10x.combined@meta.data
meta$cell<-rownames(meta)
plt<-merge(meta, umap_mat, by="cell")

plt_mean<-plt %>% group_by(CellType_refined) %>% summarize(mean_umap1=mean(UMAP_1), mean_umap2=mean(UMAP_2))
plt_mean<-as.data.frame(plt_mean)

len_x_bar<-((range(plt$UMAP_1))[2]-(range(plt$UMAP_1))[1])/10
len_y_bar<-((range(plt$UMAP_2))[2]-(range(plt$UMAP_2))[1])/10
arr <- list(x = min(plt$UMAP_1), y = min(plt$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)



###############
## Ped Healthy Output
###############
means<-read.table(here("data/cellphonedb/statistical_analysis_significant_means_06_22_2023_09:47:59.txt"), sep="\t", header=T)
means[1:5,1:10]

CCR_sig<-means[grep("CCR",means$interacting_pair), ]
CCL_sig<-means[grep("CCL3|CCL4",means$interacting_pair), ]
CCL_sig[,c(1:12)]

CCL_plt<-melt(CCL_sig, id=colnames(CCL_sig)[1:12])
CCL_plt$variable<-as.character(CCL_plt$variable)

CCL_plt$Cell1<-sapply(1:nrow(CCL_plt), function(x) strsplit(CCL_plt$variable[x], "[.]")[[1]][1])
CCL_plt$Cell2<-sapply(1:nrow(CCL_plt), function(x) strsplit(CCL_plt$variable[x], "[.]")[[1]][2])
CCL_plt<-CCL_plt[which(!(is.na(CCL_plt$value))),]

## fix cell labels
CCL_plt$Cell1<-as.factor(CCL_plt$Cell1)
levels(CCL_plt$Cell1)<-c("CD3+ T-cells", "CLNK T-cells", "Cycling Myeloid",
                         "Cycling T-cells","Doublet","gd T-cells",  
                         "KC Like","Low_Quality",  "Macrophage\n(CLEC9A high)",
                         "Macrophage\n(MHCII high)","Myeloid Erythrocytes\n(phagocytosis)",
                         "NK-like cells","Plasma cells", "Platelets","RR Myeloid")
CCL_plt$Cell2<-as.factor(CCL_plt$Cell2)
levels(CCL_plt$Cell2)<-c("CD3+ T-cells", "CLNK T-cells", "Cycling Myeloid",
                         "Cycling T-cells","Doublet", 
                         "KC Like","Low_Quality", 
                         "Macrophage\n(MHCII high)","Myeloid Erythrocytes\n(phagocytosis)",
                         "NK-like cells","RR Myeloid")


CCL_plt<-merge(CCL_plt,plt_mean, by.x="Cell1", by.y="CellType_refined")
colnames(CCL_plt)[which(colnames(CCL_plt)%in%c("mean_umap1","mean_umap2"))]<-c("Cell1x","Cell1y")
CCL_plt<-merge(CCL_plt,plt_mean, by.x="Cell2", by.y="CellType_refined")
colnames(CCL_plt)[which(colnames(CCL_plt)%in%c("mean_umap1","mean_umap2"))]<-c("Cell2x","Cell2y")


## self interactions
CCL_plt_self<-do.call(rbind,lapply(1:nrow(CCL_plt), function(x) if(CCL_plt$Cell1[x]==CCL_plt$Cell2[x]){CCL_plt[x,]}else{}))
CCL_plt_notself<-do.call(rbind,lapply(1:nrow(CCL_plt), function(x) if(CCL_plt$Cell1[x]==CCL_plt$Cell2[x]){}else{CCL_plt[x,]}))

## CCL only in KC
CCL_plt_notself_KC<-CCL_plt_notself[which(CCL_plt_notself$Cell1=="KC Like" & CCL_plt_notself$gene_a%in%c("CCL3","CCL4")),]

# Plotting the network
interacting_UMAP<-ggplot() +   
  annotate("segment", 
                      x = arr$x, xend = arr$x + c(arr$x_len, 0), 
                      y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
                      arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")+xlab("UMAP 1")+ylab("UMAP 2")+
  geom_point(aes(UMAP_1,UMAP_2), data=plt, size = 0.6, colour= "black", stroke = 1)+
  geom_point(aes(UMAP_1,UMAP_2, color=CellType_refined), data=plt,size=0.5)+xlab("UMAP 1")+ylab("UMAP 2")+
  #geom_rect(aes(xmin=range(plt$UMAP_1)[1]*1.1, xmax=range(plt$UMAP_1)[2]*1.1, ymin=range(plt$UMAP_2)[1]*1.1, ymax=range(plt$UMAP_2)[2]*1.1), fill="white", alpha=0.8) +
  geom_curve(
    data = CCL_plt_notself_KC[which(CCL_plt_notself_KC$interacting_pair=="CCL4_CCR5"),],
    aes(x = Cell1x, y = Cell1y, xend = Cell2x, yend = Cell2y), 
    color = "red",curvature = 0.2,
    lineend = "round") +
  geom_curve(
    data = CCL_plt_notself_KC[which(CCL_plt_notself_KC$interacting_pair=="CCL3_CCR5"),],
    aes(x = Cell1x, y = Cell1y, xend = Cell2x, yend = Cell2y), 
    color = "blue",curvature = 0.4,
    lineend = "round") +
  geom_curve(
    data = CCL_plt_notself_KC[which(CCL_plt_notself_KC$interacting_pair=="CCL3_CCR1"),],
    aes(x = Cell1x, y = Cell1y, xend = Cell2x, yend = Cell2y), 
    color = "grey",curvature = -0.4,
    lineend = "round") +
  geom_label(aes(mean_umap1, mean_umap2, label=CellType_refined, fill=CellType_refined), data=plt_mean, size=1.25, color="black")+
  geom_text(aes(x=range(plt$UMAP_1)[2]*0.9, y=range(plt$UMAP_2)[1]*0.9, label=unique(CCL_plt_notself_KC$interacting_pair)[1]), size=2, color="red")+
  geom_text(aes(x=range(plt$UMAP_1)[2]*0.9, y=range(plt$UMAP_2)[1]*0.95, label=unique(CCL_plt_notself_KC$interacting_pair)[2]), size=2, color="blue")+
  geom_text(aes(x=range(plt$UMAP_1)[2]*0.9, y=range(plt$UMAP_2)[1], label=unique(CCL_plt_notself_KC$interacting_pair)[3]), size=2, color="grey")+
  colscale_cellType+fillscale_cellType

save_plts(interacting_UMAP, "interacting_UMAP_KC_CCL34", w=5,h=5)




###################################################################################################################


##############
## Differential expression with age and IFALD in RR
##############
load(here("data","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

d10x.combined_myeloid<-subset(d10x.combined, subset = CellType_refined %in% c("RR Myeloid","Macrophage\n(MHCII high)","KC Like","Macrophage\n(CLEC9A high)","Cycling Myeloid","Myeloid Erythrocytes\n(phagocytosis)"))
rm(d10x.combined)
gc()
d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)
d10x.combined_myeloid <- RunUMAP(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindNeighbors(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindClusters(d10x.combined_myeloid, resolution = 0.3)
cluster0<-rownames(d10x.combined_myeloid@meta.data)[(d10x.combined_myeloid@meta.data$seurat_clusters==0)]


## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))

load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

d10x_raw_RR<-subset(d10x, subset = CellType_refined %in% c("RR Myeloid"))
d10x_raw_mhcII<-subset(d10x, subset = CellType_refined %in% c("Macrophage\n(MHCII high)"))

Idents(d10x_raw_RR)<-d10x_raw_RR$age_condition
Idents(d10x_raw_mhcII)<-d10x_raw_mhcII$age_condition

d10x$cellname<-colnames(d10x)
d10x_raw_myeloidcluster0 <- subset(d10x, subset = cellname %in% cluster0)

d10x_raw_myeloidcluster0$age_condition<-paste(d10x_raw_myeloidcluster0$AgeGroup, d10x_raw_myeloidcluster0$Treatment, sep=" ")
Idents(d10x_raw_myeloidcluster0)<-d10x_raw_myeloidcluster0$age_condition

## age differential
de_RR<-FindMarkers(d10x_raw_RR, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
de_KC<-FindMarkers(d10x_raw_myeloidcluster0, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
de_MHCII<-FindMarkers(d10x_raw_mhcII, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)


###############
## Ped Healthy Output
###############
means<-read.table(here("data/cellphonedb/statistical_analysis_significant_means_06_22_2023_09:47:59.txt"), sep="\t", header=T)
means[1:5,1:10]

differentialexp_cpdbsig<-means[which(means$gene_a%in%rownames(de_KC) | means$gene_b%in%rownames(de_KC)), ]
CCL_sig<-means[grep("CCL3|CCL4",means$interacting_pair), ]
CCL_sig[,c(1:12)]

CCL_plt<-melt(CCL_sig, id=colnames(CCL_sig)[1:12])
CCL_plt$variable<-as.character(CCL_plt$variable)

CCL_plt$Cell1<-sapply(1:nrow(CCL_plt), function(x) strsplit(CCL_plt$variable[x], "[.]")[[1]][1])
CCL_plt$Cell2<-sapply(1:nrow(CCL_plt), function(x) strsplit(CCL_plt$variable[x], "[.]")[[1]][2])
CCL_plt<-CCL_plt[which(!(is.na(CCL_plt$value))),]

## fix cell labels
CCL_plt$Cell1<-as.factor(CCL_plt$Cell1)
levels(CCL_plt$Cell1)<-c("CD3+ T-cells", "CLNK T-cells", "Cycling Myeloid",
                         "Cycling T-cells","Doublet","gd T-cells",  
                         "KC Like","Low_Quality",  "Macrophage\n(CLEC9A high)",
                         "Macrophage\n(MHCII high)","Myeloid Erythrocytes\n(phagocytosis)",
                         "NK-like cells","Plasma cells", "Platelets","RR Myeloid")
CCL_plt$Cell2<-as.factor(CCL_plt$Cell2)
levels(CCL_plt$Cell2)<-c("CD3+ T-cells", "CLNK T-cells", "Cycling Myeloid",
                         "Cycling T-cells","Doublet", 
                         "KC Like","Low_Quality", 
                         "Macrophage\n(MHCII high)","Myeloid Erythrocytes\n(phagocytosis)",
                         "NK-like cells","RR Myeloid")


CCL_plt<-merge(CCL_plt,plt_mean, by.x="Cell1", by.y="CellType_refined")
colnames(CCL_plt)[which(colnames(CCL_plt)%in%c("mean_umap1","mean_umap2"))]<-c("Cell1x","Cell1y")
CCL_plt<-merge(CCL_plt,plt_mean, by.x="Cell2", by.y="CellType_refined")
colnames(CCL_plt)[which(colnames(CCL_plt)%in%c("mean_umap1","mean_umap2"))]<-c("Cell2x","Cell2y")


## self interactions
CCL_plt_self<-do.call(rbind,lapply(1:nrow(CCL_plt), function(x) if(CCL_plt$Cell1[x]==CCL_plt$Cell2[x]){CCL_plt[x,]}else{}))
CCL_plt_notself<-do.call(rbind,lapply(1:nrow(CCL_plt), function(x) if(CCL_plt$Cell1[x]==CCL_plt$Cell2[x]){}else{CCL_plt[x,]}))

## CCL only in KC
CCL_plt_notself_KC<-CCL_plt_notself[which(CCL_plt_notself$Cell1=="KC Like" & CCL_plt_notself$gene_a%in%c("CCL3","CCL4")),]

# Plotting the network
interacting_UMAP<-ggplot() +   
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")+xlab("UMAP 1")+ylab("UMAP 2")+
  geom_point(aes(UMAP_1,UMAP_2), data=plt, size = 0.6, colour= "black", stroke = 1)+
  geom_point(aes(UMAP_1,UMAP_2, color=CellType_refined), data=plt,size=0.5)+xlab("UMAP 1")+ylab("UMAP 2")+
  #geom_rect(aes(xmin=range(plt$UMAP_1)[1]*1.1, xmax=range(plt$UMAP_1)[2]*1.1, ymin=range(plt$UMAP_2)[1]*1.1, ymax=range(plt$UMAP_2)[2]*1.1), fill="white", alpha=0.8) +
  geom_curve(
    data = CCL_plt_notself_KC[which(CCL_plt_notself_KC$interacting_pair=="CCL4_CCR5"),],
    aes(x = Cell1x, y = Cell1y, xend = Cell2x, yend = Cell2y), 
    color = "red",curvature = 0.2,
    lineend = "round") +
  geom_curve(
    data = CCL_plt_notself_KC[which(CCL_plt_notself_KC$interacting_pair=="CCL3_CCR5"),],
    aes(x = Cell1x, y = Cell1y, xend = Cell2x, yend = Cell2y), 
    color = "blue",curvature = 0.4,
    lineend = "round") +
  geom_curve(
    data = CCL_plt_notself_KC[which(CCL_plt_notself_KC$interacting_pair=="CCL3_CCR1"),],
    aes(x = Cell1x, y = Cell1y, xend = Cell2x, yend = Cell2y), 
    color = "grey",curvature = -0.4,
    lineend = "round") +
  geom_label(aes(mean_umap1, mean_umap2, label=CellType_refined, fill=CellType_refined), data=plt_mean, size=1.25, color="black")+
  geom_text(aes(x=range(plt$UMAP_1)[2]*0.9, y=range(plt$UMAP_2)[1]*0.9, label=unique(CCL_plt_notself_KC$interacting_pair)[1]), size=2, color="red")+
  geom_text(aes(x=range(plt$UMAP_1)[2]*0.9, y=range(plt$UMAP_2)[1]*0.95, label=unique(CCL_plt_notself_KC$interacting_pair)[2]), size=2, color="blue")+
  geom_text(aes(x=range(plt$UMAP_1)[2]*0.9, y=range(plt$UMAP_2)[1], label=unique(CCL_plt_notself_KC$interacting_pair)[3]), size=2, color="grey")+
  colscale_cellType+fillscale_cellType

save_plts(interacting_UMAP, "interacting_UMAP_KC_CCL34", w=5,h=5)


