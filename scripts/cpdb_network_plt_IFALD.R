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

load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

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





##############
## Differential expression with age
##############

## age differential
de_RR<-read.csv(here("data","differential_IFALD_RR.csv"))
de_KC<-read.csv(here("data","differential_IFALD_KC.csv"))
de_MHCII<-read.csv(here("data","differential_IFLAD_MHCII.csv"))


plt_median<-plt %>% group_by(CellType_refined) %>% summarize(mean_umap1=median(UMAP_1), mean_umap2=median(UMAP_2))
plt_median<-as.data.frame(plt_median)

###############
## Ped IFALD Output
###############
means<-read.table(here("data/cellphonedb/statistical_analysis_significant_means_03_05_2024_15:10:49.txt"), sep="\t", header=T)
means[1:5,1:10]


differentialexp_cpdbsig_a<-means[which(means$gene_a%in%de_KC$X), ]
differentialexp_cpdbsig_a$differential<-"adiff"
differentialexp_cpdbsig_b<-means[which(means$gene_b%in%de_KC$X), ]
differentialexp_cpdbsig_b$differential<-"bdiff"
differentialexp_cpdbsig<-rbind(differentialexp_cpdbsig_a, differentialexp_cpdbsig_b)

differentialexp_plt<-melt(differentialexp_cpdbsig, id=colnames(differentialexp_cpdbsig)[c(1:12,which(colnames(differentialexp_cpdbsig)=="differential"))])
differentialexp_plt$variable<-as.character(differentialexp_plt$variable)

differentialexp_plt$Cell1<-sapply(1:nrow(differentialexp_plt), function(x) strsplit(differentialexp_plt$variable[x], "[.]")[[1]][1])
differentialexp_plt$Cell2<-sapply(1:nrow(differentialexp_plt), function(x) strsplit(differentialexp_plt$variable[x], "[.]")[[1]][2])
differentialexp_plt<-differentialexp_plt[which(!(is.na(differentialexp_plt$value))),]



## fix cell labels
differentialexp_plt$Cell1<-as.factor(differentialexp_plt$Cell1)
levels(differentialexp_plt$Cell1)<-c("CD3+ T-cells","CDC1","Cholangiocytes", "Cycling Myeloid", "Cycling Plasma",
                                     "Doublet","Erythrocytes" ,"gd T-cells",  "Hepatocytes","HSC",
                                     "KC Like", "LSEC", 
                                     "Macrophage\n(MHCII high)","Mature B cells","Mono-Mac","Myeloid Erythrocytes\n(phagocytosis)",
                                     "Neutrophil", "NK-like cells","Plasma cells", "Platelets")
differentialexp_plt$Cell2<-as.factor(differentialexp_plt$Cell2)
levels(differentialexp_plt$Cell2)<-c("CD3+ T-cells","CDC1","Cholangiocytes", "Cycling Myeloid", "Cycling Plasma",
                                     "Doublet","Erythrocytes" ,"gd T-cells",  "Hepatocytes","HSC",
                                     "KC Like", "LSEC", 
                                     "Macrophage\n(MHCII high)","Mature B cells","Mono-Mac","Myeloid Erythrocytes\n(phagocytosis)",
                                     "Neutrophil", "NK-like cells","Plasma cells", "Platelets")

differentialexp_plt<-merge(differentialexp_plt,plt_median, by.x="Cell1", by.y="CellType_refined")
colnames(differentialexp_plt)[which(colnames(differentialexp_plt)%in%c("mean_umap1","mean_umap2"))]<-c("Cell1x","Cell1y")
differentialexp_plt<-merge(differentialexp_plt,plt_median, by.x="Cell2", by.y="CellType_refined")
colnames(differentialexp_plt)[which(colnames(differentialexp_plt)%in%c("mean_umap1","mean_umap2"))]<-c("Cell2x","Cell2y")


## Just sig in KC for now
#differentialexp_plt<-differentialexp_plt[grep("KC_Like",differentialexp_plt$variable),]
differentialexp_plta<-differentialexp_plt[which(differentialexp_plt$differential=="adiff" & differentialexp_plt$Cell1=="KC Like"),]
differentialexp_pltb<-differentialexp_plt[which(differentialexp_plt$differential=="bdiff" & differentialexp_plt$Cell2=="KC Like"),]
differentialexp_plt<-rbind(differentialexp_plta, differentialexp_pltb)



## self interactions
differentialexp_plt_self<-do.call(rbind,lapply(1:nrow(differentialexp_plt), function(x) if(differentialexp_plt$Cell1[x]==differentialexp_plt$Cell2[x]){differentialexp_plt[x,]}else{}))
differentialexp_plt_notself<-do.call(rbind,lapply(1:nrow(differentialexp_plt), function(x) if(differentialexp_plt$Cell1[x]==differentialexp_plt$Cell2[x]){}else{differentialexp_plt[x,]}))

unique(differentialexp_plt_notself$interacting_pair)


###############
## Color by receptor ligand expression
###############


save_plts(plot_gene_UMAP_2gene_network_notblend(d10x.combined,c("CXCL12","CXCR4"),"CXCL12_CXCR4",0.9), "UMAP_CXCL12_CXCR4_IFALD_map_celltype_network", w=9,h=6)
save_plts(plot_gene_UMAP_2gene_network_notblend(d10x.combined,c("THBS1","CD36"),"THBS1_CD36",0.9), "UMAP_THBS1_CD36_IFALD_map_celltype_network", w=9,h=6)
save_plts(plot_gene_UMAP_2gene_network_notblend(d10x.combined,c("CXCL12","DPP4"),"CXCL12_DPP4",0.9), "UMAP_CXCL12_DPP4_IFALD_map_celltype_network", w=9,h=6)

