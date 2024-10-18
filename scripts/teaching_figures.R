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


source("scripts/00_pretty_plots.R")
source("scripts/00_plot_gene_exp.R")
source("scripts/00_fanciest_UMAP.R")

load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1)-2, y = min(plt_myeloid$UMAP_2)-2, x_len = len_x_bar, y_len = len_y_bar)


forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=seurat_clusters),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)

fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=seurat_clusters),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  theme_classic()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     legend.position = "none")

## cell count
fanciest_UMAP <- fanciest_UMAP + annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar)-2, label=paste0("n = ",comma(ncol(d10x.combined))), size=4)

### cluster labels
plt_median<-plt_myeloid %>% group_by(seurat_clusters) %>% summarize(mean_umap1=median(UMAP_1), mean_umap2=median(UMAP_2))
plt_median<-as.data.frame(plt_median)

fanciest_UMAP <- fanciest_UMAP+geom_label(aes(mean_umap1, mean_umap2, label=seurat_clusters), data=plt_median,size=2, color="black")
fanciest_UMAP

save_plts(fanciest_UMAP, "teaching_UMAP_clusters", w=5, h=5)


## all grey
fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(color="grey",size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  theme_classic()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                        legend.position = "none")+
  annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar)-2, label=paste0("n = ",comma(ncol(d10x.combined))), size=4)

fanciest_UMAP

save_plts(fanciest_UMAP, "teaching_UMAP_clusters_grey", w=5, h=5)


#### by age group
forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=age_condition),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  theme_bw()+fillscale_agecondition+
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=10))
nice_legend<-get_leg(forlegned_plot)

fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+colscale_agecondition+
  geom_point(aes(color=age_condition),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  theme_classic()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                        legend.position = "none")+
  annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar)-2, label=paste0("n = ",comma(ncol(d10x.combined))), size=4)

fanciest_UMAP

save_plts(plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2)), "teaching_UMAP_age", w=7, h=5)




## Clusters not axis
fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=seurat_clusters),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=10,hjust = 0.05),
                     axis.title.y = element_text(size=10,hjust = 0.05,angle = 90),
                     legend.position = "none")+  
  annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar)-2, label=paste0("n = ",comma(ncol(d10x.combined))), size=4)+
  geom_label(aes(mean_umap1, mean_umap2, label=seurat_clusters), data=plt_median,size=2, color="black")

save_plts(fanciest_UMAP, "teaching_UMAP_axis", w=5, h=5)



## cell type

forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=CellType_refined),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  theme_bw()+fillscale_cellType+
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=10))
nice_legend<-get_leg(forlegned_plot)

fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=CellType_refined),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=10,hjust = 0.05),
                     axis.title.y = element_text(size=10,hjust = 0.05,angle = 90),
                     legend.position = "none")+  colscale_cellType+
  annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar)-2, label=paste0("n = ",comma(ncol(d10x.combined))), size=4)

save_plts(plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2)), "celltype_UMAP", w=7, h=5)






##################
## equivalent UMAP
##################

d10x.combined <- RunUMAP(d10x.combined, reduction = "pca", dims = 1:30, seed.use = 7, n.neighbors = 20L)

DimPlot(d10x.combined, group.by="seurat_clusters", label=T)


umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1)-2, y = min(plt_myeloid$UMAP_2)-2, x_len = len_x_bar, y_len = len_y_bar)

fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=seurat_clusters),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  theme_classic()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                        legend.position = "none")

fanciest_UMAP <- fanciest_UMAP + annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar)-2, label=paste0("n = ",comma(ncol(d10x.combined))), size=4)
plt_median<-plt_myeloid %>% group_by(seurat_clusters) %>% summarize(mean_umap1=median(UMAP_1), mean_umap2=median(UMAP_2))
plt_median<-as.data.frame(plt_median)
fanciest_UMAP <- fanciest_UMAP+geom_label(aes(mean_umap1, mean_umap2, label=seurat_clusters), data=plt_median,size=2, color="black")
fanciest_UMAP

save_plts(fanciest_UMAP, "teaching_UMAP_clusters_alternate", w=5, h=5)




############ dotplot


####################################
## Fancy dot plot
####################################

d10x<-readRDS(file =here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_d10x_adult_ped_raw.rds"))
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

## add cell type labels
load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

## markers
Macrophage_genes<-c( "PTPRC", "CD74")
LEC_genes<-c("CALCRL","RAMP2")
Hepatocyte_genes<-c("ALB", "CYP3A4")
Cholangiocytes_genes<-c( "EPCAM", "KRT7")
HSCs_genes<-c( "IGFBP7",  "SPARC")
T_genes<-c("CD3D","CD8A")
NK_genes<-c("NKG7","CD7")
CDC1_genes<-c("CLEC9A","XCR1")
gd_genes<-c("GNLY")
RBC<-c("HBA1","FCGR3A")
MAST<-c("TPSAB1", "AREG")
recent_recruit_myeloid<-c("S100A8","S100A9","CD68","LYZ")
kuffer_signature<-c("MARCO","VSIG4","CD5L")
neutro_gene<-c("CSF3R","FCGR3B")
MHCII<-c("HLA-DRA","HLA-DPB1")
b_genes_noIG<-c("MS4A1", "CD79B")
immunoglobins<-c("IGKC","IGHG1")
platelet_genes<-c("PPBP","NRGN")
cycle_genes<-c("MKI67","TOP2A")



gene_exp<-FetchData(d10x, vars=c(Macrophage_genes,LEC_genes,Hepatocyte_genes,Cholangiocytes_genes,HSCs_genes,T_genes,NK_genes,gd_genes,RBC,
                                 recent_recruit_myeloid, kuffer_signature, neutro_gene, MHCII,CDC1_genes, b_genes_noIG, immunoglobins, 
                                 platelet_genes,cycle_genes))
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)

plt<-merge(gene_exp, d10x@meta.data, by.x="cell", by.y="index")


## summarize
scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}

plt_summary<-plt %>% group_by(seurat_clusters, variable) %>% 
  summarise(mn=mean(value), count=length(which(value>0)), percent_exp=(length(which(value>0))/length(value))*100)
plt_summary <- plt_summary %>% group_by(variable) %>%
  dplyr::mutate(scaled = scale_this(mn))
plt_summary<-as.data.frame(plt_summary)

# remove dots where 0 cell expressing marker
plt_summary<-plt_summary[(which(plt_summary$count>0)),]

plt_summary$variable<-factor(plt_summary$variable, levels=rev(c(b_genes_noIG, immunoglobins, T_genes,gd_genes,NK_genes, 
                                                                Macrophage_genes,recent_recruit_myeloid,kuffer_signature,MHCII, CDC1_genes,
                                                                RBC, neutro_gene, platelet_genes,
                                                                Cholangiocytes_genes,LEC_genes,HSCs_genes,Hepatocyte_genes,cycle_genes)))


fancy_dotplot<-plot_grid(
  ggplot(plt_summary, aes(seurat_clusters, variable, color=scaled, size=percent_exp))+geom_point()+
    th+theme_classic()+
    scale_color_continuous_sequential(palette = "Oslo", rev=F, name="Scaled\nMean\nExpression")+
    scale_size(name="Percent\nCells\nExpressing")+
    theme(axis.text.x = element_blank(),axis.title = element_blank(),axis.ticks.x = element_blank(),
          plot.margin = margin(0.25,0.25,0,0.25,"cm"))+
    geom_hline(yintercept = 34.5, color="grey70")+ 
    geom_hline(yintercept = 32.5, color="grey70")+ 
    geom_hline(yintercept = 29.5, color="grey70")+ 
    geom_hline(yintercept = 27.5, color="grey70")+ 
    geom_hline(yintercept = 23.5, color="grey70")+ 
    geom_hline(yintercept = 20.5, color="grey70")+ 
    geom_hline(yintercept = 18.5, color="grey70")+ 
    geom_hline(yintercept = 16.5, color="grey70")+
    geom_hline(yintercept = 14.5, color="grey70")+
    geom_hline(yintercept = 12.5, color="grey70")+
    geom_hline(yintercept = 10.5, color="grey70")+  
    geom_hline(yintercept = 8.5, color="grey70")+
    geom_hline(yintercept = 6.5, color="grey70")+
    geom_hline(yintercept = 4.5, color="grey70")+ 
    geom_hline(yintercept = 2.5, color="grey70")  ,
  ggplot(plt_summary, aes(seurat_clusters, y=1, fill=seurat_clusters))+geom_tile(color="black")+
    th+theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",axis.line  = element_blank(),
          plot.margin = margin(0,0,1,1,"cm")),
  ncol=1, rel_heights = c(10,1), align = "v", axis="lr")
fancy_dotplot
save_plts(fancy_dotplot, "teaching_dot_plot_cluster", w=9,h=10)
