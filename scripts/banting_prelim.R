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
source("scripts/00_plot_gene_exp.R")
source("scripts/00_fanciest_UMAP.R")

load(here("data","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

d10x.combined_hsc<-subset(d10x.combined, subset = CellType_rough %in% c("HSC"))
rm(d10x.combined)
gc()
d10x.combined_hsc<-subset(d10x.combined_hsc, subset = Treatment %in% c("Healthy"))

d10x.combined_hsc <- RunPCA(d10x.combined_hsc, npcs = 30, verbose = FALSE)
d10x.combined_hsc <- RunUMAP(d10x.combined_hsc, reduction = "pca", dims = 1:30)
d10x.combined_hsc <- FindNeighbors(d10x.combined_hsc, reduction = "pca", dims = 1:30)
d10x.combined_hsc <- FindClusters(d10x.combined_hsc, resolution = 0.1)




d10x<-d10x.combined_hsc
umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1), y = min(plt_myeloid$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)


forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=AgeGroup),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  fillscale_age+theme_bw()+
  theme(legend.text = element_text(size=6),
        legend.title = element_text(size=7),
        legend.box.background = element_rect(colour = "black"))
nice_legend<-get_leg(forlegned_plot)

fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
    geom_point(size = 0.06, colour= "black", stroke = 1)+
    geom_point(aes(color=AgeGroup),size=0.05)+
  xlab("UMAP 1")+ylab("UMAP 2")+
    colscale_age+
    theme_void()+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),
                       legend.position = "none",
                       axis.title.x = element_text(size=6,hjust = 0.5, vjust = -2),
                       axis.title.y = element_text(size=6,hjust = 0.5,angle = 90, vjust = 2),
                       panel.border = element_rect(colour = "black", fill=NA))+
  annotate("text",x = min(plt_myeloid$UMAP_1)+(0.99*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.25*len_y_bar), label=paste0("n = ",comma(ncol(d10x.combined_hsc))), size=2)
fanciest_UMAP

save_plts( plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2)), "Healthy_only_HSC_UMAP", w=3,h=2)


genes<-c("CXCL12","COL1A1")
DefaultAssay(d10x) <- "RNA"
umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

cell_num_all<-as.data.frame(table(plt_myeloid$age_condition))
colnames(cell_num_all)<-c("age_condition","CellCount")

gene_exp<-FetchData(d10x, vars=genes)
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp, id="cell")
plt_myeloid<-merge(plt_myeloid, gene_exp, by='cell')

levels(plt_myeloid$age_condition)<-c("Pediatric","Pediatric","Adult")

HSC_differential_violin<-ggplot(plt_myeloid, aes(age_condition,log(value)))+
  geom_violin(fill="grey80",color="white")+geom_boxplot(aes(fill=age_condition),width=0.1)+scale_fill_manual(values=c("cornflowerblue","#D64A56"))+
  theme_bw()+th_present+xlab("Age Group")+ylab("Expression (log)")+
  theme(legend.position = "none")+facet_wrap(~variable)+
  theme(strip.background = element_rect(fill="white"))

save_plts(HSC_differential_violin, "Healthy_only_HSC_fibrosis_markers", w=5,h=3)



