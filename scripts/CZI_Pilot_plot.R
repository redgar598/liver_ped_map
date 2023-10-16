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


d10x.combined_ped$FreshorFrozen<-"Fresh"
d10x.combined_ped$FreshorFrozen[which(d10x.combined_ped$individual%in%c("C85_caud3pr", "C93_caud3pr","C64_caud5pr","C96_caud3pr"))]<-"Frozen"
d10x.combined_ped$tissue_sampletype<-paste(d10x.combined_ped$Tissue, d10x.combined_ped$FreshorFrozen)

tapply(d10x.combined_ped@meta.data$individual, d10x.combined_ped@meta.data$tissue_sampletype, function(x) length(unique(x)))
d10x.combined_ped@meta.data$tissue_sampletype<-as.factor(d10x.combined_ped@meta.data$tissue_sampletype)
levels(d10x.combined_ped@meta.data$tissue_sampletype)<-c("Fresh Unperfused Biopsies\n(4 individuals)" ,  "Fresh Perfused Caudates\n(1 individuals)" , "Frozen Perfused Caudates\n(4 individuals)", "Fresh Perfused Right Lobe\n(1 individual)")


##################
## sample type
##################

umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined_ped, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined_ped@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

plt_myeloid<-plt_myeloid[which(plt_myeloid$tissue_sampletype!="Fresh Perfused Right Lobe\n(1 individual)"),]
plt_myeloid$tissue_sampletype<-as.factor(as.character(plt_myeloid$tissue_sampletype))
plt_myeloid$tissue_sampletype<-factor(plt_myeloid$tissue_sampletype, levels=c("Fresh Unperfused Biopsies\n(4 individuals)","Fresh Perfused Caudates\n(1 individuals)", "Frozen Perfused Caudates\n(4 individuals)" ))

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1), y = min(plt_myeloid$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)


forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=CellType_refined),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  fillscale_cellType+theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)

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
cell_num_all<-as.data.frame(table(plt_myeloid$tissue_sampletype))
colnames(cell_num_all)<-c("tissue_sampletype","CellCount")
fanciest_UMAP <- fanciest_UMAP + facet_wrap(~tissue_sampletype, ncol=4)+  geom_text(aes(x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(CellCount))), cell_num_all, size=2)

fanciest_UMAP <- fanciest_UMAP+annotate("rect",xmin = -8, xmax = 1.5,ymin = -10, ymax = 2,  color = "black",size=0.25, fill=NA)

fancy_type_split<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,1))
fancy_type_split
save_plts(fancy_type_split, "IFALD_liver_PBMC_split_pedonly_sampletype", w=10,h=3)



              

#####
## Tissue UMAP
####
umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined_ped, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined_ped@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

plt_myeloid$Tissue<-as.factor(plt_myeloid$Tissue)
levels(plt_myeloid$Tissue)<-c("Biopsy","Caudate","Right Lobe")
plt_myeloid$Tissue<-as.character(plt_myeloid$Tissue)

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1), y = min(plt_myeloid$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)

forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=Tissue),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  scale_fill_manual(values=c("cornflowerblue","grey","forestgreen"))+theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)

fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=Tissue),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  scale_color_manual(values=c("cornflowerblue","grey","forestgreen"))+
  annotate("segment",
           x = arr$x, xend = arr$x + c(arr$x_len, 0),
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")

fanciest_UMAP <- fanciest_UMAP + annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(nrow(plt_myeloid))), size=2)


fancy_type_tissue<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
save_plts(fancy_type_tissue, "IFALD_liver_sample_type_tissue", w=6,h=4)



## split by tissue and color by cell type
forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=CellType_refined),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  theme_bw()+fillscale_cellType+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)

fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=CellType_refined),size=0.01)+xlab("UMAP 1")+ylab("UMAP 2")+
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
cell_num_all<-as.data.frame(table(plt_myeloid$Tissue))
colnames(cell_num_all)<-c("Tissue","CellCount")
fanciest_UMAP <- fanciest_UMAP + facet_wrap(~Tissue, ncol=3)+  geom_text(aes(x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)-1, label=paste0("n = ",comma(CellCount))), cell_num_all, size=2)

fancy_type_split<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
fancy_type_split


fancy_type_tissue<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(8,2))
save_plts(fancy_type_tissue, "IFALD_liver_sample_type_celltype", w=12,h=3)





######################
## T cells only
######################
d10x.combined_ped_tcell<-subset(d10x.combined_ped, subset = CellType_refined %in% c("NK-like cells","CD3+ T-cells","CLNK T-cells","Cycling T-cells","gd T-cells"))
d10x.combined_ped_tcell <- RunPCA(d10x.combined_ped_tcell, npcs = 30, verbose = FALSE)
d10x.combined_ped_tcell <- RunUMAP(d10x.combined_ped_tcell, reduction = "pca", dims = 1:30)
d10x.combined_ped_tcell <- FindNeighbors(d10x.combined_ped_tcell, reduction = "pca", dims = 1:30)
d10x.combined_ped_tcell <- FindClusters(d10x.combined_ped_tcell, resolution = 0.1)

rm(d10x.combined_ped)
gc()

DimPlot(d10x.combined_ped_tcell, split.by = "tissue_sampletype")
DimPlot(d10x.combined_ped_tcell, split.by = "tissue_sampletype", group.by="CellType_refined", ncol=2) + colscale_cellType

## different proportion
table(d10x.combined_ped_tcell$CellType_refined, d10x.combined_ped_tcell$tissue_sampletype)

cell_counts<-d10x.combined_ped_tcell@meta.data %>% 
  group_by(individual, CellType_refined, tissue_sampletype) %>% 
  summarise(count=length(unique(cell))) %>% 
  group_by(individual) %>%
  mutate(countT= sum(count)) %>%
  group_by(CellType_refined, add=TRUE) %>%
  mutate(per=100*count/countT)

tcell_composistion<-ggplot(cell_counts, aes(individual, per))+geom_bar(aes(fill=CellType_refined),stat = "identity", color="black")+
  theme_bw()+th+fillscale_cellType+xlab("Individual ID")+ylab("Percent of Cells in Sample")+
  facet_grid(.~tissue_sampletype, scale="free_x", space = "free")
tcell_composistion


## fancy umap
umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined_ped_tcell, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined_ped_tcell@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

plt_myeloid<-plt_myeloid[which(plt_myeloid$tissue_sampletype!="Fresh Perfused Right Lobe\n(1 individual)"),]
plt_myeloid$tissue_sampletype<-as.factor(as.character(plt_myeloid$tissue_sampletype))
plt_myeloid$tissue_sampletype<-factor(plt_myeloid$tissue_sampletype, levels=c("Fresh Unperfused Biopsies\n(4 individuals)","Fresh Perfused Caudates\n(1 individuals)", "Frozen Perfused Caudates\n(4 individuals)" ))

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1), y = min(plt_myeloid$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)


forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=CellType_refined),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  fillscale_cellType+theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)

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
cell_num_all<-as.data.frame(table(plt_myeloid$tissue_sampletype))
colnames(cell_num_all)<-c("tissue_sampletype","CellCount")
fanciest_UMAP <- fanciest_UMAP + facet_wrap(~tissue_sampletype, ncol=4)+  geom_text(aes(x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(CellCount))), cell_num_all, size=2)

fanciest_UMAP <- fanciest_UMAP+
  annotate("rect",xmin = 7, xmax = 12,ymin = -3, ymax = 2,  color = "black",size=0.25, fill=NA)+
  annotate("rect",xmin = -9, xmax = -2,ymin = 5, ymax = 10,  color = "black",size=0.25, fill=NA)


fancy_type_split<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,1))
fancy_type_split
save_plts(fancy_type_split, "IFALD_liver_Tcells_split_pedonly_sampletype", w=10,h=3)
