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
source("scripts/00_plot_gene_exp.R")



#############
## Merged Only
#############

#d10x_PBMC<-readRDS(file = here("/media/redgar/Seagate Portable Drive/processed_data","IFALD_d10x_adult_ped_raw_PBMC.rds"))
d10x_PBMC<-readRDS(file = here("../../../projects/macparland/RE/PediatricAdult/processed_data","IFALD_d10x_adult_ped_raw_PBMC.rds"))

d10x_liver<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))


d10x <- merge(d10x_liver, y= d10x_PBMC, merge.data=TRUE, project = "PBMC_liver_merged")#add.cell.ids = alldata_names2,
d10x

d10x <- NormalizeData(d10x)
d10x <- FindVariableFeatures(d10x, selection.method = "vst", nfeatures = 2000)
d10x <- ScaleData(d10x) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))

# dimension reduction
d10x <- RunPCA(d10x, ndims.print = 1:10, nfeatures.print = 10)
d10x <- RunUMAP(d10x, dims = 1:30)
d10x <- RunTSNE(d10x, dims = 1:30)

saveRDS(d10x, file = here("../../../projects/macparland/RE/PediatricAdult/processed_data","IFALD_adult_ped_PBMC_merged.rds"))

d10x<-readRDS(here("/media/redgar/Seagate Portable Drive/processed_data","IFALD_adult_ped_PBMC_merged.rds"))

#some IFLAD didn't get labelled biopsy
d10x@meta.data$Tissue[which(d10x@meta.data$individual%in%c("IFALD030", "IFALD006"))]<-"Biopsy"


SCT_cluster_umap<-DimPlot(d10x, reduction = "umap", pt.size=0.25, label=T)
save_plts(SCT_cluster_umap, "IFALD_rPCA_cluster_umap_PBMC_merged", w=6,h=4)

SCT_cluster_tsne<-DimPlot(d10x, reduction = "tsne", pt.size=0.25, label=T)
save_plts(SCT_cluster_tsne, "IFALD_rPCA_cluster_tsne_PBMC_merged", w=6,h=4)

individual_umap_sct<-DimPlot(d10x, reduction = "umap", group.by = "Tissue", pt.size=0.5)+scale_color_manual(values=c("cornflowerblue","grey","red"))
save_plts(individual_umap_sct, "IFALD_individual_rPCA_UMAP_PBMC_merged", w=6,h=5)


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
cell_label_integrated<-cell_label_integrated[match(colnames(d10x), cell_label_integrated$index),]
identical(colnames(d10x), cell_label_integrated$index)

d10x <- AddMetaData(d10x, metadata = cell_label_integrated)

DimPlot(d10x, reduction = "umap", group.by = "CellType_refined", pt.size=0.25)+colscale_cellType



fancy_073<-fanciest_UMAP(d10x, NA,F)
fancy_073

fanciest_UMAP(d10x, "pDC",F)
fanciest_UMAP(d10x, "pre B-cell",F)



umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x@meta.data
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

fanciest_UMAP <- fanciest_UMAP + annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(ncol(d10x))), size=2)


fancy_type<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
fancy_type
save_plts(fancy_type, "IFALD_all_liver_PBMC_merged", w=6,h=4)





umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x@meta.data
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
cell_num_all<-as.data.frame(table(d10x@meta.data$individual))
colnames(cell_num_all)<-c("individual","CellCount")
fanciest_UMAP <- fanciest_UMAP + facet_wrap(~individual, ncol=4)+  geom_text(aes(x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(CellCount))), cell_num_all, size=2)


fancy_type_split<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
fancy_type_split
save_plts(fancy_type_split, "IFALD_liver_PBMC_split_merged", w=12,h=8)



#####
## Tissue
####
umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1), y = min(plt_myeloid$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)

forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=Tissue),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  scale_fill_manual(values=c("cornflowerblue","grey","red"))+theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)
#guides(colour = guide_legend(override.aes = list(size=0.5),byrow = TRUE))

fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=Tissue),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  scale_color_manual(values=c("cornflowerblue","grey","red"))+
  annotate("segment",
           x = arr$x, xend = arr$x + c(arr$x_len, 0),
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")

fanciest_UMAP <- fanciest_UMAP + annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(ncol(d10x))), size=2)


fancy_type_tissue<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
fancy_type_tissue
save_plts(fancy_type_tissue, "IFALD_liver_PBMC_tissue_merged", w=6,h=4)

#########
## Split by tissue
#########
umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x@meta.data
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
cell_num_all<-as.data.frame(table(d10x@meta.data$Tissue))
colnames(cell_num_all)<-c("Tissue","CellCount")
fanciest_UMAP <- fanciest_UMAP + facet_wrap(~Tissue, ncol=3)+  geom_text(aes(x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(CellCount))), cell_num_all, size=2)

fancy_type_tissue<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(6,2))
fancy_type_tissue
save_plts(fancy_type_tissue, "IFALD_liver_PBMC_tissue_split_merged", w=14,h=4)



############
## Tissue just biopsies
############


sample_split_UMAP<-function(samples){
  d10x_biopsy<-subset(d10x, subset = individual %in% samples)
  umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x_biopsy, reduction = "umap"))#
  umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
  meta_myeloid<-d10x_biopsy@meta.data
  meta_myeloid$cell<-rownames(meta_myeloid)
  plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

  len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
  len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
  arr <- list(x = min(plt_myeloid$UMAP_1), y = min(plt_myeloid$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)

  forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
    geom_point(aes(fill=Tissue),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
    scale_fill_manual(values=c("grey","red"))+theme_bw()+
    theme(legend.text = element_text(size=5),
          legend.title = element_text(size=6))
  nice_legend<-get_leg(forlegned_plot)
  #guides(colour = guide_legend(override.aes = list(size=0.5),byrow = TRUE))

  fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
    geom_point(size = 0.06, colour= "black", stroke = 1)+
    geom_point(aes(color=Tissue),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
    scale_color_manual(values=c("grey","red"))+
    annotate("segment",
             x = arr$x, xend = arr$x + c(arr$x_len, 0),
             y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
             arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
    theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                       axis.title.x = element_text(size=5,hjust = 0.05),
                       axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                       legend.position = "none")

  fanciest_UMAP <- fanciest_UMAP + annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(ncol(d10x_biopsy))), size=2)
  fancy_type_tissue<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
  save_plts(fancy_type_tissue, paste("IFALD_liver_PBMC_tissue_biopsies_merged",samples[2], sep=""), w=6,h=4)}

sample_split_UMAP(c("IFALD073_PBMC","IFALD073"))
sample_split_UMAP(c("IFALD073_PBMC","IFALD030"))
sample_split_UMAP(c("IFALD073_PBMC","IFALD006"))
sample_split_UMAP(c("IFALD073_PBMC","C104_bx5pr"))





##############
## B cells and paried PBMC clusters only
##############

d10x.combined_bcell<-subset(d10x, subset = CellType_refined %in% c("DC","Naive CD4 T-cells","NK cells","Plasma cells","pDC","Mature B-cells","pre B-cell"))
rm(d10x)
gc()
d10x.combined_bcell <- RunPCA(d10x.combined_bcell, npcs = 30, verbose = FALSE)
d10x.combined_bcell <- RunUMAP(d10x.combined_bcell, reduction = "pca", dims = 1:30)
d10x.combined_bcell <- FindNeighbors(d10x.combined_bcell, reduction = "pca", dims = 1:30)
d10x.combined_bcell <- FindClusters(d10x.combined_bcell, resolution = 0.3)

fancy_PBMC_bcell<-fanciest_UMAP(d10x.combined_bcell, NA,F)
fancy_PBMC_bcell
save_plts(fancy_PBMC_bcell, "IFALD_liver_PBMC_bcell_merged", w=7,h=5)





umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined_bcell, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined_bcell@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1), y = min(plt_myeloid$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)

forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=Tissue),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  scale_fill_manual(values=c("grey","cornflowerblue","red"))+theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)
#guides(colour = guide_legend(override.aes = list(size=0.5),byrow = TRUE))

fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=Tissue),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  scale_color_manual(values=c("grey","cornflowerblue","red"))+
  annotate("segment",
           x = arr$x, xend = arr$x + c(arr$x_len, 0),
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")

fanciest_UMAP <- fanciest_UMAP + annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(ncol(d10x.combined_bcell))), size=2)


fancy_type_tissue<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
fancy_type_tissue
save_plts(fancy_type_tissue, "IFALD_liver_PBMC_bcell_tissue_merged", w=5,h=4)


#####################################
## fancy UMAP individual
#####################################
d10x.combined_bcell$age_id<-paste(d10x.combined_bcell$individual, d10x.combined_bcell$AgeGroup)

umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined_bcell, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined_bcell@meta.data
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
cell_num_all<-as.data.frame(table(d10x.combined_bcell@meta.data$age_id))
colnames(cell_num_all)<-c("age_id","CellCount")
fanciest_UMAP <- fanciest_UMAP + facet_wrap(~age_id, ncol=4)+  geom_text(aes(x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(CellCount))), cell_num_all, size=2)



fancy_individual<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
fancy_individual
save_plts(fancy_individual, "IFALD_individual_fancy_Bcells_merged", w=15,h=10)


FeaturePlot(d10x.combined_bcell, features=c("IGLL1","MS4A1"), split.by = "Tissue")

###########
### plot individual genes
###########
d10x.combined_bcell$age_condition<-paste(d10x.combined_bcell$AgeGroup, d10x.combined_bcell$Treatment, sep=" ")

b_markers<-plot_grid(plot_gene_UMAP(d10x.combined_bcell,"IL3RA", 0),
                     plot_gene_UMAP(d10x.combined_bcell,"PLAC8", 0),
                     plot_gene_UMAP(d10x.combined_bcell,"IGLL1", 0),
                     plot_gene_UMAP(d10x.combined_bcell,"MS4A1", 0))
b_markers
save_plts(b_markers, "IFALD_bcell_diff_genes_bcellonly_tidy_merged", w=7,h=5)
