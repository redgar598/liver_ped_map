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


source("scRNA_seq_scripts/00_pretty_plots.R")
source("scRNA_seq_scripts/00_entropy_d10x.R")
source("scRNA_seq_scripts/00_fanciest_UMAP.R")
source("scRNA_seq_scripts/00_plot_gene_exp.R")






###########
## Integrated
###########
d10x_PBMC_liver<-readRDS(file = here("../../../projects/macparland/RE/PediatricAdult/processed_data","IFALD_adult_ped_PBMC_integrated.rds"))
#d10x_PBMC_liver<-readRDS(file = here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_PBMC_integrated.rds"))

d10x_PBMC_liver$Tissue<-as.factor(d10x_PBMC_liver$Tissue)
levels(d10x_PBMC_liver$Tissue)<-c("Biopsy Caudate","Biopsy Right Lobe", "Perfused Caudate","Perfused Right Lobe", "PBMC" )
d10x_PBMC_liver$Tissue<-as.character(d10x_PBMC_liver$Tissue)

SCT_cluster_umap<-DimPlot(d10x_PBMC_liver, reduction = "umap", pt.size=0.25, label=T)
save_plts(SCT_cluster_umap, "IFALD_rPCA_cluster_umap_PBMC", w=6,h=4)

SCT_cluster_tsne<-DimPlot(d10x_PBMC_liver, reduction = "tsne", pt.size=0.25, label=T)
save_plts(SCT_cluster_tsne, "IFALD_rPCA_cluster_tsne_PBMC", w=6,h=4)


individual_umap_sct<-DimPlot(d10x_PBMC_liver, reduction = "umap", group.by = "Tissue", pt.size=0.25)+
  scale_color_manual(values=c("lightskyblue3","darkolivegreen3","firebrick3","slategray2","forestgreen"))
individual_umap_sct
save_plts(individual_umap_sct, "IFALD_individual_rPCA_UMAP_PBMC", w=8,h=5)


###########
## Add cell type
###########
load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
#load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

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



fancy_073<-fanciest_UMAP(d10x_PBMC_liver, NA,F)
fancy_073

fanciest_UMAP(d10x_PBMC_liver, "pDC",F)
fanciest_UMAP(d10x_PBMC_liver, "pre B-cell",F)



umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x_PBMC_liver, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x_PBMC_liver@meta.data
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

fanciest_UMAP <- fanciest_UMAP + annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(ncol(d10x_PBMC_liver))), size=2)


fancy_type<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
save_plts(fancy_type, "IFALD_all_liver_PBMC", w=6,h=4)





umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x_PBMC_liver, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x_PBMC_liver@meta.data
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
cell_num_all<-as.data.frame(table(d10x_PBMC_liver@meta.data$individual))
colnames(cell_num_all)<-c("individual","CellCount")
fanciest_UMAP <- fanciest_UMAP + facet_wrap(~individual, ncol=4)+  geom_text(aes(x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(CellCount))), cell_num_all, size=2)


fancy_type_split<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
save_plts(fancy_type_split, "IFALD_liver_PBMC_split", w=12,h=8)



#####
## Tissue
####
umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x_PBMC_liver, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x_PBMC_liver@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1), y = min(plt_myeloid$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)

forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=Tissue),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  scale_fill_manual(values=c("lightskyblue3","darkolivegreen3","firebrick3","slategray2","forestgreen"))+theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)
#guides(colour = guide_legend(override.aes = list(size=0.5),byrow = TRUE))

fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=Tissue),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  scale_color_manual(values=c("lightskyblue3","darkolivegreen3","firebrick3","slategray2","forestgreen"))+
  annotate("segment",
           x = arr$x, xend = arr$x + c(arr$x_len, 0),
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")

fanciest_UMAP <- fanciest_UMAP + annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(ncol(d10x_PBMC_liver))), size=2)


fancy_type_tissue<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
save_plts(fancy_type_tissue, "IFALD_liver_PBMC_tissue", w=6,h=4)

#########
## Split by tissue
#########
umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x_PBMC_liver, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x_PBMC_liver@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1), y = min(plt_myeloid$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)

plt_myeloid$Tissue<-factor(plt_myeloid$Tissue, levels=c("Biopsy Caudate","Biopsy Right Lobe", "Perfused Caudate","Perfused Right Lobe", "PBMC" ))

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
cell_num_all<-as.data.frame(table(d10x_PBMC_liver@meta.data$Tissue))
colnames(cell_num_all)<-c("Tissue","CellCount")
fanciest_UMAP <- fanciest_UMAP + facet_wrap(~Tissue, ncol=5)+  geom_text(aes(x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(CellCount))), cell_num_all, size=2)

fancy_type_tissue<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(6,1))
save_plts(fancy_type_tissue, "IFALD_liver_PBMC_tissue_split", w=22,h=4)



############
## Tissue just biopsies
############


sample_split_UMAP<-function(samples){
  d10x_PBMC_liver_biopsy<-subset(d10x_PBMC_liver, subset = individual %in% samples)
  umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x_PBMC_liver_biopsy, reduction = "umap"))#
  umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
  meta_myeloid<-d10x_PBMC_liver_biopsy@meta.data
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

  fanciest_UMAP <- fanciest_UMAP + annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(ncol(d10x_PBMC_liver_biopsy))), size=2)
  fancy_type_tissue<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
  save_plts(fancy_type_tissue, paste("IFALD_liver_PBMC_tissue_biopsies",samples[2], sep=""), w=6,h=4)}

sample_split_UMAP(c("IFALD073_PBMC","IFALD073"))
sample_split_UMAP(c("IFALD073_PBMC","IFALD030"))
sample_split_UMAP(c("IFALD073_PBMC","IFALD006"))
sample_split_UMAP(c("IFALD073_PBMC","C104_bx5pr"))


######
## Entrophy
######
plt_entropy_Tissue<-entropy_d10(d10x_PBMC_liver, "Tissue")
entropy_Tissue<-entropy_plt(plt_entropy_Tissue, "Tissue", d10x_PBMC_liver)
entropy_Tissue
save_plts(entropy_Tissue, "IFALD_entropy_Tissue_PBMC_allclusters", w=15,h=10)


##############
## B cells and paried PBMC clusters only
##############
d10x.combined_bcell<-subset(d10x_PBMC_liver, subset = CellType_refined %in% c("DC","Plasma cells","Cycling Plasma","pDC","Mature B-cells","Mature B-cells (High MT)","pre B-cell","Cycling B-cells"))
rm(d10x_PBMC_liver)
gc()
d10x.combined_bcell <- RunPCA(d10x.combined_bcell, npcs = 30, verbose = FALSE)
d10x.combined_bcell <- RunUMAP(d10x.combined_bcell, reduction = "pca", dims = 1:30)
d10x.combined_bcell <- FindNeighbors(d10x.combined_bcell, reduction = "pca", dims = 1:30)
d10x.combined_bcell <- FindClusters(d10x.combined_bcell, resolution = 0.3)

fancy_PBMC_bcell<-fanciest_UMAP(d10x.combined_bcell, NA,F)
save_plts(fancy_PBMC_bcell, "IFALD_liver_PBMC_bcell", w=6,h=4)





umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined_bcell, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined_bcell@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

#plt_myeloid$Tissue<-factor(plt_myeloid$Tissue, levels=c("Biopsy Caudate","Biopsy Right Lobe", "Perfused Caudate","Perfused Right Lobe", "PBMC" ))

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1), y = min(plt_myeloid$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)

forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=Tissue),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  scale_fill_manual(values=c("lightskyblue3","darkolivegreen3","firebrick3","slategray2","forestgreen"))+theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)
#guides(colour = guide_legend(override.aes = list(size=0.5),byrow = TRUE))

fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=Tissue),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  scale_color_manual(values=c("lightskyblue3","darkolivegreen3","firebrick3","slategray2","forestgreen"))+
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
save_plts(fancy_type_tissue, "IFALD_liver_PBMC_bcell_tissue", w=6,h=4)

save_plts(nice_legend, "legend_tissue", w=4,h=3)

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
save_plts(fancy_individual, "IFALD_individual_fancy_Bcells", w=15,h=10)





###########
### plot individual genes
###########
d10x.combined_bcell$age_condition<-paste(d10x.combined_bcell$AgeGroup, d10x.combined_bcell$Treatment, sep=" ")

b_markers<-plot_grid(plot_gene_UMAP(d10x.combined_bcell,"IL3RA", 0),
                     plot_gene_UMAP(d10x.combined_bcell,"PLAC8", 0),
                     plot_gene_UMAP(d10x.combined_bcell,"IGLL1", 0),
                     plot_gene_UMAP(d10x.combined_bcell,"MS4A1", 0))
b_markers
save_plts(b_markers, "IFALD_bcell_diff_genes_bcellonly_tidy_integrated_PBMC", w=7,h=5)
