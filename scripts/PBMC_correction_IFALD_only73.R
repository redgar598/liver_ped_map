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


d10x_PBMC<-readRDS(file = here("/media/redgar/Seagate Portable Drive/processed_data","IFALD_d10x_adult_ped_raw_PBMC.rds"))
d10x_liver<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))

d10x_liver_IFALD073<-subset(d10x_liver, subset = individual %in% c("IFALD073"))

#############
## Merged Only
#############
d10x <- merge(d10x_liver_IFALD073, y= d10x_PBMC, merge.data=TRUE, project = "PBMC_liver073")#add.cell.ids = alldata_names2,
d10x

d10x <- NormalizeData(d10x)
d10x <- FindVariableFeatures(d10x, selection.method = "vst", nfeatures = 2000)
d10x <- ScaleData(d10x) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))

# dimension reduction
d10x <- RunPCA(d10x, ndims.print = 1:10, nfeatures.print = 10)
d10x <- RunUMAP(d10x, dims = 1:30)
d10x <- RunTSNE(d10x, dims = 1:30)


#############
## Integration
#############
d10x.list<- list(PBMC = d10x_PBMC, liver = d10x_liver_IFALD073)

## run integration across donor and hopefully that will also smooth out differences with chemistry?
## data is already split by donor

# normalize, identify variable features and score cell cycle for each dataset independently
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

d10x.list <- lapply(X = d10x.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = d10x.list)
d10x.list <- lapply(X = d10x.list, FUN = function(x) {
  #x <- ScaleData(x, features = features, verbose = FALSE)
  x <- ScaleData(x, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})



## Identify anchors
chem.anchors <- FindIntegrationAnchors(object.list = d10x.list, anchor.features = features, reduction = "rpca")
d10x.combined <- IntegrateData(anchorset = chem.anchors)

DefaultAssay(d10x.combined) <- "integrated"

print("INTEGRATED")


# Run the standard workflow for visualization and clustering
d10x.combined <- ScaleData(d10x.combined, verbose = FALSE)
d10x.combined <- RunPCA(d10x.combined, npcs = 30, verbose = FALSE)
d10x.combined <- RunUMAP(d10x.combined, reduction = "pca", dims = 1:30)
d10x.combined <- RunTSNE(d10x.combined, dims = 1:30)

d10x.combined <- FindNeighbors(d10x.combined, reduction = "pca", dims = 1:30)
d10x.combined <- FindClusters(d10x.combined, resolution = 0.5)

d10x.combined


###########
## Visualize integration
###########
SCT_cluster_umap<-DimPlot(d10x.combined, reduction = "umap", pt.size=0.25, label=T)
save_plts(SCT_cluster_umap, "IFALD_rPCA_cluster_umap_PBMC073_only", w=6,h=4)

SCT_cluster_tsne<-DimPlot(d10x.combined, reduction = "tsne", pt.size=0.25, label=T)
save_plts(SCT_cluster_tsne, "IFALD_rPCA_cluster_tsne_PBMC073_only", w=6,h=4)

individual_umap_sct<-DimPlot(d10x.combined, reduction = "umap", group.by = "individual", pt.size=0.25)
save_plts(individual_umap_sct, "IFALD_individual_rPCA_UMAP_PBMC073_only", w=6,h=4)


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

cell_label<-cell_label[match(colnames(d10x_liver_IFALD073), cell_label$index),]
identical(colnames(d10x_liver_IFALD073), cell_label$index)

load(here("data","IFALD_PBMC_cell_labels.rds"))
cell_label_PBMC$index<-rownames(cell_label_PBMC)
identical(colnames(d10x_PBMC), cell_label_PBMC$index)
cell_label_PBMC<-cell_label_PBMC[c("index","CellType_refined")]

cell_label_integrated<-rbind(cell_label, cell_label_PBMC)
cell_label_integrated<-cell_label_integrated[match(colnames(d10x.combined), cell_label_integrated$index),]
identical(colnames(d10x.combined), cell_label_integrated$index)

d10x.combined <- AddMetaData(d10x.combined, metadata = cell_label_integrated)

DimPlot(d10x.combined, reduction = "umap", group.by = "CellType_refined", pt.size=0.25)+colscale_cellType

## save 073 integration
saveRDS(d10x.combined, file = here("/media/redgar/Seagate Portable Drive/processed_data","IFALD_adult_ped_PBMC073_only_integrated.rds"))

d10x.combined<-readRDS(here("/media/redgar/Seagate Portable Drive/processed_data","IFALD_adult_ped_PBMC073_only_integrated.rds"))

fancy_073<-fanciest_UMAP(d10x.combined, NA,F,T)
fancy_073
save_plts(fancy_073, "IFALD_073_liver_PBMC", w=7,h=5)


fanciest_UMAP(d10x.combined, "pDC",F)
fanciest_UMAP(d10x.combined, "pre B-cell",F)


fanciest_UMAP(d10x.combined, "NK cells",F)
fanciest_UMAP(d10x.combined, "Naive CD4 T-cells",F)

    
    
#####
## Tissue
####
umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined@meta.data
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

fanciest_UMAP <- fanciest_UMAP + annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(ncol(d10x.combined))), size=2)


fancy_type_tissue<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
save_plts(fancy_type_tissue, "IFALD_073_liver_PBMC_tissue", w=5,h=4)

#########
## Split by tissue
#########
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
cell_num_all<-as.data.frame(table(d10x.combined@meta.data$Tissue))
colnames(cell_num_all)<-c("Tissue","CellCount")
fanciest_UMAP <- fanciest_UMAP + facet_wrap(~Tissue, ncol=2)+  geom_text(aes(x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(CellCount))), cell_num_all, size=2)

fancy_type_tissue<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
save_plts(fancy_type_tissue, "IFALD_073_liver_PBMC_tissue_split", w=10,h=4)


##############
## B cells and paried PBMC clusters only
##############

d10x.combined_bcell<-subset(d10x.combined, subset = CellType_refined %in% c("DC","Naive CD4 T-cells","NK cells","Plasma cells","pDC","Mature B-cells","pre B-cell"))
rm(d10x.combined)
gc()
d10x.combined_bcell <- RunPCA(d10x.combined_bcell, npcs = 30, verbose = FALSE)
d10x.combined_bcell <- RunUMAP(d10x.combined_bcell, reduction = "pca", dims = 1:30)
d10x.combined_bcell <- FindNeighbors(d10x.combined_bcell, reduction = "pca", dims = 1:30)
d10x.combined_bcell <- FindClusters(d10x.combined_bcell, resolution = 0.3)

fancy_PBMC_bcell<-fanciest_UMAP(d10x.combined_bcell, NA,F)
save_plts(fancy_PBMC_bcell, "IFALD_073_liver_PBMC_bcell", w=7,h=5)





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

fanciest_UMAP <- fanciest_UMAP + annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(ncol(d10x.combined_bcell))), size=2)


fancy_type_tissue<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
save_plts(fancy_type_tissue, "IFALD_073_liver_PBMC_bcell_tissue", w=5,h=4)
