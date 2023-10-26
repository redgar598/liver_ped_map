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

d10x.combined_tcell<-subset(d10x.combined, subset = CellType_rough %in% c("NK and T cells"))
rm(d10x.combined)
gc()
d10x.combined_tcell <- RunPCA(d10x.combined_tcell, npcs = 30, verbose = FALSE)
d10x.combined_tcell <- RunUMAP(d10x.combined_tcell, reduction = "pca", dims = 1:30)
d10x.combined_tcell <- FindNeighbors(d10x.combined_tcell, reduction = "pca", dims = 1:30)
d10x.combined_tcell <- FindClusters(d10x.combined_tcell, resolution = 0.3)

tcell_subtype<-DimPlot(d10x.combined_tcell, label=T)
tcell_subtype

tcell_cluster_umap<-DimPlot(d10x.combined_tcell, reduction = "umap", pt.size=0.25, label=T, group.by = "CellType_refined")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-14, y=-12, label=paste0("n = ",comma(ncol(d10x.combined_tcell))))
tcell_cluster_umap

cell_num_bell<-as.data.frame(table(d10x.combined_tcell$age_condition))
colnames(cell_num_bell)<-c("age_condition","CellCount")
tcell_cluster_umap<-DimPlot(d10x.combined_tcell, reduction = "umap", pt.size=0.25, label=F,split.by = "age_condition", group.by = "CellType_refined", ncol=2)+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  geom_text(aes(y=-14, x=-15,label=paste0("n = ",comma(CellCount))),cell_num_bell, hjust=-0.1, size=3)
tcell_cluster_umap

tcell_cluster_umap<-DimPlot(d10x.combined_tcell, reduction = "umap", pt.size=0.25, label=F,split.by = "chemistry", group.by = "CellType_refined", ncol=4)+
  colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-14, y=-12, label=paste0("n = ",comma(ncol(d10x.combined_tcell))))
tcell_cluster_umap

fancy_tcell<-fanciest_UMAP(d10x.combined_tcell, NA,T)

##############
## DGE
##############
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))

######
## add cell type labels
######
load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
d10x_raw_tcell<-subset(d10x, subset = CellType_rough %in% c("NK and T cells"))
rm(d10x)
gc()


Idents(d10x_raw_tcell)<-d10x_raw_tcell$chemistry

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x_raw_tcell <- NormalizeData(d10x_raw_tcell,scale.factor = 10000, normalization.method = "LogNormalize")

d10x_raw_tcell_gd<-subset(d10x_raw_tcell, subset = CellType_refined %in% c("gd T-cells"))
## ped only cause adult mostly 3'
d10x_raw_tcell_gd<-subset(d10x_raw_tcell_gd, subset = AgeGroup %in% c("Ped"))


table(d10x_raw_tcell_gd$chemistry)

## age differential
de<-FindMarkers(d10x_raw_tcell_gd, ident.1 = "3pr", ident.2 = "5pr", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]
sig_de[which(sig_de$avg_log2FC>0),]
sig_de[which(sig_de$avg_log2FC<0),]


#########
## pathway adult/IFALD versus healthy ped
#########
source("scripts/00_GSEA_function.R")
GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

de$gene<-rownames(de)
gene_list = de$avg_log2FC
names(gene_list) = de$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

res = GSEA(gene_list, GO_file, pval = 0.05)

plt_path<-res$Results
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Enrichment_Cell<-"Up-regulated in \n3pr"
plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \n5pr"

plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))

plt_path$direction_label<-as.factor(plt_path$Enrichment)
levels(plt_path$direction_label)<-c(0.1,-0.1)
plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))

# top and bottom 15
plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])

tcell_GSEA<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_hline(yintercept=16.5, color="grey")
tcell_GSEA



###########
## Composition tcell "types"
###########
d10x.combined_tcell@meta.data$CellType_refined<-as.character(d10x.combined_tcell@meta.data$CellType_refined)
d10x.combined_tcell@meta.data$CellType_refined[which(d10x.combined_tcell@meta.data$seurat_clusters%in%c("1","2","5"))]<-"Plasma cells"
d10x.combined_tcell@meta.data$CellType_refined[which(d10x.combined_tcell@meta.data$seurat_clusters%in%c("8"))]<-"Mature B-cells (104)"
d10x.combined_tcell@meta.data$CellType_refined[which(d10x.combined_tcell@meta.data$seurat_clusters%in%c("7"))]<-"pre B-cell"
d10x.combined_tcell@meta.data$CellType_refined[which(d10x.combined_tcell@meta.data$seurat_clusters%in%c("0"))]<-"Mature B-cells"
d10x.combined_tcell@meta.data$CellType_refined[which(d10x.combined_tcell@meta.data$seurat_clusters%in%c("10"))]<-"pDC"

DimPlot(d10x.combined_tcell, reduction = "umap", pt.size=0.25, label=T, group.by = "CellType_refined")


d10x.combined_tcell<-subset(d10x.combined_tcell, subset = CellType_refined %in% c("Mature B-cells","Mature B-cells (104)","Plasma cells","pDC","pre B-cell"))
d10x.combined_tcell <- RunPCA(d10x.combined_tcell, npcs = 30, verbose = FALSE)
d10x.combined_tcell <- RunUMAP(d10x.combined_tcell, reduction = "pca", dims = 1:30)


table(d10x.combined_tcell$CellType_refined, d10x.combined_tcell$age_condition)

cell_counts<-d10x.combined_tcell@meta.data %>% 
  group_by(individual, CellType_refined, age_condition,Age) %>% 
  summarise(count=length(unique(cell))) %>% 
  group_by(individual) %>%
  mutate(countT= sum(count)) %>%
  group_by(CellType_refined, add=TRUE) %>%
  mutate(per=100*count/countT)


cell_counts$label<-sapply(1:nrow(cell_counts), function(x){
  if(length(grep("NPC", cell_counts$individual[x]))==1){
    paste(cell_counts$Age[x], "\n(", strsplit(cell_counts$individual[x],"_")[[1]][1]," NPC)", sep="")
  }else{if(length(grep("TLH", cell_counts$individual[x])==1)){
    paste(cell_counts$Age[x], "\n(", strsplit(cell_counts$individual[x],"_")[[1]][1]," TLH)", sep="")
  }else{
    paste(cell_counts$Age[x], "\n(", strsplit(cell_counts$individual[x],"_")[[1]][1],")", sep="")}}
})

cell_counts$label<-factor(cell_counts$label, c("2\n(C104)","11\n(C85)","12\n(C93)", "17\n(C64)", "17\n(C96)",
                                               "2\n(IFALD073)", "3\n(IFALD030)", "13\n(IFALD006)", 
                                               "26\n(C82)","48\n(C70)", "57\n(C97)","61\n(C68)",
                                               "65\n(C39 NPC)", "65\n(C39 TLH)", "67\n(C54)","69\n(C88)"))

cell_counts_min<-cell_counts[,c("individual","age_condition","Age","countT","label")][!duplicated(cell_counts[,c("individual","age_condition","Age","countT","label")]),]


cell_counts$CellType_refined<-factor(cell_counts$CellType_refined, levels=c("pre B-cell","Mature B-cells","Mature B-cells (104)","Plasma cells", "pDC"))
tcell_composistion<-ggplot(cell_counts, aes(label, per))+geom_bar(aes(fill=CellType_refined),stat = "identity", color="black")+
  theme_bw()+th+fillscale_cellType+xlab("Age\n(Sample ID)")+ylab("Percent of Cells in Sample")+
  facet_grid(.~age_condition, scale="free_x", space = "free")+
  geom_text(aes(label=countT, y=102), data=cell_counts_min)
tcell_composistion
save_plts(tcell_composistion, "tcell_composistion", w=16,h=10)


tcell_cluster_umap<-DimPlot(d10x.combined_tcell, reduction = "umap", pt.size=0.25, label=F, group.by = "CellType_refined")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=10, y=-15, label=paste0("n = ",comma(ncol(d10x.combined_tcell))))
tcell_cluster_umap
save_plts(tcell_cluster_umap, "IFALD_tcell_map", w=7,h=5)

cell_num_bell<-as.data.frame(table(d10x.combined_tcell$age_condition))
colnames(cell_num_bell)<-c("age_condition","CellCount")
tcell_cluster_umap<-DimPlot(d10x.combined_tcell, reduction = "umap", pt.size=0.25, label=F,split.by = "age_condition", group.by = "CellType_refined", ncol=2)+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  geom_text(aes(x=7, y=-15,label=paste0("n = ",comma(CellCount))),cell_num_bell, hjust=-0.1, size=3)
tcell_cluster_umap
save_plts(tcell_cluster_umap, "IFALD_tcell_map_split", w=8,h=6)

cell_num_bell<-as.data.frame(table(d10x.combined_tcell$age_id))
colnames(cell_num_bell)<-c("age_id","CellCount")
tcell_cluster_umap<-DimPlot(d10x.combined_tcell, reduction = "umap", pt.size=0.25, label=F,split.by = "age_id", group.by = "CellType_refined", ncol=4)+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  geom_text(aes(x=7, y=-15,label=paste0("n = ",comma(CellCount))),cell_num_bell, hjust=-0.1, size=3)
tcell_cluster_umap
save_plts(tcell_cluster_umap, "IFALD_tcell_map_split_individual", w=12,h=10)


# fancy_tcell<-fanciest_UMAP(d10x.combined_tcell,"KC Like")
# save_plts(fancy_tcell, "IFALD_tcell_diff_cluster0_highlight", w=4,h=3)

fancy_tcell<-fanciest_UMAP(d10x.combined_tcell, NA,F)
fancy_tcell
save_plts(fancy_tcell, "IFALD_tcell_UMAP", w=4,h=3)

fancy_tcell<-fanciest_UMAP(d10x.combined_tcell, NA,T)
save_plts(fancy_tcell, "IFALD_tcell_UMAP_split", w=8,h=6)


###########
### plot individual genes
###########
b_markers<-plot_grid(plot_gene_UMAP(d10x.combined_tcell,"IL3RA", 0),
                     plot_gene_UMAP(d10x.combined_tcell,"PLAC8", 0),
                     plot_gene_UMAP(d10x.combined_tcell,"IGLL1", 0),
                     plot_gene_UMAP(d10x.combined_tcell,"MS4A1", 0))
b_markers
save_plts(b_markers, "IFALD_tcell_diff_genes_tcellonly_tidy", w=7,h=5)
