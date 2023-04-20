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


source("scripts/00_pretty_plots.R")


load(here("data","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

d10x.combined_bcell<-subset(d10x.combined, subset = CellType_rough %in% c("B-cells"))
rm(d10x.combined)
gc()
d10x.combined_bcell <- RunPCA(d10x.combined_bcell, npcs = 30, verbose = FALSE)
d10x.combined_bcell <- RunUMAP(d10x.combined_bcell, reduction = "pca", dims = 1:30)
d10x.combined_bcell <- FindNeighbors(d10x.combined_bcell, reduction = "pca", dims = 1:30)
d10x.combined_bcell <- FindClusters(d10x.combined_bcell, resolution = 0.3)

DimPlot(d10x.combined_bcell, label=T)

bcell_cluster_umap<-DimPlot(d10x.combined_bcell, reduction = "umap", pt.size=0.25, label=T, group.by = "CellType_refined")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-14, y=-12, label=paste0("n = ",comma(ncol(d10x.combined_bcell))))
bcell_cluster_umap

cell_num_bell<-as.data.frame(table(d10x.combined_bcell$age_condition))
colnames(cell_num_bell)<-c("age_condition","CellCount")
bcell_cluster_umap<-DimPlot(d10x.combined_bcell, reduction = "umap", pt.size=0.25, label=F,split.by = "age_condition", group.by = "CellType_refined", ncol=2)+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  geom_text(aes(y=-14, x=-15,label=paste0("n = ",comma(CellCount))),cell_num_bell, hjust=-0.1, size=3)
bcell_cluster_umap
save_plts(bcell_cluster_umap, "IFALD_Bcell_map", w=7,h=6)

bcell_cluster_umap<-DimPlot(d10x.combined_bcell, reduction = "umap", pt.size=0.25, label=T,split.by = "age_id", group.by = "CellType_refined", ncol=4)+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-14, y=-12, label=paste0("n = ",comma(ncol(d10x.combined_bcell))))
bcell_cluster_umap



##############
## markers
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
d10x_raw_bcell<-subset(d10x, subset = CellType_rough %in% c("B-cells"))

identical(colnames(d10x_raw_bcell), colnames(d10x.combined_bcell))
d10x_raw_bcell <- AddMetaData(d10x_raw_bcell, metadata = d10x.combined_bcell@meta.data)

Idents(d10x_raw_bcell)<-d10x_raw_bcell$integrated_snn_res.0.3

de_4_0<-FindMarkers(d10x_raw_bcell, ident.1 = "4", ident.2 = "0", test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
sig_de_4_0<-de_4_0[which(de_4_0$p_val_adj < 0.005 & abs(de_4_0$avg_log2FC) > 1),]
sig_de_4_0[which(sig_de_4_0$avg_log2FC>0),]
sig_de_4_0[which(sig_de_4_0$avg_log2FC<0),]

de_6<-FindMarkers(d10x_raw_bcell, ident.1 = "6",  test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
sig_de_6<-de_6[which(de_6$p_val_adj < 0.005 & abs(de_6$avg_log2FC) > 1),]
head(sig_de_6)

b_markers<-FeaturePlot(d10x.combined_bcell, features = c("IRF8","IGLL1","MS4A1","CD79B"), min.cutoff = "q9", pt.size=1)
save_plts(b_markers, "IFALD_bcell_diff_genes", w=7,h=6)

# 4 is Pre Bcell? ^ maybe another type of a doublet?
# VPREB1
# https://www.ncbi.nlm.nih.gov/gene/7441
# Immature B cells: CHCHD10+, CD79a+, CD79b+, CD19+, MS4A1–/low, CD74–, Mki67+, Stmn1+
# https://www.cell.com/trends/immunology/fulltext/S1471-4906(22)00003-5
# IGLL1
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5017602/







#########
## pathway adult/IFALD versus healthy ped
#########
source("scripts/00_GSEA_function.R")
GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

de_4_0$gene<-rownames(de_4_0)
de<-de_4_0
gene_list = de$avg_log2FC
names(gene_list) = de$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

res = GSEA(gene_list, GO_file, pval = 0.05)

plt_path<-res$Results
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Enrichment_Cell<-"Up-regulated in \nIFALD"
plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \nHealthy"

plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))

plt_path$direction_label<-as.factor(plt_path$Enrichment)
levels(plt_path$direction_label)<-c(0.1,-0.1)
plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))

# top and bottom 15
#plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])

bcell_GSEA<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_hline(yintercept=16.5, color="grey")
save_plts(bcell_GSEA, "GSEA_IFALD_bcell", w=15,h=3)




###########
## Composition Bcell "types"
###########
d10x.combined_bcell@meta.data$CellType_rough<-as.character(d10x.combined_bcell@meta.data$CellType_rough)
d10x.combined_bcell@meta.data$CellType_rough[which(d10x.combined_bcell@meta.data$seurat_clusters%in%c("1","2","3"))]<-"Plasma"
d10x.combined_bcell@meta.data$CellType_rough[which(d10x.combined_bcell@meta.data$seurat_clusters%in%c("5"))]<-"Mature B cell (104)"
d10x.combined_bcell@meta.data$CellType_rough[which(d10x.combined_bcell@meta.data$seurat_clusters%in%c("4"))]<-"pre-B"
d10x.combined_bcell@meta.data$CellType_rough[which(d10x.combined_bcell@meta.data$seurat_clusters%in%c("0"))]<-"Mature B cell"
d10x.combined_bcell@meta.data$CellType_rough[which(d10x.combined_bcell@meta.data$seurat_clusters%in%c("6"))]<-"IRF8 B cell"


table(d10x.combined_bcell$CellType_rough, d10x.combined_bcell$age_condition)

cell_counts<-d10x.combined_bcell@meta.data %>% 
  group_by(individual, CellType_rough, age_condition,Age) %>% 
  summarise(count=length(unique(cell))) %>% 
  group_by(individual) %>%
  mutate(countT= sum(count)) %>%
  group_by(CellType_rough, add=TRUE) %>%
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
                                               "NA\n(IFALD006)", "NA\n(IFALD030)", "NA\n(IFALD073)", 
                                               "26\n(C82)","48\n(C70)", "57\n(C97)","61\n(C68)",
                                               "65\n(C39 NPC)", "65\n(C39 TLH)", "67\n(C54)","69\n(C88)"))

cell_counts_min<-cell_counts[,c("individual","age_condition","Age","countT","label")][!duplicated(cell_counts[,c("individual","age_condition","Age","countT","label")]),]

Bcell_composistion<-ggplot(cell_counts, aes(label, per))+geom_bar(aes(fill=CellType_rough),stat = "identity", color="black")+
  theme_bw()+th+scale_fill_manual(values=c("#2b103d","#9e5dc9","#e1ceed","#5f2585","#3e105c"), name="Cell Type")+xlab("Age\n(Sample ID)")+ylab("Percent of Cells in Sample")+
  facet_grid(.~age_condition, scale="free_x", space = "free")+
  geom_text(aes(label=countT, y=102), data=cell_counts_min)
Bcell_composistion
save_plts(Bcell_composistion, "Bcell_composistion", w=15,h=8)
