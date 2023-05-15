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

d10x.combined_hsc<-subset(d10x.combined, subset = CellType_rough %in% c("HSC"))
rm(d10x.combined)
gc()
d10x.combined_hsc <- RunPCA(d10x.combined_hsc, npcs = 30, verbose = FALSE)
d10x.combined_hsc <- RunUMAP(d10x.combined_hsc, reduction = "pca", dims = 1:30)
d10x.combined_hsc <- FindNeighbors(d10x.combined_hsc, reduction = "pca", dims = 1:30)
d10x.combined_hsc <- FindClusters(d10x.combined_hsc, resolution = 0.1)

DimPlot(d10x.combined_hsc, label=T)

DimPlot(d10x.combined_hsc, reduction = "umap", pt.size=0.25, label=T, group.by = "CellType_refined")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")

cell_num_HSC<-as.data.frame(table(d10x.combined_hsc$age_condition))
colnames(cell_num_HSC)<-c("age_condition","CellCount")
HSC_map<-DimPlot(d10x.combined_hsc, reduction = "umap", pt.size=0.25, label=F,split.by = "age_condition", group.by = "CellType_refined", ncol=2)+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  geom_text(aes(y=-11, x=6,label=paste0("n = ",comma(CellCount))),cell_num_HSC, hjust=-0.1, size=3)
HSC_map
save_plts(HSC_map, "IFALD_HSC_map", w=7,h=6)

DimPlot(d10x.combined_hsc, reduction = "umap", pt.size=0.25, label=T,split.by = "age_id", group.by = "CellType_refined", ncol=4)+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")

d10x.combined_hsc@meta.data$CellType_rough<-as.character(d10x.combined_hsc@meta.data$CellType_rough)
d10x.combined_hsc@meta.data$CellType_rough[which(d10x.combined_hsc@meta.data$seurat_clusters%in%c("0"))]<-"healthy_ped_HSC"
d10x.combined_hsc@meta.data$CellType_rough[which(d10x.combined_hsc@meta.data$seurat_clusters%in%c("1","2"))]<-"adult_IFALD_HSC"
d10x.combined_hsc@meta.data$CellType_rough[which(d10x.combined_hsc@meta.data$seurat_clusters%in%c("3"))]<-"Outlier HSC"

fancy_HSC<-fanciest_UMAP(d10x.combined_hsc, NA,F)
fancy_HSC
save_plts(fancy_HSC, "IFALD_HSC_UMAP", w=4,h=3)

fancy_HSC<-fanciest_UMAP(d10x.combined_hsc, NA,T)
save_plts(fancy_HSC, "IFALD_HSC_UMAP_split", w=8,h=6)


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
d10x_raw_hsc<-subset(d10x, subset = CellType_rough %in% c("HSC"))

identical(colnames(d10x_raw_hsc), colnames(d10x.combined_hsc))
d10x_raw_hsc <- AddMetaData(d10x_raw_hsc, metadata = d10x.combined_hsc@meta.data)


Idents(d10x_raw_hsc)<-d10x_raw_hsc$CellType_rough

de_0<-FindMarkers(d10x_raw_hsc, ident.1 = "adult_IFALD_HSC", ident.2="healthy_ped_HSC", test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
de_0sig<-de_0[which(de_0$p_val_adj < 0.005 & abs(de_0$avg_log2FC) > 1),]
head(de_0sig[which(de_0sig$avg_log2FC>0),], n=10)
head(de_0sig[which(de_0sig$avg_log2FC<0),])

Fib_markers<-FeaturePlot(d10x.combined_hsc, features = c("PDGFRA","CXCL12","COL1A1","IGFBP3"), min.cutoff = "q9", pt.size=0.25)
save_plts(Fib_markers, "IFALD_HSC_diff_genes", w=7,h=6)

# Altogether our findings support a profibrotic role of PDGFR-α in HSCs
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7527859/

# Our data demonstrates a novel role of CXCL12 in stellate cell contraction
# contraction is healing abd fibrosis
# https://pubmed.ncbi.nlm.nih.gov/23812037/

# IGFBP3 role more complicated
# https://www.cmghjournal.org/article/S2352-345X(20)30095-3/fulltext

# COL1A1 a main major component in fibrotic tissues
# https://www.frontiersin.org/articles/10.3389/fcell.2021.765616/full#:~:text=Activated%20hepatic%20stellate%20cells%20(HSCs,major%20component%20in%20fibrotic%20tissues.


FeaturePlot(d10x.combined_hsc, features = c("MYH11","SPARCL1","ADIRF","MCAM"), min.cutoff = "q9", pt.size=0.25)


#########
## pathway adult/IFALD versus healthy ped
#########
source("scripts/00_GSEA_function.R")
GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

de_0$gene<-rownames(de_0)
de<-de_0
  gene_list = de$avg_log2FC
  names(gene_list) = de$gene
  gene_list = sort(gene_list, decreasing = TRUE)
  gene_list = gene_list[!duplicated(names(gene_list))]
  
  res = GSEA(gene_list, GO_file, pval = 0.05)
  
  plt_path<-res$Results
  plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
  plt_path$Enrichment_Cell<-"Up-regulated in \nAdult and IFALD"
  plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \n Healthy Pediatric"
  
  plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))
  
  plt_path$direction_label<-as.factor(plt_path$Enrichment)
  levels(plt_path$direction_label)<-c(0.1,-0.1)
  plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))
  
  # top and bottom 15
  plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])
  
  HSC_GSEA<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
    theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
    geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
    geom_hline(yintercept=16.5, color="grey")
  save_plts(HSC_GSEA, "GSEA_IFALD_HSC", w=15,h=6)


###########
## Composition fibrotic HSC
###########
table(d10x_raw_hsc$CellType_rough, d10x_raw_hsc$age_condition)

cell_counts<-d10x_raw_hsc@meta.data %>% 
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
cell_counts$CellType_rough<-factor(cell_counts$CellType_rough, c("healthy_ped_HSC","Outlier HSC","adult_IFALD_HSC"))
levels(cell_counts$CellType_rough)<-c("Healthy Ped HSC", "Outlier HSC" ,"Adult/IFALD HSC")
cell_counts_min<-cell_counts[,c("individual","age_condition","Age","countT","label")][!duplicated(cell_counts[,c("individual","age_condition","Age","countT","label")]),]
  
HSC_composistion<-ggplot(cell_counts, aes(label, per))+geom_bar(aes(fill=CellType_rough),stat = "identity", color="black")+
  theme_bw()+th+scale_fill_manual(values=c("#97bade","#5c7996","#11508f"), name="Cell Type")+xlab("Age\n(Sample ID)")+ylab("Percent of Cells in Sample")+
    facet_grid(.~age_condition, scale="free_x", space = "free")+
    geom_text(aes(label=countT, y=102), data=cell_counts_min)
HSC_composistion
save_plts(HSC_composistion, "HSC_composistion", w=15,h=8)

  