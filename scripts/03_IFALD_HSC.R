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

fancyUMAP_all<-fanciest_UMAP(d10x.combined,"HSC",F)
save_plts(fancyUMAP_all, "IFALD_HSC_highlight_umap_fancy", w=6,h=4)

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
d10x.combined_hsc@meta.data$CellType_rough[which(d10x.combined_hsc@meta.data$seurat_clusters%in%c("0","2"))]<-"healthy_ped_HSC"
d10x.combined_hsc@meta.data$CellType_rough[which(d10x.combined_hsc@meta.data$seurat_clusters%in%c("1","3"))]<-"adult_IFALD_HSC"
d10x.combined_hsc@meta.data$CellType_rough[which(d10x.combined_hsc@meta.data$seurat_clusters%in%c("4"))]<-"Outlier HSC"

save(d10x.combined_hsc, file=here("data/HSC_integrated.RData"))



table(d10x.combined_hsc@meta.data$CellType_rough, d10x.combined_hsc$age_condition)
#Adult IFALD 70% Adult/IFLAD
#Healthy Ped 92% healthy ped


fancy_HSC<-fanciest_UMAP(d10x.combined_hsc, NA,F)
fancy_HSC
save_plts(fancy_HSC, "IFALD_HSC_UMAP", w=4,h=3)

fancy_HSC<-fanciest_UMAP(d10x.combined_hsc, NA,T)
save_plts(fancy_HSC, "IFALD_HSC_UMAP_split", w=8,h=6)


## colored by age_condition
umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined_hsc, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined_hsc@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1)-2, y = min(plt_myeloid$UMAP_2)-2, x_len = len_x_bar, y_len = len_y_bar)

forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=age_condition),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  fillscale_agecondition+theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)

fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=age_condition),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  colscale_agecondition+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")+
  annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar)-1, y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar)-2, label=paste0("n = ",comma(ncol(d10x.combined_hsc))), size=2)
save_plts(plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2)), "Healthy_and_IFALD_HSC_UMAP", w=3, h=2)


##############
## markers
##############
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_d10x_adult_ped_raw.rds"))

load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))


cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

d10x_raw_hsc<-subset(d10x, subset = CellType_rough %in% c("HSC"))
rm(d10x)
gc()

load(here("data/HSC_integrated.RData"))
identical(colnames(d10x_raw_hsc), colnames(d10x.combined_hsc))
d10x_raw_hsc <- AddMetaData(d10x_raw_hsc, metadata = d10x.combined_hsc@meta.data)


Idents(d10x_raw_hsc)<-d10x_raw_hsc$CellType_rough

de_0<-FindMarkers(d10x_raw_hsc, ident.1 = "adult_IFALD_HSC", ident.2="healthy_ped_HSC", test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
de_0sig<-de_0[which(de_0$p_val_adj < 0.005 & abs(de_0$avg_log2FC) > 1),]
head(de_0sig[which(de_0sig$avg_log2FC>0),], n=10)
head(de_0sig[which(de_0sig$avg_log2FC<0),])

Fib_markers<-FeaturePlot(d10x.combined_hsc, features = c("PDGFRA","CXCL12","COL1A1","IGFBP3"), min.cutoff = "q9", pt.size=0.25)
save_plts(Fib_markers, "IFALD_HSC_diff_genes", w=7,h=6)


### Diffential
write.csv(de_0sig, file=here("data","differential_HSC_types.csv"))



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

###########
### plot individual genes
###########
HSC_markers<-plot_grid(plot_gene_UMAP(d10x.combined_hsc,"PDGFRA", 0),
                     plot_gene_UMAP(d10x.combined_hsc,"CXCL12", 0),
                     plot_gene_UMAP(d10x.combined_hsc,"COL1A1", 0),
                     plot_gene_UMAP(d10x.combined_hsc,"IGFBP3", 0), ncol=4)
HSC_markers
save_plts(HSC_markers, "IFALD_HSC_diff_genes_fancy", w=14,h=2.5)

save_plts(plot_gene_UMAP(d10x.combined_hsc,"PDGFRA", 0), "IFALD_HSC_diff_PDGFRA_fancy", w=3,h=2.5)


#############
## Heat Map of differential genes
#############
gene<-c("PDGFRA","CXCL12","COL1A1","IGFBP3")

meta_HSC<-d10x.combined_hsc@meta.data
meta_HSC$cell<-rownames(meta_HSC)

gene_exp<-FetchData(d10x_raw_hsc, vars=gene)
gene_exp$cell<-rownames(gene_exp)

meta_HSC<-meta_HSC[,c("age_condition","CellType_rough","cell")]
plt_hsc<-merge(meta_HSC, gene_exp, by='cell')

melt_exp<-melt(plt_hsc)


### scale each genen across cells then take mean for gene in each cell type
scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}

plt_exp_summary <- melt_exp %>% 
  group_by(variable, CellType_rough) %>%
  dplyr::summarize(Mean = mean(value, na.rm=TRUE))

plt_exp_scaled <- plt_exp_summary %>% group_by(variable) %>%
  dplyr::mutate(scaled = scale_this(Mean))
plt_exp_summary<-as.data.frame(plt_exp_scaled)

plt_exp_summary$variable<-factor(plt_exp_scaled$variable, levels=rev(gene))

plt_exp_summary$CellType_rough<-as.factor(plt_exp_summary$CellType_rough)
levels(plt_exp_summary$CellType_rough)<-c("Adult \nand\nIFALD\nHSC","Healthy\nPed\nHSC","Outlier\nHSC")

de_0sig$variable<-rownames(de_0sig)
sig<-de_0sig[which(de_0sig$variable%in%gene),]
sig$label<-"*"
sig$CellType_rough<-"Adult \nand\nIFALD\nHSC"
HSC_fibrosis<-ggplot()+
  geom_tile(aes( CellType_rough,variable, fill=scaled), plt_exp_summary)+
  th+theme_classic()+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8, 'RdBu')), space='Lab')(100), name="Scaled\nMean\nExpression")+
  ylab("")+xlab("")+geom_text(data=sig, aes(CellType_rough,variable, label=label))
save_plts(HSC_fibrosis, "HSC_Fibrosis_heat", w=4.5,h=4)




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
  save_plts(HSC_GSEA, "GSEA_IFALD_HSC", w=20,h=6)


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

cell_counts$label<-factor(cell_counts$label, c("2\n(C104)","9\n(C105)","11\n(C85)","11\n(C115)","12\n(C93)", "16\n(C102)","16\n(C113)", "17\n(C64)", "17\n(C96)",
                                               "0.33\n(IFALD030)","0.58\n(IFALD073)", "9\n(IFALD006)", 
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
save_plts(HSC_composistion, "HSC_composistion", w=16,h=8)


#########################
## healthy age HSC
#########################
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_d10x_adult_ped_raw.rds"))

load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))


cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

d10x_raw_hsc<-subset(d10x, subset = CellType_rough %in% c("HSC"))


Idents(d10x_raw_hsc)<-d10x_raw_hsc$age_condition
table(d10x_raw_hsc$age_condition)


## age differential
de<-FindMarkers(d10x_raw_hsc, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]
sig_de[which(sig_de$avg_log2FC>0),]
sig_de[which(sig_de$avg_log2FC<0),]

write.csv(sig_de, file=here("data","differential_ageHSC.csv"))

### IFALD differential
de_IFALD<-FindMarkers(d10x_raw_hsc, ident.1 = "Ped IFALD", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de_IFALD<-de_IFALD[which(de_IFALD$p_val_adj < 0.005 & abs(de_IFALD$avg_log2FC) > 1),]
sig_de_IFALD[which(sig_de_IFALD$avg_log2FC>0),]
sig_de_IFALD[which(sig_de_IFALD$avg_log2FC<0),]

write.csv(sig_de_IFALD, file=here("data","differential_IFALDHSC.csv"))



############
#### CHECK claim that fibrosis would have been missed with only IFALD and adults
############
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_d10x_adult_ped_raw.rds"))

load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

d10x_raw_hsc<-subset(d10x, subset = CellType_rough %in% c("HSC"))
rm(d10x)
gc()


Idents(d10x_raw_hsc)<-d10x_raw_hsc$age_condition

de_0<-FindMarkers(d10x_raw_hsc, ident.1 = "Ped IFALD", ident.2="Adult Healthy", test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
de_0sig<-de_0[which(de_0$p_val_adj < 0.005 & abs(de_0$avg_log2FC) > 1),]
head(de_0sig[which(de_0sig$avg_log2FC>0),], n=10)
head(de_0sig[which(de_0sig$avg_log2FC<0),])

de_0sig[which(rownames(de_0sig) %in% c("PDGFRA","CXCL12","COL1A1","IGFBP3")),]
de_0[which(rownames(de_0) %in% c("PDGFRA","CXCL12","COL1A1","IGFBP3")),]


load(here("data/HSC_integrated.RData"))
d10x.combined_hsc_adult_IFALD<-subset(d10x.combined_hsc, subset = age_condition %in% c("Ped IFALD","Adult Healthy"))

DimPlot(d10x.combined_hsc_adult_IFALD, group.by = "age_condition")
Fib_markers<-FeaturePlot(d10x.combined_hsc_adult_IFALD, features = c("PDGFRA","CXCL12","COL1A1","IGFBP3"), min.cutoff = "q9", pt.size=0.25)
Fib_markers
VlnPlot(d10x_raw_hsc, features = c("PDGFRA","CXCL12","COL1A1","IGFBP3"))





################
# Check Identity of outlier HSC
################
Idents(d10x_raw_hsc)<-d10x_raw_hsc$CellType_rough

de_outlier<-FindMarkers(d10x_raw_hsc, ident.1 = "Outlier HSC", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
de_outlier_sig<-de_outlier[which(de_outlier$p_val_adj < 0.005 & abs(de_outlier$avg_log2FC) > 1),]
head(de_outlier_sig[which(de_outlier_sig$avg_log2FC>0),], n=10)
head(de_outlier_sig[which(de_outlier_sig$avg_log2FC<0),])

FeaturePlot(d10x.combined_hsc, features = c("LGI4","ERBB3","SOX10","S100B"), min.cutoff = "q9", pt.size=0.25)

HSC_markers<-plot_grid(plot_gene_UMAP(d10x.combined_hsc,"NRXN1", 0),
                       plot_gene_UMAP(d10x.combined_hsc,"ERBB3", 0),
                       plot_gene_UMAP(d10x.combined_hsc,"SOX10", 0),
                       plot_gene_UMAP(d10x.combined_hsc,"S100B", 0), ncol=4)
HSC_markers
save_plts(HSC_markers, "HSC_Glial_like_top_genes", w=14,h=2.5)

write.csv(de_outlier_sig[which(de_outlier_sig$avg_log2FC>0),], file="data/Glial_outlier_top_genes.csv")

de_outlier[grep("SPON|RSPO|ADAMTSL|NGFR|IGFBP3|SOX10",rownames(de_outlier)),]
FeaturePlot(d10x.combined_hsc, features = c("RSPO3","ADAMTSL2","NGFR","IGFBP3"), min.cutoff = "q9", pt.size=0.25)
VlnPlot(d10x_raw_hsc,features = c("IGFBP3","SOX10"))
VlnPlot(d10x_raw_hsc,features = c("AOX1","SOX10","NGFR"))

HSCs_genes<-c( "IGFBP7",  "SPARC")
FeaturePlot(d10x.combined_hsc, features = c("COL1A2",  "SPARC"), min.cutoff = "q9", pt.size=0.25)
FeaturePlot(d10x.combined_hsc, features = c("MYH11",  "RGS5", "LUM", "COL1A1","PLA2G2A"), min.cutoff = "q9", pt.size=0.25)

HSC_markers<-plot_grid(plot_gene_UMAP(d10x.combined_hsc,"COL1A2", 0),
                       plot_gene_UMAP(d10x.combined_hsc,"SPARC", 0),
                       plot_gene_UMAP(d10x.combined_hsc,"ID3", 0),
                       plot_gene_UMAP(d10x.combined_hsc,"IGFBP7", 0), ncol=4)
HSC_markers
save_plts(HSC_markers, "HSC_Glial_like_HSC_marker", w=14,h=2.5)



DimPlot(d10x.combined_hsc, group.by = "age_condition")
table(d10x.combined_hsc$CellType_rough, d10x.combined_hsc$individual)


# These ovelap with the genes in glial cells?
#   https://www.science.org/doi/10.1126/science.abo0510
# https://cellxgene.cziscience.com/e/d0c12af4-c0e4-4c7b-873a-70752b449689.cxg/
#   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3142460/
# schwann cells?

source("scripts/00_GSEA_function.R")
GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

de_outlier$gene<-rownames(de_outlier)
de<-de_outlier
gene_list = de$avg_log2FC
names(gene_list) = de$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

res = GSEA(gene_list, GO_file, pval = 0.05)

plt_path<-res$Results
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Enrichment_Cell<-"Up-regulated in \nOutlier"
plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \n General HSC"

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
save_plts(HSC_GSEA, "GSEA_outlier_HSC", w=20,h=6)



####################
## colored by outlier
####################
umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined_hsc, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined_hsc@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1)-2, y = min(plt_myeloid$UMAP_2)-2, x_len = len_x_bar, y_len = len_y_bar)

plt_myeloid$CellType_rough<-as.factor(plt_myeloid$CellType_rough)
levels(plt_myeloid$CellType_rough)<-c("Fibrotic HSC","HSC", "Glial Like")

forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=CellType_rough),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  scale_fill_manual(values=c("#2c7fb8","#7fcdbb","#d95f0e"))+
  theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)

fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=CellType_rough),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  scale_color_manual(values=c("#2c7fb8","#7fcdbb","#d95f0e"))+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")+
  annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar)-1, y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar)-2, label=paste0("n = ",comma(ncol(d10x.combined_hsc))), size=2)
plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
save_plts(plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2)), "Glial_HSC_UMAP", w=3, h=2)


