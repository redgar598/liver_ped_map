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

load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

d10x.combined_bcell<-subset(d10x.combined, subset = CellType_refined %in% c("Mature B-cells","Plasma cells","Cycling Plasma"))
rm(d10x.combined)
gc()
d10x.combined_bcell <- RunPCA(d10x.combined_bcell, npcs = 30, verbose = FALSE)
d10x.combined_bcell <- RunUMAP(d10x.combined_bcell, reduction = "pca", dims = 1:30)
d10x.combined_bcell <- FindNeighbors(d10x.combined_bcell, reduction = "pca", dims = 1:30)
d10x.combined_bcell <- FindClusters(d10x.combined_bcell, resolution = 0.3)

bcell_subtype<-DimPlot(d10x.combined_bcell, label=T)
save_plts(bcell_subtype, "IFALD_Bcell_map_clusters", w=7,h=6)

bcell_cluster_umap<-DimPlot(d10x.combined_bcell, reduction = "umap", pt.size=0.25, label=T, group.by = "CellType_refined")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-14, y=-12, label=paste0("n = ",comma(ncol(d10x.combined_bcell))))
bcell_cluster_umap

cell_num_bell<-as.data.frame(table(d10x.combined_bcell$age_condition))
colnames(cell_num_bell)<-c("age_condition","CellCount")
bcell_cluster_umap<-DimPlot(d10x.combined_bcell, reduction = "umap", pt.size=0.25, label=F,split.by = "age_condition", group.by = "CellType_refined", ncol=2)+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  geom_text(aes(y=-14, x=-15,label=paste0("n = ",comma(CellCount))),cell_num_bell, hjust=-0.1, size=3)
bcell_cluster_umap
save_plts(bcell_cluster_umap, "IFALD_Bcell_map", w=7,h=6)

bcell_cluster_umap<-DimPlot(d10x.combined_bcell, reduction = "umap", pt.size=0.25, label=F,split.by = "age_id", group.by = "CellType_refined", ncol=4)+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-14, y=-12, label=paste0("n = ",comma(ncol(d10x.combined_bcell))))
bcell_cluster_umap



## colored by age_condition
umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined_bcell, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined_bcell@meta.data
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

fanciest_UMAP<-ggplot(plt_myeloid[sample(nrow(plt_myeloid)),], aes(UMAP_1,UMAP_2))+
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
  annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar)-1, y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar)-2, label=paste0("n = ",comma(ncol(d10x.combined_bcell))), size=2)
save_plts(plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2)), "Healthy_and_IFALD_BCell_UMAP", w=3, h=2)





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


d10x_raw_bcell<-subset(d10x, subset = CellType_refined %in% c("Mature B-cells","Plasma cells","Cycling Plasma"))
rm(d10x)
gc()


identical(colnames(d10x_raw_bcell), colnames(d10x.combined_bcell))
d10x_raw_bcell <- AddMetaData(d10x_raw_bcell, metadata = d10x.combined_bcell@meta.data)

Idents(d10x_raw_bcell)<-d10x_raw_bcell$integrated_snn_res.0.3



de_4_0<-FindMarkers(d10x_raw_bcell, ident.1 = "4", ident.2 = "0", test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
sig_de_4_0<-de_4_0[which(de_4_0$p_val_adj < 0.005 & abs(de_4_0$avg_log2FC) > 1),]
sig_de_4_0[which(sig_de_4_0$avg_log2FC>0),]
sig_de_4_0[which(sig_de_4_0$avg_log2FC<0),]

de_5<-FindMarkers(d10x_raw_bcell, ident.1 = "5", ident.2 = "0",  test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
sig_de_5<-de_5[which(de_5$p_val_adj < 0.005 & abs(de_5$avg_log2FC) > 1),]
head(sig_de_5)
sig_de_5[which(sig_de_5$avg_log2FC > 1),]# high MT
FeaturePlot(d10x.combined_bcell, features = c("percent.mt","MT-ND4L","RPS26","MT-ATP8"), min.cutoff = "q9", pt.size=1)

de_6<-FindMarkers(d10x_raw_bcell, ident.1 = "6",  test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
sig_de_6<-de_6[which(de_6$p_val_adj < 0.005 & abs(de_6$avg_log2FC) > 1),]
head(sig_de_6)
sig_de_6[which(sig_de_6$avg_log2FC > 1),]


b_markers<-FeaturePlot(d10x.combined_bcell, features = c("IL3RA","PLAC8","IGLL1","MS4A1"), min.cutoff = "q9", pt.size=1)
save_plts(b_markers, "IFALD_bcell_diff_genes", w=7,h=6)

# 4 is Pre Bcell? ^ maybe another type of a doublet?
# VPREB1
# https://www.ncbi.nlm.nih.gov/gene/7441
# Immature B cells: CHCHD10+, CD79a+, CD79b+, CD19+, MS4A1–/low, CD74–, Mki67+, Stmn1+
# https://www.cell.com/trends/immunology/fulltext/S1471-4906(22)00003-5

# IGLL1
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5017602/
# PLAC8 and IL3RA are Plasmacytoid dendritic cell markers
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6861135/

## B cell type markers from:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6861135/
b_markers<-FeaturePlot(d10x.combined_bcell, features = c("JCHAIN","IGLL1","CD79B","TCL1A","IGKC","MS4A1","CD19","PLAC8","IL3RA"), min.cutoff = "q9", pt.size=1)
b_markers
save_plts(b_markers, "IFALD_bcell_markersfrom_fetalpaper", w=7,h=6)

DefaultAssay(d10x.combined_bcell) <- "RNA"
DotPlot(object = d10x.combined_bcell, features = c("JCHAIN","IGLL1","CD79B","TCL1A","IGKC","MS4A1","CD19"))+xlab("B Cell Marker")
DotPlot(object = d10x.combined_bcell, features = c("JCHAIN","IGKC","PLAC8","IL3RA","CD1C"))+xlab("B Cell Marker")
DefaultAssay(d10x.combined_bcell) <- "integrated"




###########
## Composition Bcell "types"
###########
d10x.combined_bcell@meta.data$CellType_refined<-as.character(d10x.combined_bcell@meta.data$CellType_refined)
d10x.combined_bcell@meta.data$CellType_refined[which(d10x.combined_bcell@meta.data$seurat_clusters%in%c("5"))]<-"Mature B-cells (High MT)"
d10x.combined_bcell@meta.data$CellType_refined[which(d10x.combined_bcell@meta.data$seurat_clusters%in%c("4"))]<-"pre B-cell"
d10x.combined_bcell@meta.data$CellType_refined[which(d10x.combined_bcell@meta.data$seurat_clusters%in%c("6"))]<-"pDC"

DimPlot(d10x.combined_bcell, reduction = "umap", pt.size=0.25, label=T, group.by = "CellType_refined")

cell_label_bcell<-d10x.combined_bcell@meta.data
save(cell_label_bcell, file=paste(here("data/"),"IFALD_B_cell_labels.rds", sep=""))





table(d10x.combined_bcell$CellType_refined, d10x.combined_bcell$age_condition)

cell_counts<-d10x.combined_bcell@meta.data %>% 
  group_by(individual, CellType_refined, age_condition,Age) %>% 
  summarise(count=length(unique(cell))) %>% 
  group_by(individual) %>%
  mutate(countT= sum(count)) %>%
  group_by(CellType_refined, add=TRUE) %>%
  mutate(per=100*count/countT)




count_plt<-as.data.frame(table(d10x.combined_bcell@meta.data$CellType_refined, d10x.combined_bcell$age_condition))


count_plt_percentage <- count_plt %>%
  group_by(Var1, Var2) %>%
  summarise(count = sum(Freq), .groups = 'drop') %>%
  group_by(Var1) %>%
  mutate(total_count = sum(count),
         percentage = (count / total_count) * 100) %>%
  ungroup()

count_percent<-ggplot(count_plt_percentage[which(count_plt_percentage$Var1%in%c("pDC","pre B-cell")),], aes(Var1,percentage))+
  geom_bar(aes(fill=Var2),stat="identity", color="black")+
  fillscale_agecondition+theme_bw()+
  xlab("")+ylab("Percent of Cells")+guides(fill=guide_legend(ncol=3))+coord_flip()+
  theme(legend.position="none")+
  theme(
    axis.text = element_text(size=12),
    strip.text = element_text(size=15),
    axis.title.y = element_text(size=15))
count_percent
save_plts(count_percent, "Bcell_bar_plot", w=4,h=1)





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

cell_counts_min<-cell_counts[,c("individual","age_condition","Age","countT","label")][!duplicated(cell_counts[,c("individual","age_condition","Age","countT","label")]),]


cell_counts$CellType_refined<-factor(cell_counts$CellType_refined, levels=c("pre B-cell","Mature B-cells","Mature B-cells (High MT)","Plasma cells", "pDC","Cycling Plasma"))
Bcell_composistion<-ggplot(cell_counts, aes(label, per))+geom_bar(aes(fill=CellType_refined),stat = "identity", color="black")+
  theme_bw()+th+fillscale_cellType+xlab("Age\n(Sample ID)")+ylab("Percent of Cells in Sample")+
  facet_grid(.~age_condition, scale="free_x", space = "free")+
  geom_text(aes(label=countT, y=102), data=cell_counts_min)
Bcell_composistion
save_plts(Bcell_composistion, "Bcell_composistion", w=19,h=10)


bcell_cluster_umap<-DimPlot(d10x.combined_bcell, reduction = "umap", pt.size=0.25, label=F, group.by = "CellType_refined")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=10, y=-15, label=paste0("n = ",comma(ncol(d10x.combined_bcell))))
bcell_cluster_umap
save_plts(bcell_cluster_umap, "IFALD_Bcell_map", w=7,h=5)

cell_num_bell<-as.data.frame(table(d10x.combined_bcell$age_condition))
colnames(cell_num_bell)<-c("age_condition","CellCount")
bcell_cluster_umap<-DimPlot(d10x.combined_bcell, reduction = "umap", pt.size=0.25, label=F,split.by = "age_condition", group.by = "CellType_refined", ncol=2)+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  geom_text(aes(x=7, y=-15,label=paste0("n = ",comma(CellCount))),cell_num_bell, hjust=-0.1, size=3)
bcell_cluster_umap
save_plts(bcell_cluster_umap, "IFALD_Bcell_map_split", w=8,h=6)

cell_num_bell<-as.data.frame(table(d10x.combined_bcell$age_id))
colnames(cell_num_bell)<-c("age_id","CellCount")
bcell_cluster_umap<-DimPlot(d10x.combined_bcell, reduction = "umap", pt.size=0.25, label=F,split.by = "age_id", group.by = "CellType_refined", ncol=4)+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  geom_text(aes(x=7, y=-15,label=paste0("n = ",comma(CellCount))),cell_num_bell, hjust=-0.1, size=3)
bcell_cluster_umap
save_plts(bcell_cluster_umap, "IFALD_Bcell_map_split_individual", w=12,h=10)


# fancy_bcell<-fanciest_UMAP(d10x.combined_bcell,"KC Like")
# save_plts(fancy_bcell, "IFALD_bcell_diff_cluster0_highlight", w=4,h=3)

fancy_bcell<-fanciest_UMAP(d10x.combined_bcell, NA,F)
fancy_bcell
save_plts(fancy_bcell, "IFALD_bcell_UMAP", w=3,h=2)

fancy_bcell<-fanciest_UMAP(d10x.combined_bcell, NA,T)
save_plts(fancy_bcell, "IFALD_bcell_UMAP_split", w=8,h=6)


###########
### plot individual genes
###########
b_markers<-plot_grid(plot_gene_UMAP(d10x.combined_bcell,"IL3RA", 0),
          plot_gene_UMAP(d10x.combined_bcell,"PLAC8", 0),
          plot_gene_UMAP(d10x.combined_bcell,"IGLL1", 0),
          plot_gene_UMAP(d10x.combined_bcell,"MS4A1", 0))
b_markers
save_plts(b_markers, "IFALD_bcell_diff_genes_bcellonly_tidy", w=7,h=5)

b_markers<-plot_grid(plot_gene_UMAP(d10x.combined_bcell,"IGLL1", 0),
                     plot_gene_UMAP(d10x.combined_bcell,"MS4A1", 0))
b_markers
save_plts(b_markers, "IFALD_prebcell_diff_genes_bcellonly_tidy", w=7,h=2.5)





####################################
## Fancy dot plot
####################################

DefaultAssay(d10x.combined_bcell)<-"RNA"

## markers
b_genes_noIG<-c("MS4A1","CD19","CD79B","TCL1A","IGLL1","JCHAIN","IGKC","IGHG1","MKI67","TOP2A","PLAC8","IL3RA")

gene_exp<-FetchData(d10x.combined_bcell, vars=c(b_genes_noIG))
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)

plt<-merge(gene_exp, d10x.combined_bcell@meta.data, by="cell")


## summarize
scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}

plt_summary<-plt %>% group_by(CellType_refined, variable) %>% 
  summarise(mn=mean(value), count=length(which(value>0)), percent_exp=(length(which(value>0))/length(value))*100)
plt_summary <- plt_summary %>% group_by(variable) %>%
  dplyr::mutate(scaled = scale_this(mn))
plt_summary<-as.data.frame(plt_summary)

# remove dots where 0 cell expressing marker
plt_summary<-plt_summary[(which(plt_summary$count>0)),]

plt_summary$variable<-factor(plt_summary$variable, levels=rev(c(b_genes_noIG)))

plt_summary$CellType_refined<-factor(plt_summary$CellType_refined, levels=(c("Mature B-cells","Mature B-cells (High MT)","pre B-cell",
                                                                             "Plasma cells","Cycling Plasma","pDC" )))


fancy_dotplot<-plot_grid(
  ggplot(plt_summary, aes(CellType_refined, variable, color=scaled, size=percent_exp))+geom_point()+
    th+theme_classic()+
    scale_color_continuous_sequential(palette = "Oslo", rev=F, name="Scaled\nMean\nExpression")+
    scale_size(name="Percent\nCells\nExpressing")+
    theme(axis.text.x = element_blank(),axis.title = element_blank(),axis.ticks.x = element_blank(),
          plot.margin = margin(0.25,0.25,0,0.25,"cm"))+
    geom_hline(yintercept = 10.5, color="grey70")+  
    geom_hline(yintercept = 8.5, color="grey70")+
    geom_hline(yintercept = 7.5, color="grey70")+
    geom_hline(yintercept = 4.5, color="grey70")+ 
    geom_hline(yintercept = 2.5, color="grey70")  ,
  ggplot(plt_summary, aes(CellType_refined, y=1, fill=CellType_refined))+geom_tile(color="black")+
    th+theme_classic()+fillscale_cellType+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",axis.line  = element_blank(),
          plot.margin = margin(0,0,1,1,"cm")),
  ncol=1, rel_heights = c(1.5,1), align = "v", axis="lr")
fancy_dotplot
save_plts(fancy_dotplot, "Healthy_only_dot_plot_celltype_Bcell", w=3,h=5)

