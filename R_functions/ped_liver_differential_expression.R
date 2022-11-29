#'---
#'title: scRNAseq Differential Expression
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---



#'### Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(here)
library(ggplot2)
library(reshape2)
library(gridExtra)
#library(limma)
library(cowplot)
library(gtools)
#library(ggsignif)


options(stringsAsFactors = FALSE)

source("R_functions/pretty_plots.R")



## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","d10x_adult_ped_raw.rds"))


######
## add cell type labels
######
load(here("data","adult_ped_cellRough.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)




##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")


## testing factor
d10x$cell_age<-paste(d10x$CellType_rough, d10x$AgeGroup, sep = "_")
Idents(d10x) <- "cell_age"

table(d10x$CellType_rough, d10x$AgeGroup)


#MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene,
#consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
#and linear regression for the continuous process (i.e., the expression level). 

cell_types<-unique(d10x$CellType_rough)

contrasts_celltype_age<-do.call(rbind,lapply(1:length(cell_types), function(x){
  combinations(n = 2, r = 2, v = d10x$cell_age[grep(cell_types[x],d10x$cell_age)], repeats.allowed = FALSE)}))

contrasts_celltype_age

nrow(contrasts_celltype_age)

# this is 179,178 tests across all comparisons (29,863 genes, 6 comparisons)

diff_exp_all<-lapply(1:nrow(contrasts_celltype_age), function(x){
  de<-FindMarkers(d10x, ident.1 = contrasts_celltype_age[x,1], ident.2 = contrasts_celltype_age[x,2], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
  print(paste(contrasts_celltype_age[x,1],"vs", contrasts_celltype_age[x,2],":", nrow(de), sep=" "))
  de$gene<-rownames(de)
  rownames(de)<-NULL
  de<-de[,c(6,1:5)]
  de$cell.1<-contrasts_celltype_age[x,1]
  de$cell.2<-contrasts_celltype_age[x,2]
  de})


diff_exp_all<-do.call(rbind, diff_exp_all)

save(diff_exp_all, file=here("data","adult_ped_diff_genes.RData"))
#load(file=here("data","adult_ped_diff_genes.RData"))

#################
## Look at some interesting markers
#################
diff_exp_all[which(diff_exp_all$gene%in%c("KLRG1", "B3GAT1", "CD69","ITGAE")),]
diff_exp_all[which(diff_exp_all$gene%in%c("LYZ", "MARCO", "MRC1","PTPRC")),]

keygenes<-c("KLRG1", "B3GAT1", "CD69","ITGAE","LYZ", "MARCO", "MRC1","PTPRC")
diff_exp_all[which(diff_exp_all$gene%in%keygenes),]

Idents(d10x) <- "AgeGroup"

all_plots<-lapply(1:length(keygenes), function(y){
plots<-lapply(1:length(cell_types),function(x){
p<-VlnPlot(subset(d10x, subset = CellType_rough == cell_types[x]) , features = keygenes[y], pt.size = 0, log=T)
p<-if(length(grep(cell_types[x], diff_exp_all[which(diff_exp_all$gene==keygenes[y]),]$cell.1))!=0){
  p+theme(plot.background = element_rect(color = "black",size = 2)) +fillscale_age +xlab("") + ylab("")+ theme(legend.position="none")}else{  
    p+fillscale_age +xlab("") + ylab("")+ theme(legend.position="none")
}
p})
plot_grid(plotlist = plots, ncol=1)})


label_blank<-lapply(1:length(cell_types), function(x){
  ggplot()+geom_blank()+theme_void()+ggtitle(cell_types[x])+ theme(plot.title = element_text(hjust = 0.5,vjust = -30))  })
label_blank<-plot_grid(plotlist = label_blank, ncol=1)

plot_grid(label_blank, plot_grid(plotlist=all_plots, ncol=length(keygenes)), rel_widths=c(0.1,1))
ggsave2(here("figures", "keyGenes_adult_ped.pdf"), w=20,h=20)
ggsave2(here("figures/jpeg", "keyGenes_adult_ped.jpeg"), w=20,h=20,bg="white")



# ### Top DE genes
diff_exp_all %>%
  group_by(cell.1) %>%
  top_n(n = 10, wt = abs(avg_log2FC)) -> top10

top_DE<-as.data.frame(top10)

Idents(d10x) <- d10x$AgeGroup

label_blank<-lapply(1:length(cell_types), function(x){
  ggplot()+geom_blank()+theme_void()+ggtitle(cell_types[x])+ theme(plot.title = element_text(hjust = 0.5,vjust = -30))  })
label_blank<-plot_grid(plotlist = label_blank, ncol=1)

plot_list_top<-lapply(1:length(cell_types), function(x){
  plots <- VlnPlot(subset(d10x, subset = CellType_rough == cell_types[x]) , features = top_DE[grep(cell_types[x],top_DE$cell.1),"gene"], pt.size = 0, log=T)
  plots <- lapply(X = plots, FUN = function(p) p + fillscale_age +xlab("")+ theme(plot.title = element_text(size = 15)))
  plot_grid(plotlist = plots, nrow=1)})
top_DE_plot<-plot_grid(plotlist = plot_list_top, nrow=length(cell_types))

plot_grid(label_blank, top_DE_plot, rel_widths=c(0.1,1))

ggsave2(here("figures", "TopDE_adult_ped.pdf"), w=20,h=20)
ggsave2(here("figures/jpeg", "TopDE_adult_ped.jpeg"), w=20,h=20,bg="white")

