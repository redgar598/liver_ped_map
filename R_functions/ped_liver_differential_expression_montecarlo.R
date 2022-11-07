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

###########
## Monte carlo the DE
###########


d10x_adult<-subset(d10x, subset = AgeGroup == "Adult")
ncol(d10x_adult)
d10x_ped<-subset(d10x, subset = AgeGroup == "Ped")
ncol(d10x_ped)


### paralize
command_args <- commandArgs(trailingOnly = TRUE)
cell_type_indx <- as.numeric(command_args[1])
cell_type<-cell_types[cell_type_indx]

samp_num=10


#DE_monte_carlo<-lapply(cell_types, function(cell_type){
  
  contrasts_celltype<-contrasts_celltype_age[grep(cell_type, contrasts_celltype_age)]
  
  d10x_adult_celltype<-subset(d10x_adult, subset = CellType_rough == cell_type)
  ncol(d10x_adult_celltype)
  d10x_ped_celltype<-subset(d10x_ped, subset = CellType_rough == cell_type)
  ncol(d10x_ped_celltype)

  
  de_lists<-sapply(1:samp_num, function(x){
    set.seed(x)
    
    ## make downsampled
    if(ncol(d10x_adult_celltype)<ncol(d10x_ped_celltype)){
      ped_cells_random <- d10x_ped_celltype[, sample(colnames(d10x_ped_celltype), size = ncol(d10x_adult_celltype), replace=F)]
      d10_DE<-merge(d10x_adult_celltype, ped_cells_random)
    }else{
      adult_cells_random <- d10x_adult_celltype[, sample(colnames(d10x_adult_celltype), size = ncol(d10x_ped_celltype), replace=F)]
      d10_DE<-merge(d10x_ped_celltype, adult_cells_random)}
    
    ## run DE 
    de<-FindMarkers(d10_DE, ident.1 = contrasts_celltype[1], ident.2 = contrasts_celltype[2], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
    print(paste(contrasts_celltype[1],"vs", contrasts_celltype[2],":", nrow(de), sep=" "))
    de$gene<-rownames(de)
    rownames(de)<-NULL
    de<-de[,c(6,1:5)]
    de$cell.1<-contrasts_celltype[1]
    de$cell.2<-contrasts_celltype[2]
    
    de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]$gene
    })
  
  sig_gene_count<-unlist(de_lists)
  if(length(sig_gene_count)==0){NA}else{
  sig_gene_count<-as.data.frame(table(sig_gene_count))
  colnames(sig_gene_count)<-c("gene","sig_count")
  
  sig_gene_count$monte_carlo_sig<-sapply(1:nrow(sig_gene_count),function(x){
    1-(sig_gene_count$sig_count[x]+1)/(samp_num+1)
  })
  
  sig_gene_count$cell<-cell_type
  sig_gene_count}#})

DE_monte_carlo<-do.call(rbind, DE_monte_carlo)
DE_monte_carlo<-DE_monte_carlo[which(!(is.na(DE_monte_carlo$gene))),]

save(DE_monte_carlo, file=here("data",paste(cell_type,"adult_ped_diff_motecarlo.RData",sep="_")))



# 
# #################
# ## Look at some interesting markers
# #################
# load(here("data","adult_ped_diff_motecarlo.RData"))
# DE_monte_carlo_sig<-DE_monte_carlo[which(DE_monte_carlo$monte_carlo_sig<0.001),]
# 
# ## Signatures
# myeloid_immune_supressive<-c("CTSB","CD163","MS4A7","FOLR2","GPNMB","VSIG4","HMOX1","MSR1")
# inflammatory_macs<-c("CD74","HLA-DRA","TYROBP","C1QC","HLA-DPA1","HLA-DPB1","LYZ","S100A6")
# exhausted_tcells<-c("TOX","PDCD1","LAG3","TNFRSF9","CXCL13","ENTPD1","HAVCR2","CD38")
# 
# 
# DE_monte_carlo[which(DE_monte_carlo$gene%in%c("KLRG1", "B3GAT1", "CD69","ITGAE")),]
# DE_monte_carlo[which(DE_monte_carlo$gene%in%c("LYZ", "MARCO", "MRC1","PTPRC")),]
# 
# keygenes<-c("KLRG1", "B3GAT1", "CD69","ITGAE","LYZ", "MARCO", "MRC1","PTPRC")
# DE_monte_carlo[which(DE_monte_carlo$gene%in%keygenes),]
# 
# Idents(d10x) <- "AgeGroup"
# 
# d10x@meta.data$CellType_rough<-as.factor(d10x@meta.data$CellType_rough)
# levels(d10x@meta.data$CellType_rough)<-c("CD3+ T-cells","Cholangiocytes",
#                                                   "gd T-cells","Hepatocytes",
#                                                   "HSC","LSEC","Myeloid cells","NK-like cells")
# 
# diff_exp_all$cell<-as.factor(diff_exp_all$cell)
# levels(diff_exp_all$cell)<-c("CD3+ T-cells","Cholangiocytes",
#                                          "gd T-cells","Hepatocytes",
#                                          "HSC","LSEC","Myeloid cells","NK-like cells")
# 
# cell_types<-as.factor(cell_types)
# levels(cell_types)<-c("CD3+ T-cells","Cholangiocytes",
#                                          "gd T-cells","Hepatocytes",
#                                          "HSC","LSEC","Myeloid cells","NK-like cells")
# 
# plot_key_genes<-function(keygenes, label){
#   all_plots<-lapply(1:length(keygenes), function(y){
#     plots<-lapply(1:length(cell_types),function(x){
#       p<-VlnPlot(subset(d10x, subset = CellType_rough == cell_types[x]) , features = keygenes[y], pt.size = 0, log=T)
#       p<-if(length(grep(cell_types[x], diff_exp_all[which(diff_exp_all$gene==keygenes[y]),]$cell.1))!=0){
#         p+theme(plot.background = element_rect(color = "black",size = 2)) +fillscale_age +xlab("") + ylab("")+ theme(legend.position="none")}else{
#           p+fillscale_age +xlab("") + ylab("")+ theme(legend.position="none")
#         }
#       p})
#     plot_grid(plotlist = plots, ncol=1)})
# 
# 
#   label_blank<-lapply(1:length(cell_types), function(x){
#     ggplot()+geom_blank()+theme_void()+ggtitle(cell_types[x])+ theme(plot.title = element_text(hjust = 0.5,vjust = -30))  })
#   label_blank<-plot_grid(plotlist = label_blank, ncol=1)
# 
#   plot_grid(label_blank, plot_grid(plotlist=all_plots, ncol=length(keygenes)), rel_widths=c(0.1,1))
#   ggsave2(paste0(here("figures/"), label, "_adult_ped.pdf"), w=20,h=20)
#   ggsave2(paste0(here("figures/jpeg/"),label, "_adult_ped.jpeg"), w=20,h=20,bg="white")}
# 
# plot_key_genes(myeloid_immune_supressive, "myeloid_immune_supressive_montecarlo")
# plot_key_genes(inflammatory_macs, "inflammatory_macs_montecarlo")
# plot_key_genes(exhausted_tcells, "exhausted_tcells_montecarlo")
# 
# 
# 
# #########
# ## Plot DE Genes
# #########
# load(here("data","adult_ped_integrated.rds"))
# 
# d10x.combined_myeloid<-subset(d10x.combined, subset = CellType_rough %in% c("Myeloid"))
# d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)
# d10x.combined_myeloid <- RunUMAP(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
# 
# d10x.combined_NK_T_B<-subset(d10x.combined, subset = CellType_rough %in% c("CD3_Tcell","gdTcell","nkTcell"))
# d10x.combined_NK_T_B <- RunPCA(d10x.combined_NK_T_B, npcs = 30, verbose = FALSE)
# d10x.combined_NK_T_B <- RunUMAP(d10x.combined_NK_T_B, reduction = "pca", dims = 1:30)
# gc()
# 
# 
# umapgene<-FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = "SERPINA1", ncol = 2, split="AgeGroup")
# save_plts(umapgene, "SERPINA1_umap_tcell_agesplit", w=12,h=5)
# violingene<-VlnPlot(subset(d10x, subset = CellType_rough == c( "CD3+ T-cells", "NK-like cells",  "gd T-cells" )) , features = c("SERPINA1"), pt.size = 0, log=T, split.by = "AgeGroup",  group.by= "CellType_rough")+fillscale_age
# save_plts(violingene, "SERPINA1_violin_tcell_agesplit", w=4,h=4)
# 
# umapgene<-FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = "ADH1B", ncol = 2, split="AgeGroup")
# save_plts(umapgene, "ADH1B_umap_tcell_agesplit", w=12,h=5)
# violingene<-VlnPlot(subset(d10x, subset = CellType_rough == c( "CD3+ T-cells", "NK-like cells",  "gd T-cells" )) , features = c("ADH1B"), pt.size = 0, log=T, split.by = "AgeGroup",  group.by= "CellType_rough")+fillscale_age
# save_plts(violingene, "ADH1B_violin_tcell_agesplit", w=4,h=4)
# 
# umapgene<-FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = "SAA1", ncol = 2, split="AgeGroup")
# save_plts(umapgene, "SAA1_umap_tcell_agesplit", w=12,h=5)
# violingene<-VlnPlot(subset(d10x, subset = CellType_rough == c( "CD3+ T-cells", "NK-like cells",  "gd T-cells" )) , features = c("SAA1"), pt.size = 0, log=T, split.by = "AgeGroup",  group.by= "CellType_rough")+fillscale_age
# save_plts(violingene, "SAA1_violin_tcell_agesplit", w=4,h=4)
# 
# umapgene<-FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = "KLRG1", ncol = 2, split="AgeGroup")
# save_plts(umapgene, "KLRG1_umap_tcell_agesplit", w=12,h=5)
# violingene<-VlnPlot(subset(d10x, subset = CellType_rough == c( "CD3+ T-cells", "NK-like cells",  "gd T-cells" )) , features = c("KLRG1"), pt.size = 0, log=T, split.by = "AgeGroup",  group.by= "CellType_rough")+fillscale_age
# save_plts(violingene, "KLRG1_violin_tcell_agesplit", w=4,h=4)
# 
# 
# 
# umapgene<-FeaturePlot(d10x.combined_myeloid, reduction = "umap", features = "ALDOB", ncol = 2, split="AgeGroup")
# umapgene
# umapgene<-FeaturePlot(d10x.combined, reduction = "umap", features = "SAA1", ncol = 2, split="AgeGroup")
# umapgene
# # save_plts(umapgene, "KLRG1_umap_tcell_agesplit", w=12,h=5)
# # violingene<-VlnPlot(subset(d10x, subset = CellType_rough == c( "CD3+ T-cells", "NK-like cells",  "gd T-cells" )) , features = c("KLRG1"), pt.size = 0, log=T, split.by = "AgeGroup",  group.by= "CellType_rough")+fillscale_age
# # save_plts(violingene, "KLRG1_violin_tcell_agesplit", w=4,h=4)
# 
# 
# # #
# # #
# # #
# # # # ### Top DE genes
# # # diff_exp_all %>%
# # #   group_by(cell.1) %>%
# # #   top_n(n = 10, wt = abs(avg_log2FC)) -> top10
# # #
# # # top_DE<-as.data.frame(top10)
# # #
# # # Idents(d10x) <- d10x$AgeGroup
# # #
# # # label_blank<-lapply(1:length(cell_types), function(x){
# # #   ggplot()+geom_blank()+theme_void()+ggtitle(cell_types[x])+ theme(plot.title = element_text(hjust = 0.5,vjust = -30))  })
# # # label_blank<-plot_grid(plotlist = label_blank, ncol=1)
# # #
# # # plot_list_top<-lapply(1:length(cell_types), function(x){
# # #   plots <- VlnPlot(subset(d10x, subset = CellType_rough == cell_types[x]) , features = top_DE[grep(cell_types[x],top_DE$cell.1),"gene"], pt.size = 0, log=T)
# # #   plots <- lapply(X = plots, FUN = function(p) p + fillscale_age +xlab("")+ theme(plot.title = element_text(size = 15)))
# # #   plot_grid(plotlist = plots, nrow=1)})
# # # top_DE_plot<-plot_grid(plotlist = plot_list_top, nrow=length(cell_types))
# # #
# # # plot_grid(label_blank, top_DE_plot, rel_widths=c(0.1,1))
# # #
# # # ggsave2(here("figures", "TopDE_adult_ped.pdf"), w=20,h=20)
# # # ggsave2(here("figures/jpeg", "TopDE_adult_ped.jpeg"), w=20,h=20,bg="white")
# # #
