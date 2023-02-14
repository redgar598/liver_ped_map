### Load libraries
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


load(here("data","adult_ped_integrated_refinedlabels_withDropletQC.rds"))

d10x.combined_myeloid<-subset(d10x.combined, subset = CellType_rough %in% c("Myeloid cells"))
rm(d10x.combined)
gc()
d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)
d10x.combined_myeloid <- RunUMAP(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindNeighbors(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindClusters(d10x.combined_myeloid, resolution = 0.1)
DimPlot(d10x.combined_myeloid, label=T)

myeliod_clusters<-d10x.combined_myeloid@meta.data

## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","d10x_adult_ped_raw.rds"))
load(here("data","adult_ped_cellRefined_withDropletQC.rds"))
cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)
d10x <- AddMetaData(d10x, metadata = cell_label)
# rm(d10x)
# gc()

d10x_raw_myeloid<-subset(d10x, subset = CellType_rough %in% c("Myeloid cells"))
d10x_raw_myeloid <- NormalizeData(d10x_raw_myeloid,scale.factor = 10000, normalization.method = "LogNormalize")


# 
# 
# ########## DGE
# Idents(d10x_raw_myeloid) <- "CellType_refined"
# 
# de<-FindMarkers(d10x_raw_myeloid, ident.1 = "RR Myeloid", ident.2 = "KC Like", test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
# sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]
# 
# de_nolatent<-FindMarkers(d10x_raw_myeloid, ident.1 = "RR Myeloid", ident.2 = "KC Like", test.use = "MAST", verbose=F)
# sig_de_nolatent<-de_nolatent[which(de_nolatent$p_val_adj < 0.005 & abs(de_nolatent$avg_log2FC) > 1),]
# 
# de_negbinom<-FindMarkers(d10x_raw_myeloid, ident.1 = "RR Myeloid", ident.2 = "KC Like", test.use = "negbinom",latent.vars="nFeature_RNA", verbose=F)
# sig_de_negbinom<-de_negbinom[which(de_negbinom$p_val_adj < 0.005 & abs(de_negbinom$avg_log2FC) > 1),]
# 
# de_negbinom_match<-de_negbinom[which(rownames(de_negbinom)%in%rownames(de)),]
# de_negbinom_match<-de_negbinom_match[match(rownames(de),rownames(de_negbinom_match)),]
# identical(rownames(de), rownames(de_negbinom_match))
# 
# 
# de_negbinom_nolatent<-FindMarkers(d10x_raw_myeloid, ident.1 = "RR Myeloid", ident.2 = "KC Like", test.use = "negbinom", verbose=F)
# sig_de_negbinom_nolatent<-de_negbinom_nolatent[which(de_negbinom_nolatent$p_val_adj < 0.005 & abs(de_negbinom_nolatent$avg_log2FC) > 1),]
# 
# de_negbinom<-de_negbinom[match(rownames(de_negbinom_nolatent),rownames(de_negbinom)),]
# identical(rownames(de_negbinom_nolatent), rownames(de_negbinom))
# identical(de_negbinom_nolatent$p_val, de_negbinom$p_val)
# cor(de_negbinom_nolatent$p_val, de_negbinom$p_val)
# 
# # de_DESeq2<-FindMarkers(d10x_raw_myeloid, ident.1 = "RR Myeloid", ident.2 = "KC Like", test.use = "DESeq2",latent.vars="nFeature_RNA", verbose=F)
# # sig_de_de_DESeq2<-de_DESeq2[which(de_DESeq2$p_val_adj < 0.005 & abs(de_DESeq2$avg_log2FC) > 1),]
# # 
# # de_DESeq2_nolatent<-FindMarkers(d10x_raw_myeloid, ident.1 = "RR Myeloid", ident.2 = "KC Like", test.use = "DESeq2", verbose=F)
# # sig_de_de_DESeq2_nolatent<-de_DESeq2_nolatent[which(de_DESeq2_nolatent$p_val_adj < 0.005 & abs(de_DESeq2_nolatent$avg_log2FC) > 1),]
# 
# save(de,de_nolatent, de_negbinom, de_negbinom_nolatent, file="../../Downloads/de_tmp.Rdata")
# 
# load("../../Downloads/de_tmp.Rdata")
# 
# colnames(de_negbinom)<-paste(colnames(de_negbinom),"negbinom", sep="_")
# de_negbinom$gene<-rownames(de_negbinom)
# colnames(de)<-paste(colnames(de),"MAST", sep="_")
# de$gene<-rownames(de)
# 
# plt<-merge(de, de_negbinom, by="gene")
# cor.test(plt$avg_log2FC_negbinom, plt$avg_log2FC_MAST)
# cor.test(plt$p_val_negbinom, plt$p_val_MAST, method="spearman")
# 
# ceiling(range(c(-log10(plt$p_val_MAST),-log10(plt$p_val_negbinom )), finite=TRUE)[2])
# 
# MAST_negbinom<-ggplot(plt, aes(-log10(p_val_MAST),-log10(p_val_negbinom )))+geom_point(shape=21, fill="grey",color="black")+ylim(0,500)+xlim(0,500)+geom_abline()+theme_bw()
# save_plts(
#   plot_grid(ggplot(plt, aes(-log10(p_val_MAST)))+geom_density(fill="grey")+xlim(0,500)+theme_bw()+xlab(""),plot.new(),
#           MAST_negbinom,ggplot(plt, aes(-log10(p_val_negbinom)))+geom_density(fill="grey")+xlim(0,500)+theme_bw()+xlab("")+coord_flip(),
#           align="v", rel_widths = c(1,0.25), rel_heights = c(0.25,1)),
#   "MAST_negbinom_bothlatent", w=6,h=5)
# 
# 
# plt_playing<-plt
# plt_playing$p_val_MAST[which(plt_playing$p_val_MAST==0)]<-5e-322
# plt_playing$p_val_negbinom[which(plt_playing$p_val_negbinom==0)]<-5e-322
# MAST_negbinom<-ggplot(plt_playing, aes(-log10(p_val_adj_MAST),-log10(p_val_adj_negbinom )))+geom_point(shape=21, fill="grey",color="black")+geom_abline()+theme_bw()
# plot_grid(ggplot(plt_playing, aes(-log10(p_val_adj_MAST)))+geom_histogram(fill="grey")+theme_bw()+xlab(""),plot.new(),
#             MAST_negbinom,ggplot(plt_playing, aes(-log10(p_val_adj_negbinom)))+geom_density(fill="grey")+theme_bw()+xlab("")+coord_flip(),
#             align="v", rel_widths = c(1,0.25), rel_heights = c(0.25,1))
# 
# 
# sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]
# sig_de_negbinom<-de_negbinom[which(de_negbinom$p_val_adj < 0.005 & abs(de_negbinom$avg_log2FC) > 1),]
# nrow(sig_de)
# nrow(sig_de_negbinom)
# length(intersect(sig_de$gene, sig_de_negbinom$gene))
# 
# 
# colnames(de_nolatent)<-paste(colnames(de_nolatent),"MAST_nolatent", sep="_")
# de_nolatent$gene<-rownames(de_nolatent)
# colnames(de)<-paste(colnames(de),"MAST", sep="_")
# de$gene<-rownames(de)
# 
# plt<-merge(de_nolatent, de, by="gene")
# cor.test(plt$avg_log2FC_MAST, plt$avg_log2FC_MAST_nolatent)
# identical(plt$avg_log2FC_MAST, plt$avg_log2FC_MAST_nolatent)
# cor.test(plt$p_val_MAST, plt$p_val_MAST_nolatent)
# 
# ceiling(range(c(-log10(plt$p_val_MAST),-log10(plt$p_val_MAST_nolatent )), finite=TRUE)[2])
# 
# ggplot(plt, aes(-log10(p_val_MAST),-log10(p_val_MAST_nolatent )))+geom_point()+ylim(0,325)+xlim(0,325)+geom_abline()
# ggplot(plt, aes(-log10(p_val_MAST_nolatent),-log10(p_val_MAST)))+geom_point()+ylim(0,325)+xlim(0,325)+geom_abline()
# 
# sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]
# sig_de_nolatent<-de_nolatent[which(de_nolatent$p_val_adj < 0.005 & abs(de_nolatent$avg_log2FC) > 1),]
# nrow(sig_de)
# nrow(sig_de_nolatent)
# length(intersect(sig_de$gene, sig_de_nolatent$gene))
# 
# plt$color<-"NS"
# plt$color[which(abs(plt$avg_log2FC_MAST) > 1)]<-"Sig FC"
# 
# MAST_latent<-ggplot(plt, aes(-log10(p_val_MAST),-log10(p_val_MAST_nolatent ), color=color))+geom_point()+ylim(0,500)+xlim(0,500)+geom_abline()+
#   geom_vline(xintercept=-log10(max(sig_de$p_val_MAST)))+geom_hline(yintercept=-log10(max(sig_de_nolatent$p_val_MAST_nolatent)))+
#   scale_color_manual(values = c("grey","blue"))+theme_bw()
# save_plts(MAST_latent, "MAST_latent", w=6,h=5)
# 
# 
# 
# colnames(de_negbinom)<-paste(colnames(de_negbinom),"negbinom", sep="_")
# de_negbinom$gene<-rownames(de_negbinom)
# colnames(de_negbinom_nolatent)<-paste(colnames(de_negbinom_nolatent),"nolatent", sep="_")
# de_negbinom_nolatent$gene<-rownames(de_negbinom_nolatent)
# 
# plt<-merge(de_negbinom_nolatent, de_negbinom, by="gene")
# cor.test(plt$avg_log2FC_negbinom, plt$avg_log2FC_nolatent)
# cor.test(plt$p_val_negbinom, plt$p_val_nolatent)
# 
# ceiling(range(c(-log10(plt$p_val_nolatent),-log10(plt$p_val_negbinom )), finite=TRUE)[2])
# 
# ggplot(plt, aes(-log10(p_val_nolatent),-log10(p_val_negbinom )))+geom_point()+ylim(0,325)+xlim(0,325)+geom_abline()
# ggplot(plt, aes(-log10(p_val_adj_nolatent),-log10(p_val_adj_negbinom)))+geom_point()+ylim(0,325)+xlim(0,325)+geom_abline()
# 
# 
# 
# load("../../Downloads/de_tmp.Rdata")
# dim(de_negbinom)
# dim(de_negbinom_nolatent)
# 
# de_negbinom<-de_negbinom[match(rownames(de_negbinom_nolatent),rownames(de_negbinom)),]
# identical(rownames(de_negbinom_nolatent), rownames(de_negbinom))
# cor.test(de_negbinom_nolatent$p_val, de_negbinom$p_val)
# 
# sig_de_negbinom<-de_negbinom[which(de_negbinom$p_val_adj < 0.005 & abs(de_negbinom$avg_log2FC) > 1),]
# sig_de_negbinom_nolatent<-de_negbinom_nolatent[which(de_negbinom_nolatent$p_val_adj < 0.005 & abs(de_negbinom_nolatent$avg_log2FC) > 1),]
# identical(rownames(sig_de_negbinom_nolatent), rownames(sig_de_negbinom))
# cor.test(sig_de_negbinom_nolatent$p_val, sig_de_negbinom$p_val)
# 
# 
# colnames(sig_de_negbinom)<-paste(colnames(sig_de_negbinom),"negbinom", sep="_")
# sig_de_negbinom$gene<-rownames(sig_de_negbinom)
# colnames(sig_de_negbinom_nolatent)<-paste(colnames(sig_de_negbinom_nolatent),"nolatent", sep="_")
# sig_de_negbinom_nolatent$gene<-rownames(sig_de_negbinom_nolatent)
# 
# plt<-merge(sig_de_negbinom_nolatent, sig_de_negbinom, by="gene")
# cor.test(plt$avg_log2FC_negbinom, plt$avg_log2FC_nolatent)
# cor.test(plt$p_val_negbinom, plt$p_val_nolatent)
# 
# ceiling(range(c(-log10(plt$p_val_nolatent),-log10(plt$p_val_negbinom )), finite=TRUE)[2])
# 
# ggplot(plt, aes(-log10(p_val_nolatent),-log10(p_val_negbinom )))+geom_point()+ylim(0,325)+xlim(0,325)+geom_abline()
# ggplot(plt, aes(-log10(p_val_adj_nolatent),-log10(p_val_adj_negbinom)))+geom_point()+ylim(0,325)+xlim(0,325)+geom_abline()
# 
# ggplot(plt, aes(-log10(p_val_negbinom)))+geom_density()
# ggplot(plt, aes(p_val_negbinom))+geom_density()
# ggplot(plt, aes(p_val_nolatent))+geom_histogram()
# 
# 


#################
## Sex age and between cell types
#################

########## DGE
Idents(d10x_raw_myeloid) <- "CellType_refined"
de_cell<-FindMarkers(d10x_raw_myeloid, ident.1 = "RR Myeloid", ident.2 = "KC Like", test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
de_cell_negbinom<-FindMarkers(d10x_raw_myeloid, ident.1 = "RR Myeloid", ident.2 = "KC Like", test.use = "negbinom",latent.vars="nFeature_RNA", verbose=F)


d10x_raw_RR<-subset(d10x, subset = CellType_refined %in% c("RR Myeloid"))

Idents(d10x_raw_RR) <- "AgeGroup"
de_age<-FindMarkers(d10x_raw_RR, ident.1 = "Adult", ident.2 = "Ped", test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
de_age_negbinom<-FindMarkers(d10x_raw_RR, ident.1 = "Adult", ident.2 = "Ped", test.use = "negbinom",latent.vars="nFeature_RNA", verbose=F)

Idents(d10x_raw_RR) <- "Sex"
de_sex<-FindMarkers(d10x_raw_RR, ident.1 = "M", ident.2 = "F", test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
de_sex_negbinom<-FindMarkers(d10x_raw_RR, ident.1 = "M", ident.2 = "F", test.use = "negbinom",latent.vars="nFeature_RNA", verbose=F)


save(de_cell,de_cell_negbinom, de_age, de_age_negbinom, de_sex,de_sex_negbinom, file=here("data/DGE_compare.Rdata"))



load(here("data/DGE_compare.Rdata"))

## grab legened from plot
get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

# sig_de<-de_sex[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]
# sig_de_negbinom<-de_sex_negbinom[which(de_negbinom$p_val_adj < 0.005 & abs(de_negbinom$avg_log2FC) > 1),]



plto_DGE<-function(MAST, negbinom){
  colnames(MAST)<-paste(colnames(MAST),"MAST", sep="_")
  MAST$gene<-rownames(MAST)
  colnames(negbinom)<-paste(colnames(negbinom),"negbinom", sep="_")
  negbinom$gene<-rownames(negbinom)
  
  plt<-merge(MAST, negbinom, by="gene")
  cor.test(plt$avg_log2FC_MAST, plt$avg_log2FC_negbinom, method="spearman")
  cor.test(plt$p_val_MAST, plt$p_val_negbinom, method="spearman")
  
  lim<-ceiling(range(c(-log10(plt$p_val_MAST),-log10(plt$p_val_negbinom )), finite=TRUE)[2])
  
  sig_de<-MAST[which(MAST$p_val_adj < 0.005 & abs(MAST$avg_log2FC) > 1),]
  sig_de_negbiom<-negbinom[which(negbinom$p_val_adj < 0.005 & abs(negbinom$avg_log2FC) > 1),]
  
  both<-paste("Sig. Both (",length(intersect(sig_de$gene, sig_de_negbiom$gene)), ")", sep="")
  sigmast<-paste("Sig. MAST (",length(which(!(sig_de$gene%in%sig_de_negbiom$gene) & sig_de$gene%in%plt$gene)), ")", sep="")
  signeg<-paste("Sig. Negbinom (",length(which(!(sig_de_negbiom$gene%in%sig_de$gene))), ")", sep="")
  ns<-paste("NS (",length(which(!(plt$gene%in%c(sig_de$gene, sig_de_negbiom$gene)))), ")", sep="")
  
  
  plt$color<-sapply(1:nrow(plt), function(x){
    if(plt$gene[x]%in%sig_de$gene & plt$gene[x]%in%sig_de_negbiom$gene){both}else{
      if(plt$gene[x]%in%sig_de$gene){sigmast}else{
        if(plt$gene[x]%in%sig_de_negbiom$gene){signeg}else{ns}}
    }
  })
  
  plt<-plt[order(plt$color),]
  MAST_negbinom<-ggplot(plt, aes(-log10(p_val_MAST),-log10(p_val_negbinom), color=color))+geom_point(shape=19)+ylim(0,lim)+xlim(0,lim)+geom_abline()+
    geom_vline(xintercept=-log10(max(sig_de$p_val_MAST)), color="grey")+geom_hline(yintercept=-log10(max(sig_de_negbiom$p_val_negbinom)), color="grey")+
    scale_color_manual(values = c("lightgrey","#1b7837","#6f99ff","#ffc026"), name="Significant\n(FDR <0.005\n|FC| > 1)")+theme_bw()+
    annotate("text", x=250, y=1, label = paste0("Rs = ",signif(cor.test(plt$p_val_MAST, plt$p_val_negbinom, method="spearman")$estimate, 2)))+
    xlab("MAST\nP Value (-log10)")+ylab("Negative Binomial\nP Value (-log10)")
  
  plot_grid(ggplot(plt, aes(-log10(p_val_MAST)))+geom_density(fill="grey")+theme_bw()+xlab(""),get_leg(MAST_negbinom),
              MAST_negbinom +theme(legend.position = "none") ,
            ggplot(plt, aes(-log10(p_val_negbinom)))+geom_density(fill="grey")+theme_bw()+xlab("")+coord_flip(),
              align="v", axis = "l",rel_widths = c(1,0.25), rel_heights = c(0.25,1))}


plto_DGE(de_cell, de_cell_negbinom)
save_plts(plto_DGE(de_cell, de_cell_negbinom), "MAST_negbinom_celltype", w=8,h=8)

plto_DGE(de_sex, de_sex_negbinom)
save_plts(plto_DGE(de_sex, de_sex_negbinom), "MAST_negbinom_sex", w=8,h=8)

plto_DGE(de_age, de_age_negbinom)
save_plts(plto_DGE(de_age, de_age_negbinom), "MAST_negbinom_age", w=8,h=8)



