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


source("R_functions/pretty_plots.R")


load(here("data","adult_ped_integrated.rds"))

d10x.combined_myeloid<-subset(d10x.combined, subset = CellType_rough %in% c("Myeloid"))
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
d10x_raw_myeloid<-subset(d10x, cells = rownames(myeliod_clusters))
rm(d10x)
gc()

myeliod_clusters$index<-rownames(myeliod_clusters)
identical(colnames(d10x_raw_myeloid), myeliod_clusters$index)

d10x_raw_myeloid <- AddMetaData(d10x_raw_myeloid, metadata = myeliod_clusters)

d10x_raw_myeloid <- NormalizeData(d10x_raw_myeloid,scale.factor = 10000, normalization.method = "LogNormalize")


d10x_raw_myeloid@meta.data$CellType_rough[which(d10x_raw_myeloid@meta.data$seurat_clusters=="1")]<-"RR_myeloid"
d10x_raw_myeloid@meta.data$CellType_rough[which(d10x_raw_myeloid@meta.data$seurat_clusters%in%c("3","4"))]<-"KC_like"
d10x_raw_myeloid@meta.data$CellType_rough[which(d10x_raw_myeloid@meta.data$seurat_clusters=="7")]<-"macro_RBC"
d10x_raw_myeloid@meta.data$CellType_rough[which(d10x_raw_myeloid@meta.data$seurat_clusters=="5")]<-"bcell"
d10x_raw_myeloid@meta.data$CellType_rough[which(d10x_raw_myeloid@meta.data$seurat_clusters%in%c("2","6"))]<-"Neutrophil"
d10x_raw_myeloid@meta.data$CellType_rough[which(d10x_raw_myeloid@meta.data$seurat_clusters=="0")]<-"badmyeloid"


########## DGE
Idents(d10x_raw_myeloid) <- "CellType_rough"

de<-FindMarkers(d10x_raw_myeloid, ident.1 = "RR_myeloid", ident.2 = "KC_like", test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]

de_nolatent<-FindMarkers(d10x_raw_myeloid, ident.1 = "RR_myeloid", ident.2 = "KC_like", test.use = "MAST", verbose=F)
sig_de_nolatent<-de_nolatent[which(de_nolatent$p_val_adj < 0.005 & abs(de_nolatent$avg_log2FC) > 1),]

de_negbinom<-FindMarkers(d10x_raw_myeloid, ident.1 = "RR_myeloid", ident.2 = "KC_like", test.use = "negbinom",latent.vars="nFeature_RNA", verbose=F)
sig_de_negbinom<-de_negbinom[which(de_negbinom$p_val_adj < 0.005 & abs(de_negbinom$avg_log2FC) > 1),]

de_negbinom_match<-de_negbinom[which(rownames(de_negbinom)%in%rownames(de)),]
de_negbinom_match<-de_negbinom_match[match(rownames(de),rownames(de_negbinom_match)),]
identical(rownames(de), rownames(de_negbinom_match))


de_negbinom_nolatent<-FindMarkers(d10x_raw_myeloid, ident.1 = "RR_myeloid", ident.2 = "KC_like", test.use = "negbinom", verbose=F)
sig_de_negbinom_nolatent<-de_negbinom_nolatent[which(de_negbinom_nolatent$p_val_adj < 0.005 & abs(de_negbinom_nolatent$avg_log2FC) > 1),]

de_negbinom<-de_negbinom[match(rownames(de_negbinom_nolatent),rownames(de_negbinom)),]
identical(rownames(de_negbinom_nolatent), rownames(de_negbinom))
identical(de_negbinom_nolatent$p_val, de_negbinom$p_val)
cor(de_negbinom_nolatent$p_val, de_negbinom$p_val)

# de_DESeq2<-FindMarkers(d10x_raw_myeloid, ident.1 = "RR_myeloid", ident.2 = "KC_like", test.use = "DESeq2",latent.vars="nFeature_RNA", verbose=F)
# sig_de_de_DESeq2<-de_DESeq2[which(de_DESeq2$p_val_adj < 0.005 & abs(de_DESeq2$avg_log2FC) > 1),]
# 
# de_DESeq2_nolatent<-FindMarkers(d10x_raw_myeloid, ident.1 = "RR_myeloid", ident.2 = "KC_like", test.use = "DESeq2", verbose=F)
# sig_de_de_DESeq2_nolatent<-de_DESeq2_nolatent[which(de_DESeq2_nolatent$p_val_adj < 0.005 & abs(de_DESeq2_nolatent$avg_log2FC) > 1),]

save(de,de_nolatent, de_negbinom, de_negbinom_nolatent, file="../../Downloads/de_tmp.Rdata")

load("../../Downloads/de_tmp.Rdata")

colnames(de_negbinom)<-paste(colnames(de_negbinom),"negbinom", sep="_")
de_negbinom$gene<-rownames(de_negbinom)
colnames(de)<-paste(colnames(de),"MAST", sep="_")
de$gene<-rownames(de)

plt<-merge(de, de_negbinom, by="gene")
cor.test(plt$avg_log2FC_negbinom, plt$avg_log2FC_MAST)
cor.test(plt$p_val_negbinom, plt$p_val_MAST)

ceiling(range(c(-log10(plt$p_val_MAST),-log10(plt$p_val_negbinom )), finite=TRUE)[2])

MAST_negbinom<-ggplot(plt, aes(-log10(p_val_MAST),-log10(p_val_negbinom )))+geom_point()+ylim(0,500)+xlim(0,500)+geom_abline()
save_plts(MAST_negbinom, "MAST_negbinom_bothlatent", w=6,h=5)

ggplot(plt, aes(-log10(p_val_negbinom),-log10(p_val_MAST)))+geom_point()+ylim(0,325)+xlim(0,325)+geom_abline()

sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]
sig_de_negbinom<-de_negbinom[which(de_negbinom$p_val_adj < 0.005 & abs(de_negbinom$avg_log2FC) > 1),]
nrow(sig_de)
nrow(sig_de_negbinom)
length(intersect(sig_de$gene, sig_de_negbinom$gene))


colnames(de_nolatent)<-paste(colnames(de_nolatent),"MAST_nolatent", sep="_")
de_nolatent$gene<-rownames(de_nolatent)
colnames(de)<-paste(colnames(de),"MAST", sep="_")
de$gene<-rownames(de)

plt<-merge(de_nolatent, de, by="gene")
cor.test(plt$avg_log2FC_MAST, plt$avg_log2FC_MAST_nolatent)
identical(plt$avg_log2FC_MAST, plt$avg_log2FC_MAST_nolatent)
cor.test(plt$p_val_MAST, plt$p_val_MAST_nolatent)

ceiling(range(c(-log10(plt$p_val_MAST),-log10(plt$p_val_MAST_nolatent )), finite=TRUE)[2])

ggplot(plt, aes(-log10(p_val_MAST),-log10(p_val_MAST_nolatent )))+geom_point()+ylim(0,325)+xlim(0,325)+geom_abline()
ggplot(plt, aes(-log10(p_val_MAST_nolatent),-log10(p_val_MAST)))+geom_point()+ylim(0,325)+xlim(0,325)+geom_abline()

sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]
sig_de_nolatent<-de_nolatent[which(de_nolatent$p_val_adj < 0.005 & abs(de_nolatent$avg_log2FC) > 1),]
nrow(sig_de)
nrow(sig_de_nolatent)
length(intersect(sig_de$gene, sig_de_nolatent$gene))

plt$color<-"NS"
plt$color[which(abs(plt$avg_log2FC_MAST) > 1)]<-"Sig FC"

MAST_latent<-ggplot(plt, aes(-log10(p_val_MAST),-log10(p_val_MAST_nolatent ), color=color))+geom_point()+ylim(0,500)+xlim(0,500)+geom_abline()+
  geom_vline(xintercept=-log10(max(sig_de$p_val_MAST)))+geom_hline(yintercept=-log10(max(sig_de_nolatent$p_val_MAST_nolatent)))+
  scale_color_manual(values = c("grey","blue"))+theme_bw()
save_plts(MAST_latent, "MAST_latent", w=6,h=5)



colnames(de_negbinom)<-paste(colnames(de_negbinom),"negbinom", sep="_")
de_negbinom$gene<-rownames(de_negbinom)
colnames(de_negbinom_nolatent)<-paste(colnames(de_negbinom_nolatent),"nolatent", sep="_")
de_negbinom_nolatent$gene<-rownames(de_negbinom_nolatent)

plt<-merge(de_negbinom_nolatent, de_negbinom, by="gene")
cor.test(plt$avg_log2FC_negbinom, plt$avg_log2FC_nolatent)
cor.test(plt$p_val_negbinom, plt$p_val_nolatent)

ceiling(range(c(-log10(plt$p_val_nolatent),-log10(plt$p_val_negbinom )), finite=TRUE)[2])

ggplot(plt, aes(-log10(p_val_nolatent),-log10(p_val_negbinom )))+geom_point()+ylim(0,325)+xlim(0,325)+geom_abline()
ggplot(plt, aes(-log10(p_val_adj_nolatent),-log10(p_val_adj_negbinom)))+geom_point()+ylim(0,325)+xlim(0,325)+geom_abline()



load("../../Downloads/de_tmp.Rdata")
dim(de_negbinom)
dim(de_negbinom_nolatent)

de_negbinom<-de_negbinom[match(rownames(de_negbinom_nolatent),rownames(de_negbinom)),]
identical(rownames(de_negbinom_nolatent), rownames(de_negbinom))
cor.test(de_negbinom_nolatent$p_val, de_negbinom$p_val)

sig_de_negbinom<-de_negbinom[which(de_negbinom$p_val_adj < 0.005 & abs(de_negbinom$avg_log2FC) > 1),]
sig_de_negbinom_nolatent<-de_negbinom_nolatent[which(de_negbinom_nolatent$p_val_adj < 0.005 & abs(de_negbinom_nolatent$avg_log2FC) > 1),]
identical(rownames(sig_de_negbinom_nolatent), rownames(sig_de_negbinom))
cor.test(sig_de_negbinom_nolatent$p_val, sig_de_negbinom$p_val)


colnames(sig_de_negbinom)<-paste(colnames(sig_de_negbinom),"negbinom", sep="_")
sig_de_negbinom$gene<-rownames(sig_de_negbinom)
colnames(sig_de_negbinom_nolatent)<-paste(colnames(sig_de_negbinom_nolatent),"nolatent", sep="_")
sig_de_negbinom_nolatent$gene<-rownames(sig_de_negbinom_nolatent)

plt<-merge(sig_de_negbinom_nolatent, sig_de_negbinom, by="gene")
cor.test(plt$avg_log2FC_negbinom, plt$avg_log2FC_nolatent)
cor.test(plt$p_val_negbinom, plt$p_val_nolatent)

ceiling(range(c(-log10(plt$p_val_nolatent),-log10(plt$p_val_negbinom )), finite=TRUE)[2])

ggplot(plt, aes(-log10(p_val_nolatent),-log10(p_val_negbinom )))+geom_point()+ylim(0,325)+xlim(0,325)+geom_abline()
ggplot(plt, aes(-log10(p_val_adj_nolatent),-log10(p_val_adj_negbinom)))+geom_point()+ylim(0,325)+xlim(0,325)+geom_abline()

ggplot(plt, aes(-log10(p_val_negbinom)))+geom_density()
ggplot(plt, aes(p_val_negbinom))+geom_density()
ggplot(plt, aes(p_val_nolatent))+geom_histogram()
