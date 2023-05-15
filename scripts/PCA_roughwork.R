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
library(ggsignif)


source("scripts/00_pretty_plots.R")
source("scripts/00_Heat_scree_plot_generic.R")


d10x_ped_IFALD<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))
## add cell type labels
load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x_ped_IFALD), cell_label$index),]
identical(colnames(d10x_ped_IFALD), cell_label$index)
d10x_ped_IFALD <- AddMetaData(d10x_ped_IFALD, metadata = cell_label)


d10x.combined_myeloid<-subset(d10x_ped_IFALD, subset = CellType_refined %in% c("KC Like","RR Myeloid","Macrophage\n(CLEC9A high)","Macrophage\n(MHCII high)"))
rm(d10x_ped_IFALD)
gc

d10x.combined_myeloid <- NormalizeData(d10x.combined_myeloid)
d10x.combined_myeloid <- FindVariableFeatures(d10x.combined_myeloid, selection.method = "vst", nfeatures = 2000)
d10x.combined_myeloid <- ScaleData(d10x.combined_myeloid) 
d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)


d10x.combined_myeloid[["pca"]]

head(Embeddings(d10x.combined_myeloid, reduction = "pca")[, 1:5])
head(Loadings(d10x.combined_myeloid, reduction = "pca")[, 1:5])
head(Stdev(d10x.combined_myeloid, reduction = "pca"))



#' ## PCA for batch effect
    #Loadings<-as.data.frame(Loadings(d10x.combined_myeloid, reduction = "pca"))
embed<-as.data.frame(Embeddings(d10x.combined_myeloid, reduction = "pca"))
vars <- Stdev(d10x.combined_myeloid, reduction = "pca")^2
Importance<-vars/sum(vars)
print(Importance[1:10])


meta_categorical <- d10x.combined_myeloid@meta.data[, c("CellType_refined","age_condition","Phase")]  # input column numbers in meta that contain categorical variables
meta_continuous <- d10x.combined_myeloid@meta.data[, c("percent.mt","nFeature_RNA","Age")]  # input column numbers in meta that contain continuous variables
        # colnames(meta_categorical) <- c("Individual","Diagnosis","Segment","Gender","Sentrix ID","WNT Type","Batch")
        # colnames(meta_continuous) <- c( "Age","Passage")

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(embed, Importance, 2, 2))

cor.test(meta_continuous$Age, embed$PC_1)
scatter_df<-cbind(meta_continuous, embed)
ggplot(scatter_df, aes(PC_1,PC_2, color=Age))+geom_point()

aov(embed[, 1] ~ meta_categorical$age_condition)$coef
aov(embed[, 1] ~ meta_categorical$CellType_refined)$coef
aov(embed[, 1] ~ meta_categorical$Phase)$coef

summary(aov(embed[, 1] ~ meta_categorical$age_condition))

boxplt_df<-cbind(meta_categorical, embed)
ggplot(boxplt_df, aes(CellType_refined,PC_1))+geom_boxplot()
ggplot(boxplt_df, aes(age_condition,PC_1))+geom_boxplot()

ggplot(boxplt_df, aes(PC_1, PC_2, color=age_condition))+
  geom_point()+
  colscale_agecondition+theme_bw()

ggplot(boxplt_df, aes(PC_1, PC_2, color=age_condition))+
  geom_point(size=0.75)+
  colscale_agecondition+theme_bw()+facet_wrap(~CellType_refined)


ggplot(boxplt_df, aes(PC_1, PC_2, color=CellType_refined))+geom_point()+colscale_cellType+theme_bw()


comp_simple<-list(c("Adult Healthy","Ped Healthy"),c("Adult Healthy","Ped IFALD"),c("Ped Healthy","Ped IFALD"))
ggplot(boxplt_df, aes(age_condition,PC_2))+
  geom_violin(fill="grey80", color="white")+
  geom_boxplot(aes(fill=age_condition), width=0.1)+
  facet_wrap(~CellType_refined)+  
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")+
  fillscale_agecondition+theme_bw()

ggplot(boxplt_df, aes(age_condition,PC_2))+
  geom_violin(fill="grey80", color="white")+
  geom_boxplot(aes(fill=age_condition), width=0.1)+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")+
  fillscale_agecondition+theme_bw()

###########
## gene in PC2
###########
Loadings<-as.data.frame(Loadings(d10x.combined_myeloid, reduction = "pca"))
rownames(Loadings)[order(Loadings$PC_2)][1:10]
rownames(Loadings)[rev(order(Loadings$PC_2))][1:20]

which(rownames(Loadings)[rev(order(Loadings$PC_2))]=="CCL4")
