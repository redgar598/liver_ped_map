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
source("scripts/00_Heat_scree_plot_generic.R")


load(here("/media/redgar/Seagate Portable Drive/processed_data/Fetal_IFALD_adult_ped_integrated_myeloid_only.RData"))

d10x.combined_myeloid[["pca"]]

DimPlot(d10x.combined_myeloid[order(d10x.combined_myeloid$cell),], reduction="umap", group.by="age_condition")

head(Embeddings(d10x.combined_myeloid, reduction = "pca")[, 1:5])
head(Loadings(d10x.combined_myeloid, reduction = "pca")[, 1:5])
head(Stdev(d10x.combined_myeloid, reduction = "pca"))

#' ## PCA for batch effect
Loadings<-as.data.frame(Loadings(d10x.combined_myeloid, reduction = "pca"))
embed<-as.data.frame(Embeddings(d10x.combined_myeloid, reduction = "pca"))
vars <- Stdev(d10x.combined_myeloid, reduction = "pca")^2
Importance<-vars/sum(vars)
print(Importance[1:10])

meta_categorical <- d10x.combined_myeloid@meta.data[, c("CellType_refined","CellType_harmonized","age_condition")]  # input column numbers in meta that contain categorical variables
meta_continuous <- d10x.combined_myeloid@meta.data[, c("percent.mt","nFeature_RNA","Age")]  # input column numbers in meta that contain continuous variables



meta_continuous$Age<-as.factor(meta_continuous$Age)
# 40 weeks gestation I am calling term so -(40-x)*(1/52)
levels(meta_continuous$Age)<-c("11","-0.56", "12","-0.54",
                               "13","-0.52", "-0.5", 
                               "-0.46","17","-0.44", 
                               "2","26" , "3","48","57","61",
                               "65","67","69",
                               "-0.64", "-0.62" , "-0.60" )
meta_continuous$Age<-as.numeric(as.character(meta_continuous$Age))

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-30

suppressWarnings(heat_scree_plot(embed, Importance, 2, 2))


meta_df<-cbind(meta_continuous,meta_categorical, embed)


cor.test(meta_continuous$Age, embed$PC_2)

ggplot(meta_df, aes(PC_1,PC_2, color=Age))+geom_point()
ggplot(meta_df, aes(PC_1,PC_2, color=age_condition))+geom_point()+colscale_agecondition_fetal
ggplot(meta_df, aes(PC_1,PC_2, color=CellType_harmonized))+geom_point()+colscale_cellType_fetal_combo
ggplot(meta_df, aes(PC_1,PC_2, color=Age))+geom_point()+facet_wrap(~age_condition)

ggplot(meta_df, aes(PC_2,PC_3, color=Age))+geom_point()
ggplot(meta_df[order(meta_df$Age),], aes(PC_4,PC_5, color=age_condition))+geom_point()+colscale_agecondition_fetal





ggplot(meta_df, aes(CellType_refined,PC_2))+geom_boxplot()
ggplot(meta_df, aes(age_condition,PC_1))+geom_boxplot()
ggplot(meta_df, aes(age_condition,PC_7))+geom_boxplot()
ggplot(meta_df, aes(age_condition,PC_2))+geom_boxplot()


ggplot(meta_df, aes(PC_1, PC_2, color=CellType_refined))+geom_point()+scale_color_manual(values=combo_colors)
ggplot(meta_df, aes(PC_1, PC_2, color=CellType_harmonized))+geom_point()+scale_color_manual(values=combo_colors)


ggplot(meta_df, aes(PC_1, PC_2, fill=age_condition))+
  geom_point(shape=21)+
  fillscale_agecondition_fetal+
  theme_bw()

ggplot(meta_df, aes(PC_1, PC_2, fill=age_condition))+
  geom_point(shape=21)+
  scale_fill_manual(values=c("cornflowerblue","darkgreen","firebrick3","tomato2"))+
  theme_bw()+facet_wrap(~CellType_harmonized)



rownames(Loadings)[order(Loadings$PC_2)][1:20]
rownames(Loadings)[rev(order(Loadings$PC_2))][1:20]

rownames(Loadings)[order(Loadings$PC_3)][1:20]
rownames(Loadings)[rev(order(Loadings$PC_3))][1:20]

