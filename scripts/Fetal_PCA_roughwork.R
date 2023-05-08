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


load(here("data","Fetal_ped_IFALD_adult_PCA_myeloid.RData"))

meta_continuous$Age<-as.factor(meta_continuous$Age)
# 40 weeks gestation I am calling term so -(40-x)*(1/52)
levels(meta_continuous$Age)<-c("11","-0.56", "12","-0.54",
                               "13","-0.52", "-0.5", 
                               "-0.46","17","-0.44", 
                               "2","26" , "3","48","57","61",
                               "65","67","69",
                               "-0.64", "-0.62" , "-0.60" )
meta_continuous$Age<-as.numeric(meta_continuous$Age)

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(embed, Importance, 2, 2))


cor.test(meta_continuous$Age, embed$PC_2)
scatter_df<-cbind(meta_continuous, embed)
ggplot(scatter_df, aes(PC_1,PC_2, color=Age))+geom_point()

boxplt_df<-cbind(meta_categorical, embed)
ggplot(boxplt_df, aes(CellType_refined,PC_1))+geom_boxplot()
ggplot(boxplt_df, aes(age_condition,PC_1))+geom_boxplot()

ggplot(boxplt_df, aes(PC_1, PC_2, color=CellType_refined))+geom_point()+scale_color_manual(values=combo_colors)
ggplot(boxplt_df, aes(PC_1, PC_2, color=CellType_refined))+geom_point()+scale_color_manual(values=combo_colors)


ggplot(boxplt_df, aes(PC_1, PC_2, fill=age_condition))+
  geom_point(shape=21)+
  scale_fill_manual(values=c("cornflowerblue","darkgreen","firebrick3","tomato2"))+
  theme_bw()

ggplot(boxplt_df, aes(PC_1, PC_2, fill=age_condition))+
  geom_point(shape=21)+
  scale_fill_manual(values=c("cornflowerblue","darkgreen","firebrick3","tomato2"))+
  theme_bw()+facet_wrap(~CellType_refined)


boxplt_df$age_condition<-factor(boxplt_df$age_condition, levels=c("Fetal Healthy","Ped Healthy","Ped IFALD","Adult Healthy" ))
plot_grid(ggplot(boxplt_df[which(boxplt_df$CellType_refined%in%c("RR Myeloid","Monocyte")),], aes(PC_1, PC_2, fill=age_condition))+
  geom_point(shape=21)+
  scale_fill_manual(values=c("darkgreen","firebrick3","tomato2","cornflowerblue"))+
  theme_bw(),
  ggplot(boxplt_df[which(boxplt_df$CellType_refined%in%c("RR Myeloid","Monocyte")),], aes(age_condition, PC_2))+
    geom_violin()+
    geom_boxplot(width=0.1, aes(fill=age_condition))+
    scale_fill_manual(values=c("darkgreen","firebrick3","tomato2","cornflowerblue"))+
    theme_bw(), rel_widths = c(2,1))


ggplot(boxplt_df[which(boxplt_df$CellType_refined%in%c("RR Myeloid","Monocyte")),], aes(PC_1, PC_2, fill=age_condition))+
  geom_point(shape=21)+
  scale_fill_manual(values=c("darkgreen","firebrick3","tomato2","cornflowerblue"))+
  theme_bw()+facet_wrap(~age_condition)
  

rownames(Loadings)[order(Loadings$PC_2)][1:20]
rownames(Loadings)[rev(order(Loadings$PC_2))][1:20]


