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

## On integrated data? And only top variable
counts<-GetAssayData(object = d10x.combined_myeloid, slot = "data")
counts<-as.matrix(counts)

counts_mini<-counts[1:100,]

## PCA followed by varimax rotation
count_PCA<-prcomp(t(counts_mini))
Loadings_counts<-count_PCA$x

factors<-varimax(Loadings_counts)

factors_loadings<-as.matrix(factors$loadings)


x <- loadings(factors)
vx <- colSums(x^2)

importance<-rbind(`SS loadings` = vx,
      `Proportion Var` = vx/nrow(x),
      `Cumulative Var` = cumsum(vx/nrow(x)))

importance[2,]

#extract numeric as a data frame
factors_loadings<-data.frame(matrix(as.numeric(factors$loadings), attributes(factors$loadings)$dim, dimnames=attributes(factors$loadings)$dimnames))


Loadings<-factors_loadings
Importance<-importance[2,]
PCs_to_view<-10
ord<-1:10

meta<-d10x.combined_myeloid@meta.data
meta_mini<-meta[which(meta$cell%in%rownames(factors_loadings)),]
meta_categorical<-meta_mini[,c("seurat_clusters","AgeGroup","individual")] # include in full ,"Chemistry","FreshorFrozen"
meta_continuous<-meta_mini[,c("percent.mt","Age","nFeature_RNA","nCount_RNA")]

heat_scree_plot(Loadings,Importance, 0.25, 0.25 )


factors_loadings$cell<-rownames(factors_loadings)
meta$cell<-rownames(meta)

Loadings_meta<-merge(factors_loadings, meta, by="cell")
ggplot(Loadings_meta, aes(PC1, PC2, fill=seurat_clusters, color=AgeGroup))+
  geom_point(shape=21,size=2)+theme_bw()+scale_color_manual(values=c("black","white"))

cor.test(Loadings_meta$PC2, Loadings_meta$percent.mt, method = "spearman", na.action = na.omit, exact = FALSE )
ggplot(Loadings_meta, aes(percent.mt, PC2))+
  geom_point(aes(fill=seurat_clusters, color=AgeGroup),shape=21,size=2)+theme_bw()+scale_color_manual(values=c("black","white"))+
  stat_smooth(method="lm")

summary(aov(Loadings_meta$PC3 ~ Loadings_meta$AgeGroup))
ggplot(Loadings_meta, aes(AgeGroup, PC3))+
  geom_point(aes(fill=seurat_clusters, color=AgeGroup),shape=21,size=2)+geom_boxplot()+theme_bw()+scale_color_manual(values=c("black","white"))+
  stat_smooth(method="lm")



Loadings_counts<-as.data.frame(Loadings_counts)
Loadings_counts$cell<-rownames(Loadings_counts)
Loadings_meta_PCA<-merge(Loadings_counts, meta, by="cell")

ggplot(Loadings_meta_PCA, aes(PC1, PC2, fill=seurat_clusters, color=AgeGroup))+
  geom_point(shape=21,size=2)+theme_bw()+scale_color_manual(values=c("black","white"))

################### heat scree function

myColors_pval <- c("#084594","#4292c6","#9ecae1","#deebf7")
color_possibilities_pval<-c( "<=0.001","<=0.01","<=0.05",">0.05")
names(myColors_pval) <- color_possibilities_pval
fillscale_pval <- scale_fill_manual(name="P Value",values = myColors_pval, drop = F)


heat_scree_plot<-function(Loadings, Importance, right_marg, left_marg){
  
  if(missing(right_marg)) {
    right_marg=2.25} 
  
  if(missing(left_marg)) {
    left_marg=1} 
  
  pca_df<-data.frame(variance=Importance, PC=seq(1:length(Importance)))
  
  scree<-ggplot(pca_df[which(pca_df$PC<=(PCs_to_view)),],aes(PC,variance))+
    geom_bar(stat = "identity",color="black",fill="grey")+theme_bw()+
    theme(axis.text.y = element_text(size =10, color="black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size =12),
          plot.margin=unit(c(1.25,2.8,-0.2,2.8),"cm"))+ylab("Variance")+
    scale_x_continuous(breaks = seq(1,PCs_to_view,1))+xlab("")
  
  
  #### Heat
  ## correlate meta with PCS
  ## Run anova of each PC on each meta data variable
  
  
  aov_PC_meta <- lapply(1:ncol(meta_categorical), function(covar) sapply(1:ncol(Loadings), 
                                                                         function(PC) 
                                                                           summary(aov(Loadings[, PC] ~ meta_categorical[, covar]))[[1]]$"Pr(>F)"[1]))
  cor_PC_meta <- lapply(1:ncol(meta_continuous), function(covar) sapply(1:ncol(Loadings), 
                                                                        function(PC) (cor.test(Loadings[, PC], as.numeric(meta_continuous[, 
                                                                                                                                          covar]), alternative = "two.sided", method = "spearman", na.action = na.omit)$p.value)))
  names(aov_PC_meta) <- colnames(meta_categorical)
  names(cor_PC_meta) <- colnames(meta_continuous)
  aov_PC_meta <- do.call(rbind, aov_PC_meta)
  cor_PC_meta <- do.call(rbind, cor_PC_meta)
  aov_PC_meta <- rbind(aov_PC_meta, cor_PC_meta)
  aov_PC_meta <- as.data.frame(aov_PC_meta)
  
  #adjust
  aov_PC_meta_adjust<-aov_PC_meta[,1:ncol(aov_PC_meta)]
  
  #reshape
  avo<-aov_PC_meta_adjust[,1:PCs_to_view]
  avo_heat_num<-apply(avo,2, as.numeric)
  avo_heat<-as.data.frame(avo_heat_num)
  avo_heat$meta<-rownames(avo)
  avo_heat_melt<-melt(avo_heat, id=c("meta"))
  
  # cluster meta data
  meta_var_order<-unique(avo_heat_melt$meta)[rev(ord)]
  avo_heat_melt$meta <- factor(avo_heat_melt$meta, levels = meta_var_order)
  
  # color if sig
  avo_heat_melt$Pvalue<-sapply(1:nrow(avo_heat_melt), function(x) if(avo_heat_melt$value[x]<=0.001){"<=0.001"}else{
    if(avo_heat_melt$value[x]<=0.01){"<=0.01"}else{
      if(avo_heat_melt$value[x]<=0.05){"<=0.05"}else{">0.05"}}})
  
  levels(avo_heat_melt$variable)<-sapply(1:PCs_to_view, function(x) paste("PC",x, sep="" ))
  
  myColors_pval <- c("#084594","#4292c6","#9ecae1","#deebf7")
  color_possibilities_pval<-c( "<=0.001","<=0.01","<=0.05",">0.05")
  names(myColors_pval) <- color_possibilities_pval
  fillscale_pval <- scale_fill_manual(name="P Value",values = myColors_pval, drop = F)
  
  heat<-ggplot(avo_heat_melt, aes(variable,meta, fill = Pvalue)) +
    geom_tile(color = "black",size=0.5) +
    theme_gray(8)+fillscale_pval+
    theme(axis.text = element_text(size =9, color="black"),
          axis.title = element_text(size =11),
          legend.text = element_text(size =10),
          legend.title = element_text(size =11),
          legend.position = c(1.35, 0.75), legend.justification = c(1,1),
          plot.margin=unit(c(-0.3,right_marg,1,left_marg),"cm"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    xlab("Principal Component")+ylab(NULL)
  
  grid.arrange(scree, heat, ncol=1,heights = c(3, 4))#
}