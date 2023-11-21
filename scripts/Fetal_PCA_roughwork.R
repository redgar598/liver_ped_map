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
library(ggridges)

source("scripts/00_pretty_plots.R")
source("scripts/00_Heat_scree_plot_generic.R")
source("scripts/varimax_delaram.R")


#load(here("/media/redgar/Seagate Portable Drive/processed_data/Fetal_IFALD_adult_ped_integrated_myeloid_only.RData"))


                            

## only KC
# d10x.combined_KC<-subset(d10x.combined_myeloid, subset = CellType_harmonized %in% c("KC Like"))
# save(d10x.combined_KC, file=here("/media/redgar/Seagate Portable Drive/processed_data/Fetal_IFALD_adult_ped_integrated_KC_only.RData"))

load(here("/media/redgar/Seagate Portable Drive/processed_data/Fetal_IFALD_adult_ped_integrated_KC_only.RData"))

d10x.combined_KC <- RunPCA(d10x.combined_KC, npcs = 30, verbose = FALSE)
d10x.combined_KC <- RunUMAP(d10x.combined_KC, reduction = "pca", dims = 1:30)

loadings<-as.data.frame(Loadings(d10x.combined_KC, reduction = "pca"))


#####################
## KC PCA
#####################
d10x.combined_KC

DimPlot(d10x.combined_KC, reduction="pca")

d10x.combined_KC[["pca"]]

head(Embeddings(d10x.combined_KC, reduction = "pca")[, 1:5])
head(Loadings(d10x.combined_KC, reduction = "pca")[, 1:5])
head(Stdev(d10x.combined_KC, reduction = "pca"))

head(loadings[order(loadings$PC_1),1:5])
tail(loadings[order(loadings$PC_1),1:5])
FeaturePlot(d10x.combined_KC, reduction="pca", features = "MARCO")
FeaturePlot(d10x.combined_KC, reduction="pca", features = "ALB")

lapply(1:10, function(x) rownames(loadings)[order(loadings[,x])][1:5])
lapply(1:10, function(x) rownames(loadings)[rev(order(loadings[,x]))][1:5])


FeaturePlot(d10x.combined_KC, reduction="pca", features =c("IL1B","CXCL3","CCL4L2","CCL4"), dims=c(1,2))
FeaturePlot(d10x.combined_KC, reduction="pca", features =c("HLA-DQA1","HLA-DQA2","HLA-DPB1","CLEC1B"), dims=c(2,3))

DimPlot(d10x.combined_KC, reduction="pca", group.by="age_condition",dims=c(1,2))+colscale_agecondition_fetal
VlnPlot(d10x.combined_KC, features = 'PC_2', group.by="age_condition")
RidgePlot(d10x.combined_KC, features = 'PC_2', group.by="age_condition")


                            # gene_rank_in_PC<-function(genes,loadings){
                            #   print(sum(sample(1:nrow(loadings), length(genes))))
                            #   sapply(1:ncol(loadings), function(x){
                            #     order_PC<-rownames(loadings)[order(loadings[,x])]
                            #     sum(as.numeric(which(order_PC%in%genes)))
                            #   })
                            # }
                            # 
                            # gene_rank_in_PC(genes, loadings)

##########################
#' ## PCA for batch effect
##########################
Loadings<-as.data.frame(Loadings(d10x.combined_KC, reduction = "pca"))
embed<-as.data.frame(Embeddings(d10x.combined_KC, reduction = "pca"))
vars <- Stdev(d10x.combined_KC, reduction = "pca")^2
Importance<-vars/sum(vars)
print(Importance[1:10])


meta_categorical <- d10x.combined_KC@meta.data[, c("chemistry","Sex","age_condition","Phase")]  # input column numbers in meta that contain categorical variables
meta_continuous <- d10x.combined_KC@meta.data[, c("percent.mt","nFeature_RNA")]  # input column numbers in meta that contain continuous variables


ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-30

suppressWarnings(heat_scree_plot(embed, Importance, 2, 2))



##########################
#' varimax PCA
##########################

## function taken from Delaram
expression<-GetAssayData(d10x.combined_KC)

df<-get_varimax_rotated(as.matrix(expression), as.matrix(loadings))

lapply(1:10, function(x) rownames(df[[1]])[order(df[[1]][,x])][1:10])
lapply(1:10, function(x) rownames(df[[1]])[rev(order(df[[1]][,x]))][1:10])

DimPlot(d10x.combined_KC, reduction="pca", group.by="age_condition",dims=c(1,2))+colscale_agecondition_fetal

plt_varimax<-data.frame(as.matrix(df$rotScores))
plt_varimax_meta<-cbind(plt_varimax,d10x.combined_KC@meta.data)

plt_varimax_meta$age_condition<-factor(plt_varimax_meta$age_condition, levels=c("Fetal Healthy","Ped Healthy","Ped IFALD","Adult Healthy"))

pca_plt<-ggplot()+
  geom_point(aes(X1, X2, color=age_condition),plt_varimax_meta[order(plt_varimax_meta$age_condition),],size=0.75)+
  colscale_agecondition_fetal+
    th_present+theme_bw()+xlab("Varimax PC1")+ylab("Varimax PC2")+
  theme(legend.position = "none")+ylim(-16,5)
pca_plt

ridge_plt<-ggplot(plt_varimax_meta, aes(X2, age_condition, fill=age_condition))+geom_density_ridges()+
  fillscale_agecondition_fetal+theme_bw()+coord_flip()+ theme(axis.title.y = element_blank())+xlim(-16,5)

plot_grid(pca_plt, ridge_plt)

# ggplot(plt_varimax_meta, aes(age_condition, X2))+geom_violin(fill="lightgrey", color="lightgrey")+
#   geom_boxplot(aes(fill=age_condition),width=0.1, outlier.shape = NA)+
#   fillscale_agecondition_fetal+theme_bw()+ylab("Varimax PC2")
# 
# ggplot(plt_varimax_meta, aes(age_condition, X2))+
#   geom_violin(aes(fill=age_condition))+
#   geom_boxplot(fill="lightgrey", width=0.1,outlier.shape = NA)+
#   fillscale_agecondition_fetal+theme_bw()
# 
# ggplot(plt_varimax_meta, aes(X2))+geom_density(aes(fill=age_condition))+fillscale_agecondition_fetal+theme_bw()+facet_wrap(~age_condition)
# 
# ggplot(plt_varimax_meta, aes(age_condition, X2))+geom_violin(aes(fill=age_condition),scale="area")+fillscale_agecondition_fetal+theme_bw()
# 
# ggplot(plt_varimax_meta, aes(X1, X2, color=age_condition))+geom_point()+colscale_agecondition_fetal+facet_wrap(~age_condition)
# 


summary(aov(plt_varimax_meta$X2~plt_varimax_meta$age_condition))

identical(rownames(embed),rownames(plt_varimax_meta))
cor(embed$PC_2, plt_varimax_meta$X2)

summary(aov(embed$PC_2~plt_varimax_meta$age_condition))


### varimax heat scree
vars <- Stdev(d10x.combined_KC, reduction = "pca")^2
Importance<-vars/sum(vars)
ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10
suppressWarnings(heat_scree_plot(df[[2]], Importance, 2, 2))




#############
## GSEA of top genes in PC2
#############
source("scripts/00_GSEA_function.R")
#http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/
GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")


### PC5 loadings
gene_list = as.data.frame(unclass(df$rotLoadings))$PC_2
names(gene_list) = rownames(as.data.frame(unclass(df$rotLoadings)))
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

res = GSEA(gene_list, GO_file, pval = 0.05)

plt_path<-res$Results
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Enrichment_Cell<-"Down in PC2"
plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up in PC2"

plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))

plt_path$direction_label<-as.factor(plt_path$Enrichment)
levels(plt_path$direction_label)<-c(0.1,-0.1)
plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))

# top and bottom 15
plt_path<-rbind(plt_path[1:5,], plt_path[(nrow(plt_path)-5):(nrow(plt_path)),])

gsea_PC2<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_hline(yintercept=30.5, color="grey")+scale_fill_manual(values=c("#D64A56","cornflowerblue"))


FeaturePlot(d10x.combined_KC, reduction="umap",features = c("IL1B"), split.by = "age_condition")
FeaturePlot(d10x.combined_KC, reduction="umap",features = c("HBB"), split.by = "age_condition")


gene_loadings<-data.frame(matrix(as.numeric(df[[1]]), attributes(df[[1]])$dim, dimnames=attributes(df[[1]])$dimnames))
x <- tableGrob( d = data.frame( Gene = rownames(gene_loadings)[order(gene_loadings$PC_2)][1:10], Score = gene_loadings[order(gene_loadings$PC_2),2][1:10] ) )
plt_table <- ggdraw() +
  draw_plot( plot = x, x = 0, y = .5, width = 1, height = .5 )



plot_grid(plot_grid(pca_plt, ridge_plt), plot_grid(gsea_PC2, plt_table,rel_widths = c(3,1)), nrow=2, rel_heights = c(2,1))

plot_grid(plot_grid(pca_plt, ridge_plt,plt_table, nrow=1, rel_widths = c(2,1,0.5)), gsea_PC2, nrow=2, rel_heights = c(2,1))


plot_grid(
  plot_grid(
  plot_grid(pca_plt, ridge_plt+theme(legend.position = "none"), rel_widths = c(2,1), align="vh", axis="t"), 
  plot_grid(plt_table,get_legend(ridge_plt), get_legend(gsea_PC2), ncol=1), rel_widths = c(4,1)),
  gsea_PC2+theme(legend.position = "none"),rel_heights = c(3,1), nrow=2)
  
  


# 
# #########################
# ## random forest
# #########################
# library(caTools)
# library(randomForest)
# 
# colnames(plt_varimax) <- paste0('varPC_', 1:ncol(plt_varimax))
# identical(rownames(d10x.combined_KC@meta.data), rownames(plt_varimax))
# plt_varimax$age_group <-as.factor(d10x.combined_KC@meta.data$age_condition)
# table(plt_varimax$age_group)
# 
# plt_varimax_adult_fetal<-plt_varimax[which(plt_varimax$age_group%in%c("Adult Healthy", "Fetal Healthy")),]
# plt_varimax_fetal_ped<-plt_varimax[which(plt_varimax$age_group%in%c("Ped Healthy", "Fetal Healthy")),]
# 
# plt_varimax_adult_fetal$age_group <-as.factor(as.character(plt_varimax_adult_fetal$age_group))
# plt_varimax_fetal_ped$age_group <-as.factor(as.character(plt_varimax_fetal_ped$age_group))
# 
# 
# #############
# ## Adult Fetal RF
# #############
# 
# ### splitting the train and test data
# sample = sample.split(plt_varimax_adult_fetal$age_group, SplitRatio = .75)
# train = subset(plt_varimax_adult_fetal, sample == TRUE)
# test = subset(plt_varimax_adult_fetal, sample == FALSE)
# 
# ### training the RF model
# age_group_pred_RF <- randomForest(age_group ~ . , data = train, importance = TRUE)
# 
# 
# pred = predict(age_group_pred_RF, newdata=test[,-ncol(test)])
# ### generating a confusion matrix
# cm = table(label=test[,ncol(test)], prediction=pred)
# cm
# gridExtra::grid.table(cm)
# #### evaluating feature importance
# imp.df = data.frame(importance(age_group_pred_RF))
# imp.df[order(imp.df$MeanDecreaseAccuracy, decreasing = T),]
# dev.off()
# varImpPlot(age_group_pred_RF,main = 'age_group Prediction Based on Varimax-PCs')
# 
# ggplot(plt_varimax_meta, aes(X1, X23, color=age_condition))+geom_point()+colscale_agecondition_fetal
# ggplot(plt_varimax_meta, aes(age_condition, X4))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=age_condition))+fillscale_agecondition_fetal+theme_bw()
# 
# ggplot(plt_varimax_meta, aes(X1, X3, color=age_condition))+geom_point()+colscale_agecondition
# ggplot(plt_varimax_meta, aes(age_condition, X3))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=age_condition))+fillscale_agecondition_fetal+theme_bw()
# 
# ggplot(plt_varimax_meta, aes(X1, X5, color=age_condition))+geom_point()+colscale_agecondition
# ggplot(plt_varimax_meta, aes(age_condition, X5))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=age_condition))+fillscale_agecondition_fetal+theme_bw()
# 
# 
# ggplot(plt_varimax_meta, aes(X1, X4, color=nFeature_RNA))+geom_point()
# 
# #############
# ## Ped Fetal RF
# #############
# 
# ### splitting the train and test data
# sample = sample.split(plt_varimax_fetal_ped$age_group, SplitRatio = .75)
# train = subset(plt_varimax_fetal_ped, sample == TRUE)
# test = subset(plt_varimax_fetal_ped, sample == FALSE)
# 
# ### training the RF model
# age_group_pred_RF <- randomForest(age_group ~ . , data = train, importance = TRUE)
# 
# 
# pred = predict(age_group_pred_RF, newdata=test[,-ncol(test)])
# ### generating a confusion matrix
# cm = table(label=test[,ncol(test)], prediction=pred)
# cm
# gridExtra::grid.table(cm)
# #### evaluating feature importance
# imp.df = data.frame(importance(age_group_pred_RF))
# imp.df[order(imp.df$MeanDecreaseAccuracy, decreasing = T),]
# dev.off()
# varImpPlot(age_group_pred_RF,main = 'age_group Prediction Based on Varimax-PCs')
# 
# ggplot(plt_varimax_meta, aes(X1, X4, color=age_condition))+geom_point()+colscale_agecondition
# ggplot(plt_varimax_meta, aes(age_condition, X4))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=age_condition))+fillscale_agecondition_fetal+theme_bw()
# 
# ggplot(plt_varimax_meta, aes(X1, X3, color=age_condition))+geom_point()+colscale_agecondition
# ggplot(plt_varimax_meta, aes(age_condition, X3))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=age_condition))+fillscale_agecondition_fetal+theme_bw()
# 
# ggplot(plt_varimax_meta, aes(X1, X5, color=age_condition))+geom_point()+colscale_agecondition
# ggplot(plt_varimax_meta, aes(age_condition, X5))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=age_condition))+fillscale_agecondition_fetal+theme_bw()
# 
# 
# ggplot(plt_varimax_meta, aes(X1, X4, color=nFeature_RNA))+geom_point()
# 
# 
### GSEA PC loadings

enrichment_PC<-function(PC){
  gene_list = PC
  names(gene_list) = rownames(as.data.frame(unclass(df$rotLoadings)))
  gene_list = sort(gene_list, decreasing = TRUE)
  gene_list = gene_list[!duplicated(names(gene_list))]

  res = GSEA(gene_list, GO_file, pval = 0.05)

  plt_path<-res$Results
  plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
  plt_path$Enrichment_Cell<-"Down in PC"
  plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up in PC"

  plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))

  plt_path$direction_label<-as.factor(plt_path$Enrichment)
  levels(plt_path$direction_label)<-c(0.1,-0.1)
  plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))

  # top and bottom 15
  if(nrow(plt_path)>30){plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])}

  ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
    theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
    geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
    scale_fill_manual(values=c("#D64A56","cornflowerblue"))
}

enrichment_PC(as.data.frame(unclass(df$rotLoadings))$PC_10)#sex but no genesets
lapply(10, function(x) rownames(df[[1]])[(order(df[[1]][,x]))][1:10])
lapply(10, function(x) rownames(df[[1]])[rev(order(df[[1]][,x]))][1:10])

enrichment_PC(as.data.frame(unclass(df$rotLoadings))$PC_16)# nothing
lapply(3, function(x) rownames(df[[1]])[(order(df[[1]][,x]))][1:10])
lapply(3, function(x) rownames(df[[1]])[rev(order(df[[1]][,x]))][1:10])

enrichment_PC(as.data.frame(unclass(df$rotLoadings))$PC_20) #nothing
lapply(20, function(x) rownames(df[[1]])[(order(df[[1]][,x]))][1:10])
lapply(20, function(x) rownames(df[[1]])[rev(order(df[[1]][,x]))][1:10])

enrichment_PC(as.data.frame(unclass(df$rotLoadings))$PC_5)# IL18


lapply(16, function(x) rownames(df[[1]])[(order(df[[1]][,x]))][1:10])
lapply(16, function(x) rownames(df[[1]])[rev(order(df[[1]][,x]))][1:10])


