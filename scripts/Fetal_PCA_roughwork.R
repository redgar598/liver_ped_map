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
                            

## only KC
d10x.combined_KC<-subset(d10x.combined_myeloid, subset = CellType_harmonized %in% c("KC Like"))
save(d10x.combined_KC, file=here("/media/redgar/Seagate Portable Drive/processed_data/Fetal_IFALD_adult_ped_integrated_KC_only.RData"))

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


meta_categorical <- d10x.combined_KC@meta.data[, c("chemistry","Sex","age_condition")]  # input column numbers in meta that contain categorical variables
meta_continuous <- d10x.combined_KC@meta.data[, c("percent.mt","nFeature_RNA","Age")]  # input column numbers in meta that contain continuous variables

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

ggplot(plt_varimax_meta, aes(X1, X2, color=age_condition))+geom_point()+colscale_agecondition_fetal
ggplot(plt_varimax_meta, aes(age_condition, X2))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=age_condition))+fillscale_agecondition_fetal+theme_bw()


summary(aov(plt_varimax_meta$X2~plt_varimax_meta$age_condition))

identical(rownames(embed),rownames(plt_varimax_meta))
cor(embed$PC_2, plt_varimax_meta$X2)

summary(aov(embed$PC_2~plt_varimax_meta$age_condition))


### heat scree
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
plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])

ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_hline(yintercept=30.5, color="grey")+scale_fill_manual(values=c("#D64A56","cornflowerblue"))

# 
# #########################
# ## random forest
# #########################
# library(caTools)
# library(randomForest)
# 
# colnames(plt_varimax) <- paste0('varPC_', 1:ncol(plt_varimax))
# identical(rownames(d10x.combined_KC@meta.data), rownames(plt_varimax))
# plt_varimax$age_group <-as.factor(d10x.combined_KC@meta.data$AgeGroup)
# table(plt_varimax$age_group)
# 
# ### splitting the train and test data
# sample = sample.split(plt_varimax$age_group, SplitRatio = .75)
# train = subset(plt_varimax, sample == TRUE) 
# test = subset(plt_varimax, sample == FALSE)
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
# ggplot(plt_varimax_meta, aes(age_condition, X4))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=age_condition))+fillscale_agecondition+theme_bw()
# 
# ggplot(plt_varimax_meta, aes(X1, X3, color=age_condition))+geom_point()+colscale_agecondition
# ggplot(plt_varimax_meta, aes(age_condition, X3))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=age_condition))+fillscale_agecondition+theme_bw()
# 
# ggplot(plt_varimax_meta, aes(X1, X5, color=age_condition))+geom_point()+colscale_agecondition
# ggplot(plt_varimax_meta, aes(age_condition, X5))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=age_condition))+fillscale_agecondition+theme_bw()
# 
# 
# ggplot(plt_varimax_meta, aes(X1, X4, color=nFeature_RNA))+geom_point()
# 
# 
# ### GSEA PC loadings
# 
# enrichment_PC<-function(PC){
#   gene_list = PC
#   names(gene_list) = rownames(as.data.frame(unclass(df$rotLoadings)))
#   gene_list = sort(gene_list, decreasing = TRUE)
#   gene_list = gene_list[!duplicated(names(gene_list))]
#   
#   res = GSEA(gene_list, GO_file, pval = 0.05)
#   
#   plt_path<-res$Results
#   plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
#   plt_path$Enrichment_Cell<-"Down in PC"
#   plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up in PC"
#   
#   plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))
#   
#   plt_path$direction_label<-as.factor(plt_path$Enrichment)
#   levels(plt_path$direction_label)<-c(0.1,-0.1)
#   plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))
#   
#   # top and bottom 15
#   plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])
#   
#   ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
#     theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
#     geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
#     scale_fill_manual(values=c("#D64A56","cornflowerblue"))
# }
# 
# enrichment_PC(as.data.frame(unclass(df$rotLoadings))$PC_4)#ALB
# 
# enrichment_PC(as.data.frame(unclass(df$rotLoadings))$PC_3)# nothing
# lapply(3, function(x) rownames(df[[1]])[(order(df[[1]][,x]))][1:10])
# lapply(3, function(x) rownames(df[[1]])[rev(order(df[[1]][,x]))][1:10])
# 
# enrichment_PC(as.data.frame(unclass(df$rotLoadings))$PC_20) #nothing
# lapply(20, function(x) rownames(df[[1]])[(order(df[[1]][,x]))][1:10])
# lapply(20, function(x) rownames(df[[1]])[rev(order(df[[1]][,x]))][1:10])
# 
# enrichment_PC(as.data.frame(unclass(df$rotLoadings))$PC_5)# IL18



