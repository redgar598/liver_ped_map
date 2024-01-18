## Load Libraries
library(here)
library(Seurat)

library(SCINA)
library(reshape2)
library(dplyr)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

library(sp)
library(rgeos)
library(scales)
library(viridis)
library(colorspace)

library(caTools)
library(randomForest)

source("scripts/00_pretty_plots.R")
source("scripts/00_long_functions.R")
source("scripts/00_Heat_scree_plot_generic.R")


## the panel comes with cell type labels
cell_markers<-read.csv(here("data/Xenium_mBrain_v1.1_metadata.csv"))
load(file=here("data/cell_type_labels.RData"))



### varimax roatation fucntion (from Delaram)
get_varimax_rotated <- function(gene_exp_matrix, loading_matrix){
  
  ## gene_exp_matrix: gene expression matrix. rows named as genes and columns named as UMIs.
  ##                  attention: Data SHOULD be SCALED.
  ##                  Potentially gained from Seurat GetAssayData() function
  ## loading_matrix: the PCA loading matrix with rows named as gene names and 
  ##                 columns named as PC numbers. Potentially gained from Seurat Loadings() function
  
  ######### Varimax rotation
  initial_data <- t(gene_exp_matrix[rownames(gene_exp_matrix) %in% rownames(loading_matrix),])
  
  ## apply varimax rotation on the loadings
  varimax_res <- varimax(loading_matrix)
  rotatedLoadings <- varimax_res$loadings
  ## calculating the PC scores matrix
  invLoadings     <- t(pracma::pinv(rotatedLoadings))
  #scores          <- scale(initial_data) %*% invLoadings ## this second scaling is not necessary
  scores          <- initial_data %*% invLoadings ## this second scaling is not necessary
  ## compacting the rotated loading and score matrices in a list
  rotated_data <- list(rotLoadings=rotatedLoadings, rotScores = scores)
  return(rotated_data)
}





#####################################################

smple<-"tJWH051"


load(file=paste(here("data/"),smple, "_object_raw.RData",sep=""))
plt_umap_xenium_sample<-plt_umap_xenium[which(plt_umap_xenium$sample==smple),c("cell","reference_cluster","CellType")]
plt_umap_xenium_sample<-plt_umap_xenium_sample[!duplicated(plt_umap_xenium_sample),]
rownames(plt_umap_xenium_sample)<-sapply(1:nrow(plt_umap_xenium_sample), function(x) {
  strsplit(plt_umap_xenium_sample$cell[x], "_")[[1]][1]})
xenium.obj

plt_umap_xenium_sample<-plt_umap_xenium_sample[match(colnames(xenium.obj), rownames(plt_umap_xenium_sample)),]
identical(colnames(xenium.obj), rownames(plt_umap_xenium_sample))
plt_umap_xenium_sample$cell<-NULL

xenium.obj<-AddMetaData(xenium.obj, plt_umap_xenium_sample)

plot_freely<-cbind(xenium.obj@meta.data,GetTissueCoordinates(xenium.obj))
plot_freely$cell<-NULL

## load shapes
structure_annotation<-read.csv(paste(here("data/"),smple, " Coordinates.csv",sep=""), skip=2)
structure_annotation$Selection<-as.factor(structure_annotation$Selection)
levels(structure_annotation$Selection)

xenium.obj_stroke<-grab_cells("stroke_site",structure_annotation,xenium.obj)
xenium.obj_not_stroke<-grab_cells("contralateral_stroke",structure_annotation,xenium.obj)

xenium.obj_stroke$region<-"stroke_site"
xenium.obj_not_stroke$region<-"contralateral_stroke"
xenium.obj_regions_combined<-merge(xenium.obj_not_stroke, xenium.obj_stroke)


## onyl neurons
xenium.obj_regions_combined_neuron<-subset(xenium.obj_regions_combined, subset = CellType %in% c(
    "L2/3 IT CTX", "L2/3 IT PPP","L2 IT ENTl",   "L2 IT ENTm" ,
    "L3 IT ENT","L2/3 IT RHP",  "L2/3 IT ENTl", "L5 PT CTX","L5/6 NP CTX",
    "L4 RSP-ACA",
    "L4/5 IT CTX",
    "L5 NP CTX","L5 IT CTX", "L5 PPP",
    "L6 CT CTX","L6b/CT ENT","L6b CTX","L6 IT CTX", "L6 IT ENTl",
    "CR", "Car3","Meis2", 
    "Dentate gyrus granule cells",
    "Sst Chodl inteneurons", "Sst Chodl", "Sncg interneurons","Pvalb interneurons" ,"Pvalb",
    "Cajal-Retzius cells","CR","Lamp5 interneurons","Lamp5",
    "NP PPP",   "CA1-ProS","Sst interneurons", "Sst",
    "CA2", "CA3","Vip interneurons","Vip",
    "CT SUB",   "NP SUB",   "SUB-ProS","Sncg",
    "unknown","general_interneuron","lamp5 (low quality)","l2IT ENTl (low quality)"))

xenium.obj_regions_combined_neuron <- NormalizeData(xenium.obj_regions_combined_neuron)
xenium.obj_regions_combined_neuron <- FindVariableFeatures(xenium.obj_regions_combined_neuron, selection.method = "vst", nfeatures = 2000)
xenium.obj_regions_combined_neuron <- ScaleData(xenium.obj_regions_combined_neuron) 
xenium.obj_regions_combined_neuron <- RunPCA(xenium.obj_regions_combined_neuron, npcs = 30, verbose = FALSE)
xenium.obj_regions_combined_neuron <- RunUMAP(xenium.obj_regions_combined_neuron, reduction = "pca", dims = 1:30)

loadings<-as.data.frame(Loadings(xenium.obj_regions_combined_neuron, reduction = "pca"))


#####################
## Regular PCA
#####################

DimPlot(xenium.obj_regions_combined_neuron, reduction="pca", group.by = "region")

head(Embeddings(xenium.obj_regions_combined_neuron, reduction = "pca")[, 1:5])
head(Loadings(xenium.obj_regions_combined_neuron, reduction = "pca")[, 1:5])
head(Stdev(xenium.obj_regions_combined_neuron, reduction = "pca"))

head(loadings[order(loadings$PC_1),1:5])
tail(loadings[order(loadings$PC_1),1:5])
FeaturePlot(xenium.obj_regions_combined_neuron, reduction="pca", features = "Fn1")
FeaturePlot(xenium.obj_regions_combined_neuron, reduction="pca", features = "Nrn1")


DimPlot(xenium.obj_regions_combined_neuron, reduction="pca", group.by =  "CellType")+colscale_cellType


lapply(1:10, function(x) rownames(loadings)[order(loadings[,x])][1:5])
lapply(1:10, function(x) rownames(loadings)[rev(order(loadings[,x]))][1:5])


##########################
#' ## heat scree
##########################
Loadings<-as.data.frame(Loadings(xenium.obj_regions_combined_neuron, reduction = "pca"))
embed<-as.data.frame(Embeddings(xenium.obj_regions_combined_neuron, reduction = "pca"))
vars <- Stdev(xenium.obj_regions_combined_neuron, reduction = "pca")^2
Importance<-vars/sum(vars)
print(Importance[1:10])

meta_categorical <- xenium.obj_regions_combined_neuron@meta.data[, c("region","CellType")]  # input column numbers in meta that contain categorical variables
meta_continuous <- xenium.obj_regions_combined_neuron@meta.data[, c("nCount_Xenium","nFeature_Xenium")]  # input column numbers in meta that contain continuous variables

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-30

suppressWarnings(heat_scree_plot(embed, Importance, 2, 2))



#########################
## random forest (regualr PCA)
#########################
PCA_scores<-embed
identical(rownames(xenium.obj_regions_combined_neuron@meta.data), rownames(PCA_scores))
PCA_scores$region <-as.factor(xenium.obj_regions_combined_neuron@meta.data$region)
table(PCA_scores$region)

ggplot(PCA_scores, aes(PC_1, PC_9, color=region))+geom_point()


### splitting the train and test data
sample = sample.split(PCA_scores$region, SplitRatio = .75)
train = subset(PCA_scores, sample == TRUE) 
test = subset(PCA_scores, sample == FALSE)

### training the RF model
region_pred_RF <- randomForest(region ~ . , data = train, importance = TRUE)

pred = predict(region_pred_RF, newdata=test[,-ncol(test)])
### generating a confusion matrix
cm = table(label=test[,ncol(test)], prediction=pred)
cm
gridExtra::grid.table(cm)
#### evaluating feature importance 
imp.df = data.frame(importance(region_pred_RF))        
imp.df[order(imp.df$MeanDecreaseAccuracy, decreasing = T),]
dev.off()
varImpPlot(region_pred_RF,main = 'region Prediction Based on PCs')   

## Plot PCs
identical(rownames(embed), rownames(xenium.obj_regions_combined_neuron@meta.data))

plt_PCA_meta<-cbind(embed,xenium.obj_regions_combined_neuron@meta.data)

ggplot(plt_PCA_meta, aes(PC_1, PC_9, color=region))+geom_point()
ggplot(plt_PCA_meta, aes(region, PC_1))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=region), width=0.1)+theme_bw()
ggplot(plt_PCA_meta, aes(region, PC_9))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=region), width=0.1)+theme_bw()

ggplot(plt_PCA_meta, aes(PC_1, PC_9, color=nCount_Xenium))+geom_point()
ggplot(plt_PCA_meta, aes(PC_1, PC_9, color=CellType))+geom_point()+colscale_cellType+theme_bw()
ggplot(plt_PCA_meta, aes(PC_1, PC_9, color=CellType))+geom_point()+colscale_cellType+facet_wrap(~region, ncol=1)+theme_bw()



lapply(c(1,9), function(x) rownames(loadings)[order(loadings[,x])][1:5])
lapply(c(1,9), function(x) rownames(loadings)[rev(order(loadings[,x]))][1:5])


FeaturePlot(xenium.obj_regions_combined_neuron, features = "Grik3", reduction="pca", dims = c(1,9), pt.size=1)
FeaturePlot(xenium.obj_regions_combined_neuron, features = "Cux2", reduction="pca", dims = c(1,9), pt.size=1)

FeaturePlot(xenium.obj_regions_combined_neuron, features = "Rmst", reduction="pca", dims = c(1,9), pt.size=1)
FeaturePlot(xenium.obj_regions_combined_neuron, features = "Nr2f2", reduction="pca", dims = c(1,9), pt.size=1)
FeaturePlot(xenium.obj_regions_combined_neuron, features = "Gfap", reduction="pca", dims = c(1,9), pt.size=1)


DimPlot(xenium.obj_regions_combined_neuron, group.by = "region", reduction="pca", dims = c(1,9), pt.size=1)


ggplot(plt_PCA_meta, aes(PC_4, PC_6, color=region))+geom_point()
ggplot(plt_PCA_meta, aes(PC_4, PC_6, color=CellType))+geom_point()+colscale_cellType
ggplot(plt_PCA_meta, aes(region, PC_4))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=region))+theme_bw()
ggplot(plt_PCA_meta, aes(region, PC_6))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=region))+theme_bw()

lapply(c(4,6), function(x) rownames(loadings)[order(loadings[,x])][1:5])
lapply(c(4,6), function(x) rownames(loadings)[rev(order(loadings[,x]))][1:5])

FeaturePlot(xenium.obj_regions_combined_neuron, features = "Gfap", reduction="pca", dims = c(4,6), pt.size=1)
FeaturePlot(xenium.obj_regions_combined_neuron, features = "Rspo1", reduction="pca", dims = c(4,6), pt.size=1)

FeaturePlot(xenium.obj_regions_combined_neuron, features = "Galnt14", reduction="pca", dims = c(4,6), pt.size=1)
FeaturePlot(xenium.obj_regions_combined_neuron, features = "Cldn5", reduction="pca", dims = c(4,6), pt.size=1)
FeaturePlot(xenium.obj_regions_combined_neuron, features = "Meis2", reduction="pca", dims = c(4,6), pt.size=1)
FeaturePlot(xenium.obj_regions_combined_neuron, features = "Rspo2", reduction="pca", dims = c(4,6), pt.size=1)

##########################
#' varimax PCA
##########################

## function taken from Delaram
expression<-GetAssayData(xenium.obj_regions_combined_neuron)

df<-get_varimax_rotated(as.matrix(expression), as.matrix(loadings))

lapply(1:10, function(x) rownames(df[[1]])[order(df[[1]][,x])][1:5])
lapply(1:10, function(x) rownames(df[[1]])[rev(order(df[[1]][,x]))][1:5])


plt_varimax<-data.frame(as.matrix(df$rotScores))
plt_varimax_meta<-cbind(plt_varimax,xenium.obj_regions_combined_neuron@meta.data)
ggplot(plt_varimax_meta, aes(X1, X3, color=region))+geom_point()



### heat scree
vars <- Stdev(xenium.obj_regions_combined_neuron, reduction = "pca")^2
Importance<-vars/sum(vars)

meta_categorical <- xenium.obj_regions_combined_neuron@meta.data[, c("region","CellType")]  # input column numbers in meta that contain categorical variables
meta_continuous <- xenium.obj_regions_combined_neuron@meta.data[, c("nCount_Xenium","nFeature_Xenium")]  # input column numbers in meta that contain continuous variables

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-30
suppressWarnings(heat_scree_plot(df[[2]], Importance, 2, 2))


#########################
## random forest
#########################
library(caTools)
library(randomForest)

colnames(plt_varimax) <- paste0('varPC_', 1:ncol(plt_varimax))
identical(rownames(xenium.obj_regions_combined_neuron@meta.data), rownames(plt_varimax))
plt_varimax$region <-as.factor(xenium.obj_regions_combined_neuron@meta.data$region)
table(plt_varimax$region)

### splitting the train and test data
sample = sample.split(plt_varimax$region, SplitRatio = .75)
train = subset(plt_varimax, sample == TRUE) 
test = subset(plt_varimax, sample == FALSE)

### training the RF model
region_pred_RF <- randomForest(region ~ . , data = train, importance = TRUE)


pred = predict(region_pred_RF, newdata=test[,-ncol(test)])
### generating a confusion matrix
cm = table(label=test[,ncol(test)], prediction=pred)
cm
gridExtra::grid.table(cm)
#### evaluating feature importance 
imp.df = data.frame(importance(region_pred_RF))        
imp.df[order(imp.df$MeanDecreaseAccuracy, decreasing = T),]
dev.off()
varImpPlot(region_pred_RF,main = 'region Prediction Based on Varimax-PCs')   

ggplot(plt_varimax_meta, aes(X1, X9, color=region))+geom_point()
ggplot(plt_varimax_meta, aes(region, X1))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=region))+theme_bw()

ggplot(plt_varimax_meta, aes(X1, X9, color=nCount_Xenium))+geom_point()


lapply(1, function(x) rownames(df[[1]])[(order(df[[1]][,x]))][1:10])
lapply(1, function(x) rownames(df[[1]])[rev(order(df[[1]][,x]))][1:10])

FeaturePlot(xenium.obj_regions_combined_neuron, features = "Neurod6", reduction="pca", dims = c(1,9), pt.size=1)
FeaturePlot(xenium.obj_regions_combined_neuron, features = "Nxph3", reduction="pca", dims = c(1,9), pt.size=1)

DimPlot(xenium.obj_regions_combined_neuron, group.by = "region", reduction="pca", dims = c(1,9), pt.size=1)



###############
## cell specificity of PCs
###############
## rank sum of cell type markers in each varimax PC
loadings<-df[[1]]

gene_rank_in_PC<-function(genes,loadings){
  center_rank<-sum(round(median(1:nrow(loadings))-(length(genes)/2)):round(median(1:nrow(loadings))+(length(genes)/2)))
  sapply(1:ncol(loadings), function(x){
    order_PC<-rownames(loadings)[order(loadings[,x])]
    sum(as.numeric(which(order_PC%in%genes)))-center_rank
  })
}

celltypes<-unique(cell_markers$Annotation)

ranks_celltype<-do.call(rbind,lapply(celltypes, function(celltype){
  genes<-cell_markers[which(cell_markers$Annotation==celltype),]$Genes
  gene_rank_in_PC(genes, loadings)
}))

## negative ranks are meaningful, just at the other end of the PC
ranks_celltype<-as.data.frame(ranks_celltype)

rownames(ranks_celltype)<-celltypes
colnames(ranks_celltype)<-colnames(loadings)

ranks_celltype[1:5,1:5]
ranks_celltype[rev(order(abs(ranks_celltype$PC_1))),1:5]
ranks_celltype[rev(order(abs(ranks_celltype$PC_9))),1:5]


cell_in_regions<-ranks_celltype[which(rownames(ranks_celltype)%in%unique(plt_varimax_meta$CellType)),]
cell_in_regions[rev(order(abs(cell_in_regions$PC_1))),1:5]

ggplot(plt_varimax_meta, aes(X1, X2, color=CellType))+geom_point()+colscale_cellType+theme_bw()+facet_wrap(~region)
ggplot()+
  geom_point(aes(X1, X2), plt_varimax_meta, color="grey")+
  geom_point(aes(X1, X2, color=CellType), plt_varimax_meta[which(plt_varimax_meta$CellType%in%c("L5 NP CTX","VLMC")),])+
  colscale_cellType+theme_bw()+facet_wrap(~region)

ggplot(plt_varimax_meta, aes(X1, X2, color=region))+geom_point()+theme_bw()
ggplot(plt_varimax_meta, aes(X1, color=region))+geom_density()+theme_bw()


loadings[rev(order(loadings$PC_1)),][1:5,1:5]
loadings[(order(loadings$PC_1)),][1:5,1:5]

identical(rownames(plt_varimax_meta), colnames(xenium.obj_regions_combined_neuron))
plt_varimax_meta$gene<-FetchData(xenium.obj_regions_combined_neuron, vars = "Nrep")$Nrep
ggplot(plt_varimax_meta, aes(X1, X2, color=gene))+geom_point()+theme_bw()
cor(plt_varimax_meta$X1, plt_varimax_meta$gene)


cell_in_regions[rev(order(abs(cell_in_regions$PC_9))),6:9]
ggplot()+
  geom_point(aes(X1, X9), plt_varimax_meta, color="grey")+
  geom_point(aes(X1, X9, color=CellType), plt_varimax_meta[which(plt_varimax_meta$CellType%in%c("Lamp5 interneurons","Oligodendrocytes")),])+
  colscale_cellType+theme_bw()
ggplot(plt_varimax_meta, aes(X9, color=region))+geom_density()+theme_bw()


loadings[rev(order(loadings$PC_9)),][1:5,6:9]
loadings[(order(loadings$PC_9)),][1:5,6:9]

plt_varimax_meta$gene<-FetchData(xenium.obj_regions_combined_neuron, vars = "Cntn6")$Cntn6
ggplot(plt_varimax_meta, aes(X1, X9, color=gene))+geom_point()+theme_bw()
cor(plt_varimax_meta$X1, plt_varimax_meta$gene)
ggplot(plt_varimax_meta, aes(X9, gene, color=gene))+geom_point()+theme_bw()


cell_markers[which(cell_markers$Genes=="Cntn6"),]




#########################
## PCA in each cell type seperately
########################
xenium.obj_regions_combined_neuron<-merge(xenium.obj_not_stroke, xenium.obj_stroke)

cell<-"VLMC"

Idents(xenium.obj_regions_combined_neuron) <- "CellType" 

celltypes_present<-names(table(xenium.obj_regions_combined_neuron$CellType)[which(table(xenium.obj_regions_combined_neuron$CellType)>25)])

region_astro_PCs<-lapply(celltypes_present, function(cell){
  print(cell)
  xenium.obj_regions_combined_neuron_cell<-subset(xenium.obj_regions_combined_neuron, idents = cell, invert = FALSE)
  xenium.obj_regions_combined_neuron_cell <- NormalizeData(xenium.obj_regions_combined_neuron_cell)
  xenium.obj_regions_combined_neuron_cell <- FindVariableFeatures(xenium.obj_regions_combined_neuron_cell, selection.method = "vst", nfeatures = 2000)
  xenium.obj_regions_combined_neuron_cell <- ScaleData(xenium.obj_regions_combined_neuron_cell) 
  xenium.obj_regions_combined_neuron_cell <- RunPCA(xenium.obj_regions_combined_neuron_cell, npcs = 30, verbose = FALSE)
  
  loadings<-as.data.frame(Loadings(xenium.obj_regions_combined_neuron_cell, reduction = "pca"))
  
  #' varimax PCA
  ## function taken from Delaram
  expression<-GetAssayData(xenium.obj_regions_combined_neuron_cell)
  df<-get_varimax_rotated(as.matrix(expression), as.matrix(loadings))
  
  loadings<-df[[1]]
  
  ranks_celltype<-do.call(rbind,lapply(celltypes, function(celltype){
    genes<-cell_markers[which(cell_markers$Annotation==celltype),]$Genes
    gene_rank_in_PC(genes, loadings)
  }))
  
  ## negative ranks are meaningful, just at the other end of the PC
  ranks_celltype<-as.data.frame(ranks_celltype)
  
  rownames(ranks_celltype)<-celltypes
  colnames(ranks_celltype)<-colnames(loadings)
  # 
  # unlist(lapply(1:30, function(x){ 
  #   if("Astrocytes"%in%rownames(ranks_celltype)[order(ranks_celltype[,x])][1:10]){
  #     print(paste(cell, ":", x))
  #   }
  #  }))
  
  astorcyte_PC<-unlist(lapply(1:30, function(x){ 
    if("Astrocytes"%in%rownames(ranks_celltype)[order(ranks_celltype[,x])][1:10]){
      print(paste("varPC_", x, sep=""))
    }
  }))
  
  ## random forest
  plt_varimax<-data.frame(as.matrix(df$rotScores))
  plt_varimax_meta<-cbind(plt_varimax,xenium.obj_regions_combined_neuron_cell@meta.data)
  colnames(plt_varimax) <- paste0('varPC_', 1:ncol(plt_varimax))
  identical(rownames(xenium.obj_regions_combined_neuron_cell@meta.data), rownames(plt_varimax))
  plt_varimax$region <-as.factor(xenium.obj_regions_combined_neuron_cell@meta.data$region)
  
  ### splitting the train and test data
  sample = sample.split(plt_varimax$region, SplitRatio = .75)
  train = subset(plt_varimax, sample == TRUE) 
  test = subset(plt_varimax, sample == FALSE)
  
  ### training the RF model
  region_pred_RF <- randomForest(region ~ . , data = train, importance = TRUE)
  pred = predict(region_pred_RF, newdata=test[,-ncol(test)])
  #### evaluating feature importance 
  imp.df = data.frame(importance(region_pred_RF))        
  seperates_region<-rownames(imp.df[order(imp.df$MeanDecreaseAccuracy, decreasing = T),][1:10,])
  
  if(length(intersect(seperates_region, astorcyte_PC))!=0){
    data.frame(cell=cell, pc=intersect(seperates_region, astorcyte_PC))}
})

region_astro_PCs

##
# L4/5 IT CTX varPC_8
# general_interneuron varPC_5
# L6 CT CTX varPC_8
# L6 CT CTX varPC_5
# Lamp5 interneurons varPC_15
# Lamp5 interneurons varPC_25
# Pvalb interneurons  varPC_8
# Pvalb interneurons varPC_10

## explore more indepth
cell<-"general_interneuron"
xenium.obj_regions_combined_neuron_cell<-subset(xenium.obj_regions_combined_neuron, idents = cell, invert = FALSE)
xenium.obj_regions_combined_neuron_cell <- NormalizeData(xenium.obj_regions_combined_neuron_cell)
xenium.obj_regions_combined_neuron_cell <- FindVariableFeatures(xenium.obj_regions_combined_neuron_cell, selection.method = "vst", nfeatures = 2000)
xenium.obj_regions_combined_neuron_cell <- ScaleData(xenium.obj_regions_combined_neuron_cell) 
xenium.obj_regions_combined_neuron_cell <- RunPCA(xenium.obj_regions_combined_neuron_cell, npcs = 30, verbose = FALSE)

loadings<-as.data.frame(Loadings(xenium.obj_regions_combined_neuron_cell, reduction = "pca"))
expression<-GetAssayData(xenium.obj_regions_combined_neuron_cell)
df<-get_varimax_rotated(as.matrix(expression), as.matrix(loadings))
loadings<-df[[1]]

ranks_celltype<-do.call(rbind,lapply(celltypes, function(celltype){
  genes<-cell_markers[which(cell_markers$Annotation==celltype),]$Genes
  gene_rank_in_PC(genes, loadings)}))
ranks_celltype<-as.data.frame(ranks_celltype)
rownames(ranks_celltype)<-celltypes
colnames(ranks_celltype)<-colnames(loadings)

lapply(5, function(x) rownames(df[[1]])[order(df[[1]][,x])][1:5])
lapply(5, function(x) rownames(df[[1]])[rev(order(df[[1]][,x]))][1:5])

cell_markers[which(cell_markers$Genes%in%lapply(5, function(x) rownames(df[[1]])[rev(order(df[[1]][,x]))][1:20])[[1]]),]
cell_markers[which(cell_markers$Genes%in%lapply(5, function(x) rownames(df[[1]])[(order(df[[1]][,x]))][1:20])[[1]]),]

cell_markers[which(cell_markers$Genes%in%lapply(8, function(x) rownames(df[[1]])[rev(order(df[[1]][,x]))][1:10])[[1]]),]
cell_markers[which(cell_markers$Genes%in%lapply(8, function(x) rownames(df[[1]])[(order(df[[1]][,x]))][1:10])[[1]]),]


plt_varimax<-data.frame(as.matrix(df$rotScores))
plt_varimax_meta<-cbind(plt_varimax,xenium.obj_regions_combined_neuron_cell@meta.data)

ggplot(plt_varimax_meta, aes(X1, X5, color=region))+geom_point()
ggplot(plt_varimax_meta, aes(X5, color=region))+geom_density()+coord_flip()
ggplot(plt_varimax_meta, aes(X5, color=region))+geom_density()

ggplot(plt_varimax_meta, aes(X1, X5, color=FetchData(xenium.obj_regions_combined_neuron_cell, vars = "Gli3")$Gli3))+geom_point()
ggplot(plt_varimax_meta, aes(X1, X5, color=FetchData(xenium.obj_regions_combined_neuron_cell, vars = "Mapk4")$Mapk4))+geom_point()



#########################
## random forest
#########################

colnames(plt_varimax) <- paste0('varPC_', 1:ncol(plt_varimax))
identical(rownames(xenium.obj_regions_combined_neuron_cell@meta.data), rownames(plt_varimax))
plt_varimax$region <-as.factor(xenium.obj_regions_combined_neuron_cell@meta.data$region)
table(plt_varimax$region)

### splitting the train and test data
sample = sample.split(plt_varimax$region, SplitRatio = .75)
train = subset(plt_varimax, sample == TRUE) 
test = subset(plt_varimax, sample == FALSE)

### training the RF model
region_pred_RF <- randomForest(region ~ . , data = train, importance = TRUE)


pred = predict(region_pred_RF, newdata=test[,-ncol(test)])
### generating a confusion matrix
cm = table(label=test[,ncol(test)], prediction=pred)
cm
gridExtra::grid.table(cm)
#### evaluating feature importance 
imp.df = data.frame(importance(region_pred_RF))        
imp.df[order(imp.df$MeanDecreaseAccuracy, decreasing = T),]
dev.off()
varImpPlot(region_pred_RF,main = 'region Prediction Based on Varimax-PCs')   



ggplot(plt_varimax_meta, aes(X1, X29, color=region))+geom_point()
ggplot(plt_varimax_meta, aes(X1, X10, color=region))+geom_point()

ggplot(plt_varimax_meta, aes(region, X29))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=region))+theme_bw()

ggplot(plt_varimax_meta, aes(X1, X9, color=nCount_Xenium))+geom_point()


lapply(29, function(x) rownames(df[[1]])[(order(df[[1]][,x]))][1:10])
