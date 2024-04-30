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
source("scripts/00_fanciest_UMAP.R")
source("scripts/varimax_delaram.R")


##############
## All Myeloid
##############
load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/Fetal_IFALD_adult_ped_integrated_myeloid_only.RData"))

# only KC
d10x.combined_KC<-subset(d10x.combined_myeloid, subset = CellType_harmonized %in% c("KC Like"))
save(d10x.combined_KC, file=here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/Fetal_IFALD_adult_ped_integrated_KC_only.RData"))
#save(d10x.combined_KC, file=here("../../../projects/macparland/RE/PediatricAdult/processed_data","Fetal_IFALD_adult_ped_integrated_KC_only.RData"))







##############
## Just KC
##############
load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/Fetal_IFALD_adult_ped_integrated_KC_only.RData"))

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

lapply(1:20, function(x) rownames(df[[1]])[order(df[[1]][,x])][1:10])
lapply(1:20, function(x) rownames(df[[1]])[rev(order(df[[1]][,x]))][1:10])

DimPlot(d10x.combined_KC, reduction="pca", group.by="age_condition",dims=c(1,2))+colscale_agecondition_fetal

plt_varimax<-data.frame(as.matrix(df$rotScores))
plt_varimax_meta<-cbind(plt_varimax,d10x.combined_KC@meta.data)

plt_varimax_meta$age_condition<-factor(plt_varimax_meta$age_condition, levels=c("Fetal Healthy","Ped Healthy","Ped IFALD","Adult Healthy"))




pca_plt<-ggplot()+
  geom_point(aes(X2, X19, color=age_condition),plt_varimax_meta[order(plt_varimax_meta$age_condition),],size=0.75)+
  colscale_agecondition_fetal+xlim(-15,3.5)+ylim(-9,8.5)+
    theme_bw()+xlab("Varimax PC2")+ylab("Varimax PC20")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  guides(colour = guide_legend(override.aes = list(size=5)))

ridge_plt<-ggplot(plt_varimax_meta, aes(X19, age_condition, fill=age_condition))+geom_density_ridges()+
  fillscale_agecondition_fetal+theme_void()+coord_flip()+xlim(-7,8.5)+
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.text=element_blank(),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

ridge_plt_top<-ggplot(plt_varimax_meta, aes(X2, age_condition, fill=age_condition))+geom_density_ridges()+
  fillscale_agecondition_fetal+theme_void()+xlim(-15,3.5)+
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        plot.margin = unit(c(0.1, 0, 0, 0), "cm"))

plot_grid(ridge_plt_top, get_leg(pca_plt), pca_plt+theme(legend.position = "none"), ridge_plt, align="hv", axis="lrtb", rel_widths = c(2,1), rel_heights = c(1,2))




summary(aov(plt_varimax_meta$X2~plt_varimax_meta$age_condition))
summary(aov(plt_varimax_meta$X20~plt_varimax_meta$age_condition))


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


### PC2 loadings
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
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score (PC2 Gene Scores)")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=2)+
  geom_hline(yintercept=30.5, color="grey")+scale_fill_manual(values=c("#D64A56","cornflowerblue"))+
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=12))




gene_loadings<-data.frame(matrix(as.numeric(df[[1]]), attributes(df[[1]])$dim, dimnames=attributes(df[[1]])$dimnames))
PC2_table <- tableGrob( d = data.frame( Gene = rownames(gene_loadings)[order(gene_loadings$PC_2)][1:10], Score = gene_loadings[order(gene_loadings$PC_2),2][1:10] ),rows = NULL ,theme=ttheme_default(base_size = 7))
PC2_table <- ggdraw() +  draw_plot( plot = PC2_table, x = 0, y = 0, width = 1, height = 1)
PC2_title <- ggdraw() +  draw_label("PC2 Top Genes", fontface='bold')

PC19_table <- tableGrob( d = data.frame( Gene = rownames(gene_loadings)[rev(order(gene_loadings$PC_19))][1:10], Score = gene_loadings[rev(order(gene_loadings$PC_19)),19][1:10] ), rows = NULL ,theme=ttheme_default(base_size = 7))
PC19_table <- ggdraw() +  draw_plot( plot = PC19_table, x = 0, y = 0, width = 1, height = 1 )
PC19_title <- ggdraw() +draw_label("PC19 Top Genes", fontface='bold')

PC2_table<-plot_grid(PC2_title, PC2_table, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
PC19_table<-plot_grid(PC20_title, PC19_table, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
PC2_table
PC19_table
  

#################
## plot gene varimax
#################
gene="HBA1"
percentile=0.9

varimax_gene_plot<-function(gene, percentile){
  plt_varimax_meta$cell<-rownames(plt_varimax_meta)
  DefaultAssay(d10x.combined_KC)<-"RNA"
  gene_exp<-FetchData(d10x.combined_KC, vars=gene)
  gene_exp$cell<-rownames(gene_exp)
  plt_varimax_meta<-merge(plt_varimax_meta, gene_exp, by='cell')
  
  exp_limit<-quantile(plt_varimax_meta[, which(colnames(plt_varimax_meta)==gene)], percentile)
  plt_varimax_meta$gene_exp_limited<-NA
  over_limit<-which(plt_varimax_meta[, which(colnames(plt_varimax_meta)==gene)]>exp_limit)
  plt_varimax_meta$gene_exp_limited[over_limit]<-plt_varimax_meta[over_limit, which(colnames(plt_varimax_meta)==gene)]
  plt_varimax_meta<-plt_varimax_meta[rev(order(plt_varimax_meta$gene_exp_limited)),]
  
  plt_varimax_meta<-rbind(plt_varimax_meta[which(is.na(plt_varimax_meta$gene_exp_limited)),],
                          plt_varimax_meta[which(!(is.na(plt_varimax_meta$gene_exp_limited))),][(order(plt_varimax_meta[which(!(is.na(plt_varimax_meta$gene_exp_limited))),]$gene_exp_limited)),])
  
  ggplot(plt_varimax_meta, aes(X2,X19))+
    geom_point(aes(color=gene_exp_limited),size=0.75)+xlab("Varimax PC2")+ylab("Varimax PC20")+
    scale_color_continuous_sequential(palette = "Blues 3", rev=T, 
                                      name=paste(gene, "\nExpression)"),na.value = "grey80")+
    theme_bw()}

varimax_gene_plot("HBA1", 0.8)
varimax_gene_plot("IL1B", 0.8)
varimax_gene_plot("TMCC2", 0.8)
varimax_gene_plot("CCL4", 0.8)


### combined plot
scatters<-plot_grid(ridge_plt_top, get_leg(pca_plt), pca_plt+theme(legend.position = "none"), ridge_plt, align="hv", axis="lrtb", rel_widths = c(2,1), rel_heights = c(1,2))
genes<-plot_grid(varimax_gene_plot("IL1B", 0.8), varimax_gene_plot("HBA1", 0.8))
gene_tables<-plot_grid(PC2_table, PC19_table, ncol=2)
fetal_ped<-plot_grid(plot_grid(scatters,gene_tables, ncol=2, rel_widths = c(2,1)),genes, gsea_PC2+theme(legend.position = "none"),nrow=3, rel_heights = c(2,1.2,1.2))

save_fetal_plts(fetal_ped, "fetal_ped_varimax_PCA", w=9,h=12)


###################
## RBC PC
###################
sapply(1:30, function(x) which(rownames(df[[1]])[rev(order(df[[1]][,x]))]=="HBA1"))
sapply(1:30, function(x) which(rownames(df[[1]])[(order(df[[1]][,x]))]=="HBA1"))

lapply(20, function(x) rownames(df[[1]])[order(df[[1]][,x])][1:10])
lapply(20, function(x) rownames(df[[1]])[rev(order(df[[1]][,x]))][1:10])

### PC20 loadings
gene_list = as.data.frame(unclass(df$rotLoadings))$PC_19
names(gene_list) = rownames(as.data.frame(unclass(df$rotLoadings)))
gene_list = sort(gene_list, decreasing = FALSE)
gene_list = gene_list[!duplicated(names(gene_list))]

res = GSEA(gene_list, GO_file, pval = 0.05)

res$Results




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




