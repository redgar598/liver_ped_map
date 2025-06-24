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
library(cowplot)



###########################
## FLASHMM
###########################
#devtools::install_github("https://github.com/Baderlab/FLASHMM", build_vignettes = F)
library(FLASHMM)

          #tx cancelled in some or transferred from UHN so no follow-up
          seu_cleaned_UHN_outcome_rejection<-subset(seu_cleaned_UHN_outcome, subset = Post.LT.Rejection_simplified %in% c("Yes","None","Mild"))
          seu_cleaned_UHN_outcome_rejection <- FindVariableFeatures(seu_cleaned_UHN_outcome_rejection, selection.method = "vst", nfeatures = 2000)
          
          
          #seu_cleaned_UHN_outcome_rejection <- NormalizeData(seu_cleaned_UHN_outcome_rejection,scale.factor = 10000, normalization.method = "LogNormalize")
          
          var.genes<-VariableFeatures(seu_cleaned_UHN_outcome_rejection)
          
          counts <- as.matrix(seu_cleaned_UHN_outcome_rejection@assays$RNA@counts)
          counts<-counts[var.genes,]
          metadata <- seu_cleaned_UHN_outcome_rejection@meta.data
          metadata$libsize<-colSums(counts)
          
          metadata$Post.LT.Rejection_simplified<-as.factor(metadata$Post.LT.Rejection_simplified)
          levels(metadata$Post.LT.Rejection_simplified)<-c(2,1,3)
          metadata$Post.LT.Rejection_simplified<-as.numeric(as.character(metadata$Post.LT.Rejection_simplified))
          
          Y <- log(counts + 1) 
          X <- model.matrix(~ 0 + log(libsize) + Post.LT.Rejection_simplified + Gamma.Annotation:Post.LT.Rejection_simplified, data = metadata)
          Z <- model.matrix(~ 0 + donor, data = metadata)
          d <- ncol(Z)
          
          fit <- lmmfit(Y, X, Z, d = d)
          
          test <- lmmtest(fit)
          all(t(fit$t) == test[, grep("_t", colnames(test))])
          
          fit$t[, 1:5]
          
          ### pvalues
          all(t(fit$p) == test[, grep("_p", colnames(test))])
          fit$p[, 1:5]
          


##############
## Differential expression with age in KC
##############
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_d10x_adult_ped_raw.rds"))

load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
d10x$Sex[which(d10x$individual%in%c("C113","C115"))]<-"M"

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")


d10x_raw_KC<-subset(d10x, subset = CellType_refined %in% c("KC Like"))

Idents(d10x_raw_KC)<-d10x_raw_KC$age_condition
table(d10x_raw_KC$age_condition)

## age differential
de<-FindMarkers(d10x_raw_KC, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]
sig_de[which(sig_de$avg_log2FC>0),]
sig_de[which(sig_de$avg_log2FC<0),]

#################
## FlashMM
#################
## age differential mixed model
d10x_raw_KC_healthy<-subset(d10x_raw_KC, subset = age_condition %in% c("Adult Healthy", "Ped Healthy"))

d10x_raw_KC_healthy<-FindVariableFeatures(d10x_raw_KC_healthy)
var.genes<-VariableFeatures(d10x_raw_KC_healthy)

counts <- as.matrix(d10x_raw_KC_healthy@assays$RNA@counts)
counts<-counts[var.genes,]
metadata <- d10x_raw_KC_healthy@meta.data
metadata$libsize<-colSums(counts)

metadata$age_condition<-as.character(metadata$age_condition)

Y <- log(counts + 1) 
X <- model.matrix(~ 0 + log(libsize) +  age_condition + Sex , data = metadata)
Z <- model.matrix(~ 0 + file, data = metadata)
d <- ncol(Z)

fit <- lmmfit(Y, X, Z, d = d)

### pvalues
age_KC_coef<-data.frame(gene=colnames(fit$p),pval=fit$p["age_conditionPed Healthy",], 
                        coef_ped=fit$coef["age_conditionPed Healthy",],
                        coef_adult=fit$coef["age_conditionAdult Healthy",])
age_KC_coef$fdr<-p.adjust(age_KC_coef$pval, method="fdr")

age_KC_coef$logFC<-age_KC_coef$coef_adult-age_KC_coef$coef_ped
age_KC_coef$FC <- exp(age_KC_coef$logFC) # Natural FC
age_KC_coef$log2FC <- age_KC_coef$logFC / log(2) 

age_KC_sig<-age_KC_coef[which(age_KC_coef$fdr<0.005),]

intersect(rownames(sig_de),age_KC_sig$gene)


### testing fold change
de$gene<-rownames(de)
de_test<-merge(de, age_KC_coef, by="gene")

ggplot(de_test, aes(avg_log2FC, log2FC))+geom_point()
