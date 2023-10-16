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
library(RColorBrewer)

library(edgeR)

source("scripts/00_pretty_plots.R")
source(here("scripts/00_edgeR_updated_pseudobulkfunctions.R"))

d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))

load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)



celltypes<-unique(d10x$CellType_refined)


pseudobulk_DGE_age_celltypes<-lapply(1:length(celltypes), function(x){
  d10x_cell<-subset(d10x, subset = CellType_refined %in% celltypes[x])
  print(paste("Pseudobulk DGE in",ncol(d10x_cell),celltypes[x], "cells"))
  
  # fewer than 10 cells in any group won't be tested
  if(any(tapply(d10x_cell$cell, d10x_cell$age_condition, function(y) length(unique(y)))<10)){}else{
  
    ## calculate samplewise pseudo bulk counts
    y <- Seurat2PB(d10x_cell, sample="individual", cluster="age_condition")
            #head(y$samples, n=10L)
            #summary(y$samples$lib.size)
    
    ## keep samples with a library size >5e4
    keep.samples <- y$samples$lib.size > 5e4
    table(keep.samples)
    keep.samples
    y <- y[, keep.samples]
    
    
    ## keep genes with a minimum of 10 cells in a age_condition group and 20 across all samples
    keep.genes <- filterByExpr(y, group=y$samples$cluster, min.count=10, min.total.count=20)
    table(keep.genes)
    y <- y[keep.genes, , keep=FALSE]
    
    # Normalize pseudobulk counts by library sizes
    y <- normLibSizes(y)
            #head(y$samples, n=10L)
            #head(y$counts, n=10L)
    
    
    ### Meta data for DGE
    metadata <- d10x_cell@meta.data[,c("individual","CellType_refined","Treatment","Tissue","chemistry","Sex","Age","AgeGroup", "FreshorFrozen", "Perfused", "BMI","age_condition")]
    metadata<-metadata[!duplicated(metadata),]
    metadata$individual_cell<-paste(metadata$individual, "_cluster",metadata$age_condition, sep="")
    
    metadata<-metadata[match(colnames(y$counts),metadata$individual_cell),]
    metadata$Sex <- as.factor(metadata$Sex) # ensure it is a factor
    
    metadata$age_condition <- as.factor(metadata$age_condition)
    metadata$age_condition <- relevel(metadata$age_condition, ref="Ped Healthy")
    
    # fewer than 3 samples in any healthy age group won't be tested
    sample_group_count<-as.data.frame(tapply(metadata$individual, metadata$age_condition, function(y) length(unique(y))))
    sample_group_count_fewerthan3<-which(sample_group_count[,1]<3)
    if(1 %in%  sample_group_count_fewerthan3 | 3 %in% sample_group_count_fewerthan3){}else{
    
      design <- model.matrix(~ age_condition + Sex, data = metadata)
              #head(design)
      
      y <- estimateDisp(y, design, robust=TRUE)
      fit <- glmQLFit(y, design, robust=TRUE)
      
      
      # Perform differential expression analysis
      de_results <- glmLRT(fit, coef = "age_conditionAdult Healthy") 
      top_de_genes <- topTags(de_results, n = Inf)$table  
      top_de_genes}}})

names(pseudobulk_DGE_age_celltypes)<-celltypes




##########
## Compare MAST single-cell testing to pseudobulk
#########
pseudobulk_DGE_age_celltypes[["KC Like"]][which(pseudobulk_DGE_age_celltypes[["KC Like"]]$logFC>0),][1:10,]
pseudobulk_DGE_age_celltypes[["KC Like"]][which(pseudobulk_DGE_age_celltypes[["KC Like"]]$FDR<0.005 & pseudobulk_DGE_age_celltypes[["KC Like"]]$logFC>2),]
pseudobulk_DGE_age_celltypes[["KC Like"]][which(pseudobulk_DGE_age_celltypes[["KC Like"]]$FDR<0.005 & pseudobulk_DGE_age_celltypes[["KC Like"]]$logFC<=-2),]

pseudobulk_DGE_age_celltypes[["KC Like"]][grep("ALB",rownames(pseudobulk_DGE_age_celltypes[["KC Like"]])),]

pseudobulk_DGE_age_celltypes[["KC Like"]][grep("CCL",rownames(pseudobulk_DGE_age_celltypes[["KC Like"]])),]
pseudobulk_DGE_age_celltypes[["KC Like"]][grep("IL1B",rownames(pseudobulk_DGE_age_celltypes[["KC Like"]])),]

de_KC<-read.csv(here("data","differential_age_KC.csv"))

intersect(de_KC$X, pseudobulk_DGE_age_celltypes[["KC Like"]][which(pseudobulk_DGE_age_celltypes[["KC Like"]]$FDR<0.01 & pseudobulk_DGE_age_celltypes[["KC Like"]]$logFC>0),"gene"])
intersect(de_KC$X, pseudobulk_DGE_age_celltypes[["KC Like"]][which(pseudobulk_DGE_age_celltypes[["KC Like"]]$FDR<0.01 & pseudobulk_DGE_age_celltypes[["KC Like"]]$logFC<0),"gene"])

plt_FC_cor<-merge(de_KC, pseudobulk_DGE_age_celltypes[["KC Like"]], by.x="X", by.y="gene")
ggplot(plt_FC_cor, aes(avg_log2FC, logFC))+geom_point()+geom_abline(slope=1, intercept = 0)

ggplot(plt_FC_cor, aes(-log10(p_val), -log10(PValue)))+geom_point()+
  ylab("Pseudobulk P value (-log10)")+xlab("Single-cell MAST P value (-log10)")+
  geom_hline(yintercept = -log10(max(plt_FC_cor[which(plt_FC_cor$FDR<0.005),]$PValue)))+
  geom_vline(xintercept = -log10(max(plt_FC_cor[which(plt_FC_cor$p_val_adj<0.005),]$p_val)))+
  theme_bw()
  
  
pseudo_sig<-pseudobulk_DGE_age_celltypes[["KC Like"]][which(pseudobulk_DGE_age_celltypes[["KC Like"]]$FDR<0.01 & abs(pseudobulk_DGE_age_celltypes[["KC Like"]]$logFC) > 1),]



de_KC<-read.csv(here("data","all_gene_stats_age_KC.csv"))
plt_FC_cor_KC<-merge(de_KC, pseudobulk_DGE_age_celltypes[["KC Like"]], by.x="X", by.y="gene")
plt_FC_cor_KC$cell<-"KC Like"
de_RR<-read.csv(here("data","all_gene_stats_age_RR.csv"))
plt_FC_cor_RR<-merge(de_RR, pseudobulk_DGE_age_celltypes[["RR Myeloid"]], by.x="X", by.y="gene")
plt_FC_cor_RR$cell<-"RR Myeloid"
de_MHCII<-read.csv(here("data","all_gene_stats_age_MHCII.csv"))
plt_FC_cor_MHCII<-merge(de_MHCII, pseudobulk_DGE_age_celltypes[["Macrophage\n(MHCII high)"]], by.x="X", by.y="gene")
plt_FC_cor_MHCII$cell<-"Macrophage\n(MHCII high)"

plt_FC_cor<-rbind(plt_FC_cor_MHCII, plt_FC_cor_KC, plt_FC_cor_RR)

ggplot(plt_FC_cor, aes(avg_log2FC, logFC))+geom_point()+geom_abline(slope=1, intercept = 0)

sig_lines<-data.frame(cell=c("Macrophage\n(MHCII high)", "RR Myeloid", "KC Like"), 
                      PValue=c(-log10(max(plt_FC_cor_MHCII[which(plt_FC_cor_MHCII$FDR<0.005),]$PValue)),
                               -log10(max(plt_FC_cor_RR[which(plt_FC_cor_RR$FDR<0.005),]$PValue)),
                               -log10(max(plt_FC_cor_KC[which(plt_FC_cor_KC$FDR<0.005),]$PValue))), 
                      p_val=c(-log10(max(plt_FC_cor_MHCII[which(plt_FC_cor_MHCII$p_val_adj<0.005),]$p_val)),
                              -log10(max(plt_FC_cor_RR[which(plt_FC_cor_RR$p_val_adj<0.005),]$p_val)),
                              -log10(max(plt_FC_cor_KC[which(plt_FC_cor_KC$p_val_adj<0.005),]$p_val))))

sig_counts

plt_FC_cor %>% group_by(cell) %>% function(x){data.frame(sigboth=length(which(x$FDR < 0.005 & x$p_val_adj < 0.005)),
                                             sigpseudo=length(which(x$FDR < 0.005 & x$p_val_adj > 0.005)),
                                             sigMAST=length(which(x$FDR > 0.005 & x$p_val_adj < 0.005)),
                                             notsig=length(which(x$FDR > 0.005 & x$p_val_adj > 0.005)))}

count_per_cell<-do.call(rbind,lapply(c("Macrophage\n(MHCII high)", "RR Myeloid", "KC Like"), function(cell){
  cell_cor<-plt_FC_cor[which(plt_FC_cor$cell==cell),]
  data.frame(cell=cell,
             sigboth=length(which(cell_cor$FDR < 0.005 & cell_cor$p_val_adj < 0.005)),
             sigpseudo=length(which(cell_cor$FDR < 0.005 & cell_cor$p_val_adj > 0.005)),
             sigMAST=length(which(cell_cor$FDR > 0.005 & cell_cor$p_val_adj < 0.005)),
             notsig=length(which(cell_cor$FDR > 0.005 & cell_cor$p_val_adj > 0.005)))}))



myeloid_pseudo<-ggplot(plt_FC_cor, aes(-log10(p_val), -log10(PValue)))+geom_point(size=0.5)+
  ylab("Pseudobulk P value (-log10)")+xlab("Single-cell MAST P value (-log10)")+
  geom_hline(data=sig_lines,aes(yintercept = PValue))+
  geom_vline(data=sig_lines, aes(xintercept = p_val))+
  theme_bw()+facet_wrap(~cell)+
  geom_text(x=200, y=30, aes(label=sigboth), data=count_per_cell)+
  geom_text(x=200, y=-2, aes(label=sigMAST), data=count_per_cell)+
  geom_text(x=-25, y=30, aes(label=sigpseudo), data=count_per_cell)+
  geom_text(x=-25, y=-2, aes(label=notsig), data=count_per_cell)+
  xlim(-50,300)+ylim(-5,45)

myeloid_pseudo

save_plts(myeloid_pseudo, "myeloid_MAST_vs_pseudobulk", 10, 4)

##### Fold Change
myeloid_pseudo_fc<-ggplot(plt_FC_cor, aes(avg_log2FC, logFC))+geom_point(size=0.5)+geom_abline(slope=1, intercept = 0)+
  facet_wrap(~cell)+theme_bw()+  
  geom_hline(yintercept = c(1,-1), color="grey")+
  geom_vline(xintercept = c(1,-1),  color="grey")
save_plts(myeloid_pseudo_fc, "myeloid_MAST_vs_pseudobulk_FC", 10, 4)




## log normalized counts for library size
lcpm <- cpm(y, log=TRUE)
identical(colnames(lcpm),metadata$individual_cell)
plt<-metadata

gene<-"CCL3"
gene<-"AREG"
gene<-"IL1B"

gene<-"KMT2B"
gene<-"SAA2"
gene<-"DHCR24"
gene<-"GNG10"


plt$gene<-lcpm[which(rownames(lcpm)==gene),]

ggplot(plt, aes(age_condition, gene))+geom_boxplot(outlier.shape = NA)+geom_point(position = position_jitter(w=0.25))
ggplot(plt, aes(age_condition, gene))+geom_boxplot()+geom_text(aes(label=individual))




########################################################################################################################### 
## Pseudo bulk all cells
########################################################################################################################### 

d10x_cell<-d10x

# fewer than 10 cells in any group won't be tested
any(tapply(d10x_cell$cell, d10x_cell$age_condition, function(y) length(unique(y)))<10)
    
## calculate samplewise pseudo bulk counts
y <- Seurat2PB(d10x_cell, sample="individual", cluster="age_condition")
#head(y$samples, n=10L)
#summary(y$samples$lib.size)

## keep samples with a library size >5e4
keep.samples <- y$samples$lib.size > 5e4
table(keep.samples)
keep.samples
y <- y[, keep.samples]
    

## keep genes with a minimum of 10 cells in a age_condition group and 20 across all samples
keep.genes <- filterByExpr(y, group=y$samples$cluster, min.count=10, min.total.count=20)
table(keep.genes)
y <- y[keep.genes, , keep=FALSE]

# Normalize pseudobulk counts by library sizes
y <- normLibSizes(y)
#head(y$samples, n=10L)
#head(y$counts, n=10L)

    
### Meta data for DGE
metadata <- d10x_cell@meta.data[,c("individual","CellType_refined","Treatment","Tissue","chemistry","Sex","Age","AgeGroup", "FreshorFrozen", "Perfused", "BMI","age_condition")]
metadata<-metadata[!duplicated(metadata),]
metadata$individual_cell<-paste(metadata$individual, "_cluster",metadata$age_condition, sep="")

metadata<-metadata[match(colnames(y$counts),metadata$individual_cell),]
metadata$Sex <- as.factor(metadata$Sex) # ensure it is a factor

metadata$age_condition <- as.factor(metadata$age_condition)
metadata$age_condition <- relevel(metadata$age_condition, ref="Ped Healthy")

# fewer than 3 samples in any healthy age group won't be tested
sample_group_count<-as.data.frame(tapply(metadata$individual, metadata$age_condition, function(y) length(unique(y))))
which(sample_group_count[,1]<3)

design <- model.matrix(~ age_condition + Sex, data = metadata)
#head(design)

y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)


# Perform differential expression analysis
de_results <- glmLRT(fit, coef = "age_conditionPed IFALD") 
top_de_genes <- topTags(de_results, n = Inf)$table  
top_de_genes


top_de_genes[which(top_de_genes$logFC>0),][1:10,]
top_de_genes[which(top_de_genes$FDR<0.005 & top_de_genes$logFC>2),][1:10,]
top_de_genes[which(top_de_genes$FDR<0.005 & top_de_genes$logFC<=-2),][1:10,]

top_de_genes[grep("IL1B",top_de_genes$gene),]


## log normalized counts for library size
lcpm <- cpm(y, log=TRUE)
identical(colnames(lcpm),metadata$individual_cell)
plt<-metadata
gene<-"IL1B"
plt$gene<-lcpm[which(rownames(lcpm)==gene),]

ggplot(plt, aes(age_condition, gene))+geom_boxplot(outlier.shape = NA)+geom_point(position = position_jitter(w=0.25))