#'---
#'title: scRNAseq Differential Expression
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---



#'### Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(here)
library(ggplot2)
library(reshape2)
library(gridExtra)
#library(limma)
library(cowplot)
library(gtools)
#library(ggsignif)
library(scales)


options(stringsAsFactors = FALSE)

source("R_functions/pretty_plots.R")



## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","d10x_adult_ped_raw.rds"))


######
## add cell type labels
######
load(here("data","adult_ped_cellRefined.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

d10x_RR_KC<-subset(d10x, subset = CellType_refined %in% c("RR Myeloid","KC Like"))



##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x_RR_KC <- NormalizeData(d10x_RR_KC,scale.factor = 10000, normalization.method = "LogNormalize")


## testing factor
d10x_RR_KC$cell_age<-paste(d10x_RR_KC$CellType_refined, d10x_RR_KC$AgeGroup, sep = "_")
Idents(d10x_RR_KC) <- "cell_age"

d10x_RR_KC$CellType_refined<-as.character(d10x_RR_KC$CellType_refined)

table(d10x_RR_KC$CellType_refined, d10x_RR_KC$AgeGroup)


#MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene,
#consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
#and linear regression for the continuous process (i.e., the expression level). 

cell_types<-unique(d10x_RR_KC$CellType_refined)

contrasts_celltype_age<-do.call(rbind,lapply(1:length(cell_types), function(x){
  combinations(n = 2, r = 2, v = d10x_RR_KC$cell_age[grep(cell_types[x],d10x_RR_KC$cell_age)], repeats.allowed = FALSE)}))

contrasts_celltype_age

nrow(contrasts_celltype_age)

###########
## Monte carlo the DE
###########


d10x_adult<-subset(d10x_RR_KC, subset = AgeGroup == "Adult")
ncol(d10x_adult)
d10x_ped<-subset(d10x_RR_KC, subset = AgeGroup == "Ped")
ncol(d10x_ped)


      # ### paralize
      # command_args <- commandArgs(trailingOnly = TRUE)
      # cell_type_indx <- as.numeric(command_args[1])
      # cell_type<-cell_types[cell_type_indx]
      # 



samp_num=100


DE_monte_carlo<-lapply(cell_types, function(cell_type){
  
  contrasts_celltype<-contrasts_celltype_age[grep(cell_type, contrasts_celltype_age)]
  
  d10x_adult_celltype<-subset(d10x_adult, subset = CellType_refined == cell_type)
  ncol(d10x_adult_celltype)
  d10x_ped_celltype<-subset(d10x_ped, subset = CellType_refined == cell_type)
  ncol(d10x_ped_celltype)

  
  de_lists<-sapply(1:samp_num, function(x){
    set.seed(x)
    
    ## make downsampled
    if(ncol(d10x_adult_celltype)<ncol(d10x_ped_celltype)){
      ped_cells_random <- d10x_ped_celltype[, sample(colnames(d10x_ped_celltype), size = ncol(d10x_adult_celltype), replace=F)]
      d10_DE<-merge(d10x_adult_celltype, ped_cells_random)
    }else{
      adult_cells_random <- d10x_adult_celltype[, sample(colnames(d10x_adult_celltype), size = ncol(d10x_ped_celltype), replace=F)]
      d10_DE<-merge(d10x_ped_celltype, adult_cells_random)}
    
    ## run DE 
    de<-FindMarkers(d10_DE, ident.1 = contrasts_celltype[1], ident.2 = contrasts_celltype[2], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
    print(paste(contrasts_celltype[1],"vs", contrasts_celltype[2],":", nrow(de), sep=" "))
    de$gene<-rownames(de)
    rownames(de)<-NULL
    de<-de[,c(6,1:5)]
    de$cell.1<-contrasts_celltype[1]
    de$cell.2<-contrasts_celltype[2]
    
    de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]$gene
    })
  
  sig_gene_count<-unlist(de_lists)
  if(length(sig_gene_count)==0){NA}else{
  sig_gene_count<-as.data.frame(table(sig_gene_count))
  colnames(sig_gene_count)<-c("gene","sig_count")
  
  sig_gene_count$monte_carlo_sig<-sapply(1:nrow(sig_gene_count),function(x){
    1-(sig_gene_count$sig_count[x]+1)/(samp_num+1)
  })
  
  sig_gene_count$cell<-cell_type
  sig_gene_count}})

DE_monte_carlo<-do.call(rbind, DE_monte_carlo)
DE_monte_carlo<-DE_monte_carlo[which(!(is.na(DE_monte_carlo$gene))),]

save(DE_monte_carlo, file=here("data","adult_ped_diff_motecarlo_100_KC_RRonly.RData"))



########
## interpret DGE
########

### get fold change in whole cohort for pathway analysis
## run DE 

diff_exp_all<-lapply(1:nrow(contrasts_celltype_age), function(x){
  de<-FindMarkers(d10x_RR_KC, ident.1 = contrasts_celltype_age[x,1], ident.2 = contrasts_celltype_age[x,2], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
  print(paste(contrasts_celltype_age[x,1],"vs", contrasts_celltype_age[x,2],":", nrow(de), sep=" "))
  de$gene<-rownames(de)
  rownames(de)<-NULL
  de<-de[,c(6,1:5)]
  de$cell.1<-contrasts_celltype_age[x,1]
  de$cell.2<-contrasts_celltype_age[x,2]
  de})


diff_exp_all<-do.call(rbind, diff_exp_all)


#########
## pathway adult versus ped in KC and RR
#########
source("R_functions/GSEA_function_tmp.R")
GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

pathway_plt<-function(de){
  gene_list = de$avg_log2FC
  names(gene_list) = de$gene
  gene_list = sort(gene_list, decreasing = TRUE)
  gene_list = gene_list[!duplicated(names(gene_list))]
  
  res = GSEA(gene_list, GO_file, pval = 0.05)
  
  plt_path<-res$Results
  plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
  plt_path$Enrichment_Cell<-"Up-regulated in \nAdult"
  plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \nPediatric"
  
  plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))
  
  plt_path$direction_label<-as.factor(plt_path$Enrichment)
  levels(plt_path$direction_label)<-c(0.1,-0.1)
  plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))
  
  # top and bottom 15
  plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])
  
  ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
    theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
    geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
    geom_vline(xintercept = 0, color="grey40")+scale_fill_manual(values=c("#fd8d3c","#6baed6"))+ 
    guides(fill = guide_legend(override.aes = list(size=5)))}

RR_GSEA<-pathway_plt(diff_exp_all[which(diff_exp_all$cell.1=="RR Myeloid_Adult"),])
save_plts(RR_GSEA, "GSEA_adult_ped_recently_recruited", w=20,h=10)

KC_GSEA<-pathway_plt(diff_exp_all[which(diff_exp_all$cell.1=="KC Like_Adult"),])
save_plts(KC_GSEA, "GSEA_adult_ped_KClike", w=20,h=10)


## statistically significant
load(file=here("data","adult_ped_diff_motecarlo_100_KC_RRonly.RData"))

sig_MC<-DE_monte_carlo[which(DE_monte_carlo$monte_carlo_sig<0.05),]

diff_exp_all[which(diff_exp_all$gene%in%sig_MC$gene & diff_exp_all$avg_log2FC>0),]
diff_exp_all[which(diff_exp_all$gene%in%sig_MC$gene & diff_exp_all$avg_log2FC<0),]


vol_celltype<-function(cellType){
  diff_exp_all_celltype<-diff_exp_all[grep(cellType, diff_exp_all$cell.1),]
    volcano<-data.frame(gene=diff_exp_all_celltype$gene,Pvalue=diff_exp_all_celltype$p_val, Delta_Beta=diff_exp_all_celltype$avg_log2FC)
    
    #Thresholds 
    dB<-1 #delta beta cutoff
    Pv<-1.8e-07 #Pvalue cutoff
    
    volcano<-volcano[complete.cases(volcano),]
    
    ## positive delta beta is hypomethylated (code for volcano should be right now, should colors change?)
    color3<-sapply(1:nrow(volcano), function(x) if(volcano$Pvalue[x]<=Pv){
      if(abs(volcano$Delta_Beta[x])>dB){
        if(volcano$Delta_Beta[x]>dB){"Higher Expression in Adults\n(with Potential Biological Impact)"}else{"Higher Expression in Pediatric\n(with Potential Biological Impact)"}
      }else{if(volcano$Delta_Beta[x]>0){"Higher Expression in Adults"}else{"Higher Expression in Pediatric"}}}else{"Not Significantly Different"})
    
    volcano$Interesting_CpG3<-color3
    
    
    # COLORS! define here so they are consistent between plots
    # so even if you don't have CpGs in a color catagory the pattern will be maintained
    myColors <- c(muted("red", l=80, c=30),"red",muted("blue", l=70, c=40),"blue", "grey")
    
    color_possibilities<-c("Higher Expression in Adults",
                           "Higher Expression in Adults\n(with Potential Biological Impact)",
                           "Higher Expression in Pediatric",
                           "Higher Expression in Pediatric\n(with Potential Biological Impact)",
                           "Not Significantly Different")
    
    names(myColors) <- color_possibilities
    colscale <- scale_color_manual(name = "Direction of Change",
                                   values = myColors, drop = FALSE)
    fillscale <- scale_fill_manual(name = "Direction of Change",
                                   values = myColors, drop = FALSE)
    
    
    sig_MC_celltype<-sig_MC[grep(cellType,sig_MC$cell),]
    volcano_label<-volcano[which(volcano$gene%in%sig_MC_celltype$gene),]
    volcano$sig<-"not_sig"
    volcano$sig[which(volcano$gene%in%sig_MC_celltype$gene)]<-"MC_sig"
    
    ggplot(volcano, aes(Delta_Beta, -log10(Pvalue), fill=Interesting_CpG3, color=sig))+
      geom_point(shape=21, size=1)+theme_bw()+
      fillscale+scale_color_manual(values=c("black","white"))+guides(color = "none")+
      geom_vline(xintercept=c(-dB,dB), color="grey60")+
      geom_hline(yintercept=-log10(Pv), color="grey60")+
      ylab("P Value (-log10)")+xlab("Differential Expression\n(Fold change)")+
      theme(plot.margin=unit(c(1,1,1,2),"cm"))+ 
      guides(fill = guide_legend(override.aes = list(size = 4)))+
      geom_text(aes(label=gene),volcano_label,color="black",vjust=4, hjust=1,size=3)
    
    ggsave(file=paste(here("figures/"),cellType,"_differential_Ped_adult.pdf", sep=""), h=8, w=10)
    ggsave(file=paste(here("figures/jpeg/"),cellType,"_differential_Ped_adult.jpeg", sep=""), h=8, w=10)}

vol_celltype("RR")
vol_celltype("KC")

diff_exp_all_celltype<-merge(diff_exp_all[grep("RR", diff_exp_all$cell.1),], diff_exp_all[grep("KC", diff_exp_all$cell.1),], by="gene",  suffixes = c("_RR","_KC"))


sig_MC_RR<-sig_MC[grep("RR",sig_MC$cell),]
sig_MC_KC<-sig_MC[grep("KC",sig_MC$cell),]


diff_exp_all_celltype$sig<-sapply(1:nrow(diff_exp_all_celltype), function(x){
  if(diff_exp_all_celltype$gene[x]%in%sig_MC_RR$gene & diff_exp_all_celltype$gene[x]%in%sig_MC_KC$gene){"both"}else{
    if(diff_exp_all_celltype$gene[x]%in%sig_MC_RR$gene){"RR only"}else{
      if(diff_exp_all_celltype$gene[x]%in%sig_MC_KC$gene){"KC only"}else{"NS"}
    }
  }
})

diff_exp_all_celltype_label<-diff_exp_all_celltype[which(diff_exp_all_celltype$gene%in%sig_MC$gene),]


ggplot(diff_exp_all_celltype, aes(avg_log2FC_RR, avg_log2FC_KC, color=sig))+geom_point()+th_present+theme_bw()+
  ylab("KC Differential Expression\n(Fold change)")+xlab("RR Differential Expression\n(log2 Fold change)")+
  scale_color_manual(values = c("red","#f7057d","grey","#b80783"))+
  geom_text(aes(label=gene),diff_exp_all_celltype_label,color="black",vjust=-0.75, hjust=1,size=3)+
  geom_vline(xintercept = c(-1,1), color="grey")+  geom_hline(yintercept = c(-1,1), color="grey")+
  ylim(c(min(diff_exp_all$avg_log2FC),max(diff_exp_all$avg_log2FC)))+
  xlim(c(min(diff_exp_all$avg_log2FC),max(diff_exp_all$avg_log2FC)))+
  annotate("text", x=2, y=3.6, label="Higher Expression in Adults")+
  annotate("text", x=-3.5, y=-3.8, label="Higher Expression in Peds")
ggsave(file=here("figures/FC_correlation_KC_RR_differential_Ped_adult.pdf"), h=7, w=8)
ggsave(file=here("figures/jpeg/FC_correlation_KC_RR_differential_Ped_adult.jpeg"), h=7, w=8)





diff_exp_all[which(diff_exp_all$gene%in%sig_MC_RR$gene),]
