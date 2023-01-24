#'---
#'title: scRNAseq signatures
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
library(ggsignif)
library(stats)
library(scales)
library(colorspace)



options(stringsAsFactors = FALSE)

source("scripts/00_pretty_plots.R")



## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","d10x_adult_ped_raw.rds"))


######
## add cell type labels
######
load(here("data","adult_ped_cellRough.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)




##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")
Idents(d10x) <- "AgeGroup"
cell_types<-unique(d10x$CellType_rough)


#########
## Signature Genes
#########
myeloid_immune_supressive<-c("CTSB","CD163","MS4A7","FOLR2","GPNMB","VSIG4","HMOX1","MSR1")
inflammatory_macs<-c("CD74","HLA-DRA","TYROBP","C1QC","HLA-DPA1","HLA-DPB1","LYZ","S100A6")
exhausted_tcells<-c("TOX","PDCD1","LAG3","TNFRSF9","CXCL13","ENTPD1","HAVCR2","CD38")

recent_recruit_myeloid<-c("S100A8","S100A9","CD68","LYZ")
kuffer_signature<-c("VSIG4","MARCO","CD5L","HMOX1")
#kuffer_signature<-c("CD74","HLA-DRA","TYROBP","C1QC","HLA-DPA1","HLA-DPB1","S100A6","MARCO","CD5L")



######
## Score Signatures
######
d10x <- AddModuleScore(
  object = d10x,
  features = list(myeloid_immune_supressive),
  ctrl = 5,
  name = 'myeloid_immune_supressive_score'
)

d10x <- AddModuleScore(
  object = d10x,
  features = list(inflammatory_macs),
  ctrl = 5,
  name = 'inflammatory_macs_score'
)

d10x <- AddModuleScore(
  object = d10x,
  features = list(exhausted_tcells),
  ctrl = 5,
  name = 'exhausted_tcells_score'
)

d10x <- AddModuleScore(
  object = d10x,
  features = list(recent_recruit_myeloid),
  ctrl = 5,
  name = 'recently_recruited_myeloid'
)

d10x <- AddModuleScore(
  object = d10x,
  features = list(kuffer_signature),
  ctrl = 5,
  name = 'kuffer_like_score'
)


score_data<-d10x@meta.data[,c("myeloid_immune_supressive_score1","inflammatory_macs_score1","exhausted_tcells_score1","recently_recruited_myeloid1","kuffer_like_score1")]
rm(d10x)
gc()

### load integrate for UMAP etc
load(here("data","adult_ped_integrated_refinedlabels.rds"))
score_data<-score_data[match(rownames(d10x.combined@meta.data),rownames(score_data)),]
identical(rownames(d10x.combined@meta.data),rownames(score_data))
d10x.combined <- AddMetaData(d10x.combined, metadata = score_data)



######
## plot scores
######
umap_mat<-as.data.frame(Embeddings(object = d10x.combined, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x.combined@meta.data
meta$cell<-rownames(meta)


plt<-merge(meta, umap_mat, by="cell")


d10x.combined_NK_T_B<-subset(d10x.combined, subset = CellType_rough %in% c("NK and T cells"))
d10x.combined_NK_T_B <- RunPCA(d10x.combined_NK_T_B, npcs = 30, verbose = FALSE)
d10x.combined_NK_T_B <- RunUMAP(d10x.combined_NK_T_B, reduction = "pca", dims = 1:30)
umap_mat_Tcell<-as.data.frame(Embeddings(object = d10x.combined_NK_T_B, reduction = "umap"))#
umap_mat_Tcell$cell<-rownames(umap_mat_Tcell)
meta_Tcell<-d10x.combined_NK_T_B@meta.data
meta_Tcell$cell<-rownames(meta_Tcell)
plt_Tcell<-merge(meta_Tcell, umap_mat_Tcell, by="cell")
cell_num_tcell<-as.data.frame(table(plt_Tcell$AgeGroup))
colnames(cell_num_tcell)<-c("AgeGroup","CellCount")

d10x.combined_myeloid<-subset(d10x.combined, subset = CellType_rough %in% c("Myeloid cells"))
d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)
d10x.combined_myeloid <- RunUMAP(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined_myeloid, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined_myeloid@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")
cell_num_myeloid<-as.data.frame(table(plt_myeloid$AgeGroup))
colnames(cell_num_myeloid)<-c("AgeGroup","CellCount")

cell_num_all<-as.data.frame(table(plt$AgeGroup))
colnames(cell_num_all)<-c("AgeGroup","CellCount")

rm(d10x.combined)
rm(d10x.combined_myeloid)
rm(d10x.combined_NK_T_B)
gc()




#########
## Signatures Monte Carlo
#########

#plt_Tcell
#plt_myeloid


########## 
## Myeloid signatures
########## 

plt_myeloid_adult<-plt_myeloid[which(plt_myeloid$AgeGroup=="Adult"),]
adult_cellcount<-nrow(plt_myeloid_adult)
print(paste(adult_cellcount, " adult myeloid cells",sep=""))

plt_myeloid_ped<-plt_myeloid[which(plt_myeloid$AgeGroup=="Ped"),]
ped_cellcount<-nrow(plt_myeloid_ped)
print(paste(ped_cellcount, " ped myeloid cells",sep=""))


samp_num<-10000
pval<-0.05
myeloid_pval_montecarlo<-do.call(rbind, lapply(1:samp_num, function(x){
  set.seed(x)
  plt_myeloid_ped_random<-plt_myeloid_ped[sample(ped_cellcount, adult_cellcount),]

  supressive<-t.test(plt_myeloid_ped_random$myeloid_immune_supressive_score1, plt_myeloid_adult$myeloid_immune_supressive_score1)$p.value
  inflammatory<-t.test(plt_myeloid_ped_random$inflammatory_macs_score1, plt_myeloid_adult$inflammatory_macs_score1)$p.value
  recruit<-t.test(plt_myeloid_ped_random$recently_recruited_myeloid1, plt_myeloid_adult$recently_recruited_myeloid1)$p.value
  kuffer<-t.test(plt_myeloid_ped_random$kuffer_like_score1, plt_myeloid_adult$kuffer_like_score1)$p.value
  
  data.frame(supressive=supressive, inflammatory=inflammatory, recruit=recruit ,kuffer=kuffer)
  }))

print(paste("Comparing the myeloid immune supressive score in myeloid cells, there is a sig difference between ped and adult (p <", pval,
            ") in ", samp_num, " random samples at a p value of ",
            round((length(myeloid_pval_montecarlo$supressive[which(myeloid_pval_montecarlo$supressive>pval)])+1)/(samp_num+1), 3), sep=""))
print(paste("Comparing the inflammatory macs score in myeloid cells, there is a sig difference between ped and adult (p <", pval,
            ") in ", samp_num, " random samples at a p value of ",
            round((length(myeloid_pval_montecarlo$inflammatory[which(myeloid_pval_montecarlo$inflammatory>pval)])+1)/(samp_num+1), 3), sep=""))
print(paste("Comparing the recently recruited myeloid score in myeloid cells, there is a sig difference between ped and adult (p <", pval,
            ") in ", samp_num, " random samples at a p value of ",
            round((length(myeloid_pval_montecarlo$recruit[which(myeloid_pval_montecarlo$recruit>pval)])+1)/(samp_num+1), 3), sep=""))
print(paste("Comparing the kuffer-like score in myeloid cells, there is a sig difference between ped and adult (p <", pval,
            ") in ", samp_num, " random samples at a p value of ",
            round((length(myeloid_pval_montecarlo$kuffer[which(myeloid_pval_montecarlo$kuffer>pval)])+1)/(samp_num+1), 3), sep=""))


########## 
## T cell Exhaustion
########## 
lapply(c("NK_T"), function(cell_type){
  plt_Tcell_celltype<-plt_Tcell[which(plt_Tcell$CellType_rough==cell_type),]
  
  plt_Tcell_adult<-plt_Tcell_celltype[which(plt_Tcell_celltype$AgeGroup=="Adult"),]
  adult_cellcount<-nrow(plt_Tcell_adult)
  print(paste(adult_cellcount, " adult ",cell_type,sep=""))
  
  plt_Tcell_ped<-plt_Tcell_celltype[which(plt_Tcell_celltype$AgeGroup=="Ped"),]
  ped_cellcount<-nrow(plt_Tcell_ped)
  print(paste(ped_cellcount, " ped ",cell_type,sep=""))
  
  
  samp_num<-10000
  pval<-0.05
  Tcell_pval_montecarlo<-sapply(1:samp_num, function(x){
    set.seed(x)
    plt_Tcell_ped_random<-plt_Tcell_ped[sample(ped_cellcount, adult_cellcount),]
    t.test(plt_Tcell_ped_random$exhausted_tcells_score1, plt_Tcell_adult$exhausted_tcells_score1)$p.value})
  print(paste("Comparing t cell exhaustion in ", cell_type, 
              ". Signature score is sig different between ped and adult (p <", pval, 
              ") in ", samp_num, " random samples at a p value of ", 
              round((length(Tcell_pval_montecarlo[which(Tcell_pval_montecarlo>pval)])+1)/(samp_num+1), 3), sep=""))
  (length(Tcell_pval_montecarlo[which(Tcell_pval_montecarlo>pval)])+1)/(samp_num+1)
})

                                  #   
                                  # #####
                                  # # Plot Signatures
                                  # #####
                                  # ####
                                  # ## myeloid_immune_supressive_score1
                                  # ####
                                  # ## UMAPs
                                  # immune_all<-ggplot(plt, aes(UMAP_1,UMAP_2, color=myeloid_immune_supressive_score1))+
                                  #   geom_point(size=1.5)+
                                  #   theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
                                  #   annotate("text", x = -11, y = -12, label = paste0("n = ",comma(nrow(plt))))+
                                  #   scale_colour_continuous("Immune\nSupressive\nSignature Score")
                                  # save_plts(immune_all, "immunesupressive_umap_all", w=7,h=5)
                                  # 
                                  # immune_all_agesplit<-ggplot(plt, aes(UMAP_1,UMAP_2))+
                                  #   geom_point(aes(color=myeloid_immune_supressive_score1),size=1.5)+
                                  #   theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
                                  #   facet_wrap(~AgeGroup)+
                                  #   geom_text(aes(x = -11, y = -12, label=paste0("n = ",comma(CellCount))), cell_num_all)+
                                  #   scale_colour_continuous("Immune\nSupressive\nSignature Score")
                                  # save_plts(immune_all_agesplit, "immunesupressive_umap_all_agesplit", w=12,h=5)
                                  # 
                                  # immune_myeloid<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2, color=myeloid_immune_supressive_score1))+
                                  #   geom_point(size=1.5)+
                                  #   theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
                                  #   annotate("text", x = -4, y = -12, label = paste0("n = ",comma(nrow(plt_myeloid))))+
                                  #   scale_colour_continuous("Immune\nSupressive\nSignature Score")
                                  # save_plts(immune_myeloid, "immunesupressive_umap_myeloid", w=7,h=5)
                                  # 
                                  # immune_myeloid_agesplit<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
                                  #   geom_point(aes(color=myeloid_immune_supressive_score1), size=1.5)+
                                  #   theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
                                  #   facet_wrap(~AgeGroup)+
                                  #   geom_text(aes(x = -4, y = -12, label=paste0("n = ",comma(CellCount))), cell_num_myeloid)+
                                  #   scale_colour_continuous("Immune\nSupressive\nSignature Score")
                                  # save_plts(immune_myeloid_agesplit, "immunesupressive_umap_myeloid_agesplit", w=12,h=5)
                                  # 
                                  # # BOX PLOT
                                  # plt_max<-ceiling(max(plt_myeloid$myeloid_immune_supressive_score1))
                                  # plt_min<-floor(min(plt_myeloid$myeloid_immune_supressive_score1))
                                  # 
                                  # myeloid_immune_box<-
                                  #   ggplot(plt_myeloid, aes(AgeGroup,myeloid_immune_supressive_score1))+
                                  #   geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=AgeGroup))+
                                  #   theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
                                  #   facet_grid(.~CellType_rough)+fillscale_age+
                                  #   geom_signif(stat="identity",
                                  #               data=data.frame(x=c(1, 1, 2), xend=c(1, 2, 2),
                                  #                               y=c(plt_max, plt_max+0.25, plt_max+0.25), yend=c(plt_max+0.25, plt_max+0.25, plt_max),
                                  #                               annotation=c("*")),
                                  #               aes(x=x,xend=xend, y=y, yend=yend, annotation=annotation), color="grey50")+ylim(plt_min, plt_max+1)+
                                  #   xlab("Age Group")+ylab("Immune Supressive Signature Score")+
                                  #   geom_text(aes(y=-1, x=AgeGroup,label=paste0("n = ",comma(CellCount))),cell_num_myeloid, hjust=-0.1, size=3)
                                  # save_plts(myeloid_immune_box, "immunesupressive_box_myeloid", w=4,h=4)
                                  # 
                                  # 
                                  # ####
                                  # ## inflammatory_macs_score1
                                  # ####
                                  # 
                                  # ## UMAPs
                                  # inflam_all<-ggplot(plt, aes(UMAP_1,UMAP_2, color=inflammatory_macs_score1))+
                                  #   geom_point(size=1.5)+
                                  #   theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
                                  #   annotate("text", x = -11, y = -12, label = paste0("n = ",comma(nrow(plt))))+
                                  #   scale_colour_continuous("Inflammatory\nMacrophage\nSignature Score")
                                  # save_plts(inflam_all, "inflammatory_umap_all", w=7,h=5)
                                  # 
                                  # inflam_all_agesplit<-ggplot(plt, aes(UMAP_1,UMAP_2))+
                                  #   geom_point(aes(color=inflammatory_macs_score1),size=1.5)+
                                  #   theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
                                  #   facet_wrap(~AgeGroup)+
                                  #   geom_text(aes(x = -11, y = -12, label=paste0("n = ",comma(CellCount))), cell_num_all)+
                                  #   scale_colour_continuous("Inflammatory\nMacrophage\nSignature Score")
                                  # save_plts(inflam_all_agesplit, "inflammatory_umap_all_agesplit", w=12,h=5)
                                  # 
                                  # inflam_myeloid<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2, color=inflammatory_macs_score1))+
                                  #   geom_point(size=1.5)+
                                  #   theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
                                  #   annotate("text", x = -4, y = -12, label = paste0("n = ",comma(nrow(plt_myeloid))))+
                                  #   scale_colour_continuous("Inflammatory\nMacrophage\nSignature Score")
                                  # save_plts(inflam_myeloid, "inflammatory_umap_myeloid", w=7,h=5)
                                  # 
                                  # inflam_myeloid_agesplit<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
                                  #   geom_point(aes(color=inflammatory_macs_score1), size=1.5)+
                                  #   theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
                                  #   facet_wrap(~AgeGroup)+
                                  #   geom_text(aes(x = -4, y = -12, label=paste0("n = ",comma(CellCount))), cell_num_myeloid)+
                                  #   scale_colour_continuous("Inflammatory\nMacrophage\nSignature Score")
                                  # save_plts(inflam_myeloid_agesplit, "inflammatory_umap_myeloid_agesplit", w=12,h=5)
                                  # 
                                  # 
                                  # 
                                  # # BOX PLOT
                                  # plt_max<-ceiling(max(plt_myeloid$inflammatory_macs_score1))
                                  # plt_min<-floor(min(plt_myeloid$inflammatory_macs_score1))
                                  # 
                                  # inflam_box_myeloid<-
                                  #   ggplot(plt_myeloid, aes(AgeGroup,inflammatory_macs_score1))+
                                  #   geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=AgeGroup))+
                                  #   theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
                                  #   facet_grid(.~CellType_rough)+fillscale_age+
                                  #   geom_signif(stat="identity",
                                  #               data=data.frame(x=c(1, 1, 2), xend=c(1, 2, 2),
                                  #                               y=c(plt_max, plt_max+0.25, plt_max+0.25), yend=c(plt_max+0.25, plt_max+0.25, plt_max),
                                  #                               annotation=c("*")),
                                  #               aes(x=x,xend=xend, y=y, yend=yend, annotation=annotation), color="grey50")+ylim(plt_min, plt_max+1)+
                                  #   xlab("Age Group")+ylab("Inflammatory Macrophage Signature Score")+
                                  #   geom_text(aes(y=-1.8, x=AgeGroup,label=paste0("n = ",comma(CellCount))),cell_num_myeloid, hjust=-0.1, size=3)
                                  # save_plts(inflam_box_myeloid, "inflammatory_box_myeloid", w=4,h=4)

####
## exhausted_tcells_score1
####

## UMAPs
tcell_exhaust_all<-ggplot(plt, aes(UMAP_1,UMAP_2, color=exhausted_tcells_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  annotate("text", x = -11, y = -12, label = paste0("n = ",comma(nrow(plt))))+
  scale_colour_continuous("Exhausted\nT Cell\nSignature Score")
save_plts(tcell_exhaust_all, "tcell_exhaustumap_all", w=7,h=5)

tcell_exhaust_all_agesplit<-ggplot(plt, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=exhausted_tcells_score1),size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~AgeGroup)+
  geom_text(aes(x = -11, y = -12, label=paste0("n = ",comma(CellCount))), cell_num_all)+
  scale_colour_continuous("Exhausted\nT Cell\nSignature Score")
save_plts(tcell_exhaust_all_agesplit, "tcell_exhaustumap_all_agesplit", w=12,h=5)

## T cell types for visualization
                    # plt_Tcell$CellType_rough<-as.factor(plt_Tcell$CellType_rough)
                    # levels(plt_Tcell$CellType_rough)<-c("NK_T")
                    # 
                    # tcell_color<-ggplot(plt_Tcell, aes(UMAP_1,UMAP_2, color=CellType_rough))+
                    #   geom_point(size=1.5)+
                    #   theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
                    #   annotate("text", x = -9, y = -6, label = paste0("n = ",comma(nrow(plt_Tcell))))+colscale_cellType
                    # save_plts(tcell_color, "tcell_umap", w=7,h=5)


tcell_exhaust_tcells<-ggplot(plt_Tcell, aes(UMAP_1,UMAP_2, color=exhausted_tcells_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  annotate("text", x = -9, y = -6, label = paste0("n = ",comma(nrow(plt_Tcell))))+
  scale_colour_continuous("Exhausted\nT Cell\nSignature Score")
save_plts(tcell_exhaust_tcells, "tcell_exhaustumap_tcell", w=7,h=5)


tcell_exhaust_tcells_agesplit<-ggplot(plt_Tcell, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=exhausted_tcells_score1), size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~AgeGroup)+
  geom_text(aes( x = -9, y = -6, label=paste0("n = ",comma(CellCount))), cell_num_tcell)+
  scale_colour_continuous("Exhausted\nT Cell\nSignature Score")
save_plts(tcell_exhaust_tcells_agesplit, "tcell_exhaustumap_tcell_agesplit", w=12,h=5)

# BOX PLOT
plt_max<-ceiling(max(plt_Tcell$exhausted_tcells_score1))
plt_min<-floor(min(plt_Tcell$exhausted_tcells_score1))

SignificanceDF<-data.frame(x=c(1, 1, 2), xend=c(1, 2, 2),
           y=c(plt_max, plt_max+0.1, plt_max+0.1), yend=c(plt_max+0.1, plt_max+0.1, plt_max),
           annotation=c("*"),
           CellType_rough="NK_T")

tcellexhaust_box_tcell<-ggplot(plt_Tcell, aes(AgeGroup,exhausted_tcells_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=AgeGroup))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~CellType_rough)+fillscale_age+
  geom_signif(stat="identity",
              data=SignificanceDF,
              aes(x=x,xend=xend, y=y, yend=yend, annotation=annotation), color="grey50")+ylim(plt_min+0.5, plt_max+0.25)+
  geom_text(aes(y=-0.5, x=AgeGroup,label=paste0("n = ",comma(CellCount))),cell_num_tcell, hjust=-0.05, size=3)+
  xlab("Age Group")+ylab("Exhausted T Cell Signature Score")
save_plts(tcellexhaust_box_tcell, "tcell_exhaust_box_tcell", w=8,h=4)
  
  

####
## recently_recruited_myeloid1
####
## UMAPs
recruit_all<-ggplot(plt, aes(UMAP_1,UMAP_2, color=recently_recruited_myeloid1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  annotate("text", x = -11, y = -12, label = paste0("n = ",comma(nrow(plt))))+
  scale_color_continuous_sequential(palette = "Mako") + 
  guides(color=guide_legend(title="Recently\nRecruited\nMyeloid\nSignature Score"))
save_plts(recruit_all, "recruit_umap_all", w=7,h=5)

recruit_all_agesplit<-ggplot(plt, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=recently_recruited_myeloid1),size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~AgeGroup)+
  geom_text(aes(x = -11, y = -12, label=paste0("n = ",comma(CellCount))), cell_num_all)+
  scale_color_continuous_sequential(palette = "Mako") + 
  guides(color=guide_legend(title="Recently\nRecruited\nMyeloid\nSignature Score"))
save_plts(recruit_all_agesplit, "recruit_umap_all_agesplit", w=12,h=5)

recruit_all_myeloid<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2, color=recently_recruited_myeloid1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  annotate("text", x = -5, y = -12, label = paste0("n = ",comma(nrow(plt_myeloid))))+
  scale_color_continuous_sequential(palette = "Mako") + 
  guides(color=guide_legend(title="Recently\nRecruited\nMyeloid\nSignature Score"))
save_plts(recruit_all_myeloid, "recruit_umap_myeloid", w=7,h=5)

recruit_all_myeloid_agesplit<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=recently_recruited_myeloid1), size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~AgeGroup)+
  geom_text(aes(x = -5, y = -12, label=paste0("n = ",comma(CellCount))), cell_num_myeloid)+
  scale_color_continuous_sequential(palette = "Mako") + 
  guides(color=guide_legend(title="Recently\nRecruited\nMyeloid\nSignature Score"))
save_plts(recruit_all_myeloid_agesplit, "recruit_umap_myeloid_agesplit", w=12,h=5)

# BOX PLOT
plt_max<-ceiling(max(plt_myeloid$recently_recruited_myeloid1))
plt_min<-floor(min(plt_myeloid$recently_recruited_myeloid1))+0.5

myeloid_recruit_box<-
  ggplot(plt_myeloid, aes(AgeGroup,recently_recruited_myeloid1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=AgeGroup))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~CellType_rough)+fillscale_age+
  geom_signif(stat="identity",
              data=data.frame(x=c(1, 1, 2), xend=c(1, 2, 2),
                              y=c(plt_max, plt_max+0.25, plt_max+0.25), yend=c(plt_max+0.25, plt_max+0.25, plt_max),
                              annotation=c("*")),
              aes(x=x,xend=xend, y=y, yend=yend, annotation=annotation), color="grey50")+ylim(plt_min, plt_max+1)+
  xlab("Age Group")+ylab("Recently Recruited Myeloid Signature Score")+
  geom_text(aes(y=-1, x=AgeGroup,label=paste0("n = ",comma(CellCount))),cell_num_myeloid, hjust=-0.1, size=3)
save_plts(myeloid_recruit_box, "recruit_box_myeloid", w=4,h=4)


####
## kuffer_like_score1
####
## UMAPs
kuffer_all<-ggplot(plt, aes(UMAP_1,UMAP_2, color=kuffer_like_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  annotate("text", x = -11, y = -12, label = paste0("n = ",comma(nrow(plt))))+
  scale_color_continuous_sequential(palette = "Mako") +
  guides(color=guide_legend(title="Kuffer-like\nSignature Score"))
save_plts(kuffer_all, "kuffer_umap_all", w=7,h=5)

kuffer_all_agesplit<-ggplot(plt, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=kuffer_like_score1),size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~AgeGroup)+
  geom_text(aes(x = -11, y = -12, label=paste0("n = ",comma(CellCount))), cell_num_all)+
  scale_color_continuous_sequential(palette = "Mako") +
  guides(color=guide_legend(title="Kuffer-like\nSignature Score"))
save_plts(kuffer_all_agesplit, "kuffer_umap_all_agesplit", w=12,h=5)

kuffer_all_myeloid<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2, color=kuffer_like_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  annotate("text", x = -5, y = -12, label = paste0("n = ",comma(nrow(plt_myeloid))))+
  scale_color_continuous_sequential(palette = "Mako") +
  guides(color=guide_legend(title="Kuffer-like\nSignature Score"))
save_plts(kuffer_all_myeloid, "kuffer_umap_myeloid", w=7,h=5)

kuffer_all_myeloid_agesplit<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=kuffer_like_score1), size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~AgeGroup)+
  geom_text(aes(x = -5, y = -12, label=paste0("n = ",comma(CellCount))), cell_num_myeloid)+
  scale_color_continuous_sequential(palette = "Mako") +
  guides(color=guide_legend(title="Kuffer-like\nSignature Score"))
save_plts(kuffer_all_myeloid_agesplit, "kuffer_umap_myeloid_agesplit", w=12,h=5)

# BOX PLOT
plt_max<-ceiling(max(plt_myeloid$kuffer_like_score1))
plt_min<-floor(min(plt_myeloid$kuffer_like_score1))+0.5

myeloid_kuffer_box<-
  ggplot(plt_myeloid, aes(AgeGroup,kuffer_like_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=AgeGroup))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~CellType_rough)+fillscale_age+
  geom_signif(stat="identity",
              data=data.frame(x=c(1, 1, 2), xend=c(1, 2, 2),
                              y=c(plt_max, plt_max+0.25, plt_max+0.25), yend=c(plt_max+0.25, plt_max+0.25, plt_max),
                              annotation=c("*")),
              aes(x=x,xend=xend, y=y, yend=yend, annotation=annotation), color="grey50")+ylim(plt_min, plt_max+1)+
  xlab("Age Group")+ylab("Kuffer-like Signature Score")+
  geom_text(aes(y=-1.5, x=AgeGroup,label=paste0("n = ",comma(CellCount))),cell_num_myeloid, hjust=-0.1, size=3)
save_plts(myeloid_kuffer_box, "kuffer_box_myeloid", w=4,h=4)



##############
## Age and Score associations
##############
plt_myeloid$CellType_refined<-as.character(plt_myeloid$CellType_refined)
cell_num_myeloid_type<-as.data.frame(table(plt_myeloid$AgeGroup, plt_myeloid$CellType_refined))
colnames(cell_num_myeloid_type)<-c("AgeGroup","CellType_refined","CellCount")

comp_simple<-list(c("Adult","Ped"))

ggplot(plt_myeloid, aes(AgeGroup,kuffer_like_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=AgeGroup))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  fillscale_age+
  ylim(plt_min, plt_max+1)+
  xlab("Age Group")+ylab("Kuffer-like Signature Score")+facet_wrap(~CellType_refined)+
  geom_text(aes(y=-1.5, x=AgeGroup,label=paste0("n = ",comma(CellCount))),cell_num_myeloid_type, hjust=-0.1, size=3)+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.5, vjust = 0.5,
              textsize = 5,  map_signif_level = T, color="grey60")


plt_KCRR<-plt_myeloid[which(plt_myeloid$CellType_refined %in% c("KC Like","RR Myeloid")),]
plt_KCRR$CellType_refined<-as.character(plt_KCRR$CellType_refined)
cell_num_KCRR_type<-as.data.frame(table(as.factor(plt_KCRR$Age)))
colnames(cell_num_KCRR_type)<-c("Age","CellCount")

ggplot(plt_KCRR, aes(as.factor(Age),kuffer_like_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.25,aes(fill=AgeGroup))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  fillscale_age+
  xlab("Age Group")+ylab("Kuffer-like Signature Score")+
  geom_text(aes(y=-1.5, x=as.factor(Age),label=comma(CellCount)),cell_num_KCRR_type, size=3)


plt_KCRR<-plt_myeloid[which(plt_myeloid$CellType_refined %in% c("KC Like","RR Myeloid")),]
plt_KCRR$CellType_refined<-as.character(plt_KCRR$CellType_refined)
cell_num_KCRR_type<-as.data.frame(table(as.factor(plt_KCRR$Age), plt_KCRR$CellType_refined))
colnames(cell_num_KCRR_type)<-c("Age","CellType_refined","CellCount")

ggplot(plt_KCRR, aes(as.factor(Age),kuffer_like_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.25,aes(fill=AgeGroup))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  fillscale_age+
  xlab("Age Group")+ylab("Kuffer-like Signature Score")+facet_wrap(.~CellType_refined, ncol=1)+
  geom_text(aes(y=-1.5, x=as.factor(Age),label=comma(CellCount)),cell_num_KCRR_type, size=3)




### RR myeloid
plt_myeloid$CellType_refined<-as.character(plt_myeloid$CellType_refined)
cell_num_myeloid_type<-as.data.frame(table(plt_myeloid$AgeGroup, plt_myeloid$CellType_refined))
colnames(cell_num_myeloid_type)<-c("AgeGroup","CellType_refined","CellCount")

comp_simple<-list(c("Adult","Ped"))

ggplot(plt_myeloid, aes(AgeGroup,recently_recruited_myeloid1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=AgeGroup))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  fillscale_age+
  ylim(plt_min, plt_max+1)+
  xlab("Age Group")+ylab("Recently Recruited Myeloid Signature Score")+facet_wrap(~CellType_refined)+
  geom_text(aes(y=-1.5, x=AgeGroup,label=paste0("n = ",comma(CellCount))),cell_num_myeloid_type, hjust=-0.1, size=3)+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.5, vjust = 0.5,
              textsize = 5,  map_signif_level = T, color="grey60")

plt_myeloid$CellType_refined<-as.character(plt_myeloid$CellType_refined)
cell_num_myeloid_type<-as.data.frame(table(as.factor(plt_myeloid$Age)))
colnames(cell_num_myeloid_type)<-c("Age","CellCount")

ggplot(plt_myeloid, aes(as.factor(Age),recently_recruited_myeloid1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=AgeGroup))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  fillscale_age+
  ylim(plt_min, plt_max+1)+
  xlab("Age Group")+ylab("Recently Recruited Myeloid Signature Score")+
  geom_text(aes(y=-1.5, x=as.factor(Age),label=paste0("n = ",comma(CellCount))),cell_num_myeloid_type, hjust=-0.1, size=3)+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.5, vjust = 0.5,
              textsize = 5,  map_signif_level = T, color="grey60")


plt_KCRR<-plt_myeloid[which(plt_myeloid$CellType_refined %in% c("KC Like","RR Myeloid")),]
plt_KCRR$CellType_refined<-as.character(plt_KCRR$CellType_refined)
cell_num_KCRR_type<-as.data.frame(table(as.factor(plt_KCRR$Age)))
colnames(cell_num_KCRR_type)<-c("Age","CellCount")

ggplot(plt_KCRR, aes(as.factor(Age),recently_recruited_myeloid1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.25,aes(fill=AgeGroup))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  fillscale_age+
  xlab("Age Group")+ylab("Recently Recruited Myeloid Signature Score")+
  geom_text(aes(y=-1.5, x=as.factor(Age),label=comma(CellCount)),cell_num_KCRR_type, size=3)


plt_KCRR<-plt_myeloid[which(plt_myeloid$CellType_refined %in% c("KC Like","RR Myeloid")),]
plt_KCRR$CellType_refined<-as.character(plt_KCRR$CellType_refined)
cell_num_KCRR_type<-as.data.frame(table(as.factor(plt_KCRR$Age), plt_KCRR$CellType_refined))
colnames(cell_num_KCRR_type)<-c("Age","CellType_refined","CellCount")

ggplot(plt_KCRR, aes(as.factor(Age),recently_recruited_myeloid1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.25,aes(fill=AgeGroup))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  fillscale_age+
  xlab("Age Group")+ylab("Recently Recruited Myeloid Signature Score")+facet_wrap(.~CellType_refined, ncol=1)+
  geom_text(aes(y=-1.5, x=as.factor(Age),label=comma(CellCount)),cell_num_KCRR_type, size=3)


ggplot(plt_KCRR, aes(Age,recently_recruited_myeloid1))+geom_point()+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  fillscale_age+stat_smooth(method="lm")+
  xlab("Age Group")+ylab("Recently Recruited Myeloid Signature Score")

ggplot(plt_KCRR, aes(Age,kuffer_like_score1))+geom_point()+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  fillscale_age+stat_smooth(method="lm")+
  xlab("Age Group")+ylab("Kuffer-like Signature Score")
