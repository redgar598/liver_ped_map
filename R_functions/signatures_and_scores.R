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



options(stringsAsFactors = FALSE)

source("R_functions/pretty_plots.R")



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

##########
## Load DE stats
##########
load(file=here("data","adult_ped_diff_genes.RData"))

#########
## Signatures
#########
myeloid_immune_supressive<-c("CTSB","CD163","MS4A7","FOLR2","GPNMB","VSIG4","HMOX1","MSR1")
inflammatory_macs<-c("CD74","HLA-DRA","TYROBP","C1QC","HLA-DPA1","HLA-DPB1","LYZ","S100A6")
exhausted_tcells<-c("TOX","PDCD1","LAG3","TNFRSF9","CXCL13","ENTPD1","HAVCR2","CD38")

diff_exp_all[which(diff_exp_all$gene%in%myeloid_immune_supressive & diff_exp_all$cell.1=="Myeloid_Adult"),]
diff_exp_all[which(diff_exp_all$gene%in%inflammatory_macs & diff_exp_all$cell.1=="Myeloid_Adult"),]
diff_exp_all[which(diff_exp_all$gene%in%exhausted_tcells & diff_exp_all$cell.1%in%c("CD3_Tcell","nkTcell","gdTcell")),]
#diff_exp_all[which(diff_exp_all$gene%in%exhausted_tcells & diff_exp_all$cell.1=="NK_T_B_Adult"),]

######
## plot genes
######
plot_key_genes<-function(keygenes, label){
  all_plots<-lapply(1:length(keygenes), function(y){
    plots<-lapply(1:length(cell_types),function(x){
      p<-VlnPlot(subset(d10x, subset = CellType_rough == cell_types[x]) , features = keygenes[y], pt.size = 0, log=T)
      p<-if(length(grep(cell_types[x], diff_exp_all[which(diff_exp_all$gene==keygenes[y]),]$cell.1))!=0){
        p+theme(plot.background = element_rect(color = "black",size = 2)) +fillscale_age +xlab("") + ylab("")+ theme(legend.position="none")}else{  
          p+fillscale_age +xlab("") + ylab("")+ theme(legend.position="none")
        }
      p})
    plot_grid(plotlist = plots, ncol=1)})
  
  
  label_blank<-lapply(1:length(cell_types), function(x){
    ggplot()+geom_blank()+theme_void()+ggtitle(cell_types[x])+ theme(plot.title = element_text(hjust = 0.5,vjust = -30))  })
  label_blank<-plot_grid(plotlist = label_blank, ncol=1)
  
  plot_grid(label_blank, plot_grid(plotlist=all_plots, ncol=length(keygenes)), rel_widths=c(0.1,1))
  ggsave2(paste0(here("figures/"), label, "_adult_ped.pdf"), w=20,h=20)
  ggsave2(paste0(here("figures/jpeg/"),label, "_adult_ped.jpeg"), w=20,h=20,bg="white")}

plot_key_genes(myeloid_immune_supressive, "myeloid_immune_supressive")
plot_key_genes(inflammatory_macs, "inflammatory_macs")
plot_key_genes(exhausted_tcells, "exhausted_tcells")




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

score_data<-d10x@meta.data[,c("myeloid_immune_supressive_score1","inflammatory_macs_score1","exhausted_tcells_score1")]


### load integrate for UMAP etc
load(here("data","adult_ped_integrated.rds"))
score_data<-score_data[match(rownames(d10x.combined@meta.data),rownames(score_data)),]
identical(rownames(d10x.combined@meta.data),rownames(score_data))
d10x.combined <- AddMetaData(d10x.combined, metadata = score_data)

                # d10x.combined@meta.data$CellType_rough[which(d10x.combined@meta.data$CellType_rough=="macrophage")]<-"Myeloid"
                # d10x.combined@meta.data$CellType_rough[which(d10x.combined@meta.data$seurat_clusters==5)]<-"CD3_Tcell"
                # d10x.combined@meta.data$CellType_rough[which(d10x.combined@meta.data$seurat_clusters==11)]<-"nkTcell"
                # d10x.combined@meta.data$CellType_rough[which(d10x.combined@meta.data$seurat_clusters==16)]<-"gdTcell"


######
## plot scores
######
umap_mat<-as.data.frame(Embeddings(object = d10x.combined, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x.combined@meta.data
meta$cell<-rownames(meta)


plt<-merge(meta, umap_mat, by="cell")


d10x.combined_NK_T_B<-subset(d10x.combined, subset = CellType_rough %in% c("CD3_Tcell","nkTcell","gdTcell"))
d10x.combined_NK_T_B <- RunPCA(d10x.combined_NK_T_B, npcs = 30, verbose = FALSE)
d10x.combined_NK_T_B <- RunUMAP(d10x.combined_NK_T_B, reduction = "pca", dims = 1:30)
umap_mat_Tcell<-as.data.frame(Embeddings(object = d10x.combined_NK_T_B, reduction = "umap"))#
umap_mat_Tcell$cell<-rownames(umap_mat_Tcell)
meta_Tcell<-d10x.combined_NK_T_B@meta.data
meta_Tcell$cell<-rownames(meta_Tcell)
plt_Tcell<-merge(meta_Tcell, umap_mat_Tcell, by="cell")
cell_num_tcell<-as.data.frame(table(plt_Tcell$AgeGroup))
colnames(cell_num_tcell)<-c("AgeGroup","CellCount")

d10x.combined_myeloid<-subset(d10x.combined, subset = CellType_rough %in% c("Myeloid"))
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


#######
### statistics
#######
#############myeloid_immune_supressive_score1
immune_all<-ggplot(plt, aes(UMAP_1,UMAP_2, color=myeloid_immune_supressive_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  annotate("text", x = -10, y = -12, label = paste0("n = ",comma(nrow(plt))))+
  scale_colour_continuous("Immune\nSupressive\nSignature Score")
save_plts(immune_all, "immunesupressive_umap_all", w=7,h=5)


immune_all_agesplit<-ggplot(plt, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=myeloid_immune_supressive_score1),size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~AgeGroup)+
  geom_text(aes(x=-10, y=-12, label=paste0("n = ",comma(CellCount))), cell_num_all)+
  scale_colour_continuous("Immune\nSupressive\nSignature Score")
save_plts(immune_all_agesplit, "immunesupressive_umap_all_agesplit", w=12,h=5)


immune_myeloid<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2, color=myeloid_immune_supressive_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  annotate("text", x = -6, y = -12, label = paste0("n = ",comma(nrow(plt_myeloid))))+
  scale_colour_continuous("Immune\nSupressive\nSignature Score")
save_plts(immune_myeloid, "immunesupressive_umap_myeloid", w=7,h=5)


immune_myeloid_agesplit<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=myeloid_immune_supressive_score1), size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~AgeGroup)+
  geom_text(aes(x=-6, y=-12, label=paste0("n = ",comma(CellCount))), cell_num_myeloid)+
  scale_colour_continuous("Immune\nSupressive\nSignature Score")
save_plts(immune_myeloid_agesplit, "immunesupressive_umap_myeloid_agesplit", w=12,h=5)

# Group the data by cell type do another anova
stats_plt_age<-lapply(cell_types, function(x){
  print("##################")
  cell_meta<-meta[which(meta$CellType_rough==x),]
  print(unique(cell_meta$CellType_rough))
  print(summary(aov(myeloid_immune_supressive_score1 ~ AgeGroup, data = cell_meta)))
  pairwise<-pairwise.t.test(cell_meta$myeloid_immune_supressive_score1, cell_meta$AgeGroup,p.adjust.method = "BH", pool.sd = FALSE)
  print(pairwise)
  pairwise<-melt(pairwise$p.value)
  pairwise[!is.na(pairwise$value),]
  pairwise$AgeGroup<-x
  pairwise
})

stats_plt_age<-do.call(rbind, stats_plt_age)
stats_plt_age$Var1<-as.character(stats_plt_age$Var1)
stats_plt_age$Var2<-as.character(stats_plt_age$Var2)

comp_simple<-list(c("Ped","Adult"))
### CONFIRM P vale match the pairwise above
all_immune_box<-ggplot(plt, aes(AgeGroup,myeloid_immune_supressive_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=AgeGroup))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~CellType_rough)+fillscale_age+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")#
save_plts(all_immune_box, "immunesupressive_box_all_cell", w=12,h=3)


plt_myeloid<-plt[which(plt$CellType_rough=="Myeloid"),]
myeloid_immune_box<-ggplot(plt_myeloid, aes(AgeGroup,myeloid_immune_supressive_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=AgeGroup))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~CellType_rough)+fillscale_age+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")+
  geom_text(aes(y=-1, x=AgeGroup,label=paste0("n = ",comma(CellCount))),cell_num_myeloid, hjust=-0.1, size=3)+
  xlab("Age Group")+ylab("Immune Supressive Signature Score")
save_plts(myeloid_immune_box, "immunesupressive_box_myeloid", w=4,h=4)

  


############ inflammatory_macs_score1
inflam_all<-ggplot(plt, aes(UMAP_1,UMAP_2, color=inflammatory_macs_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  annotate("text", x = -10, y = -12, label = paste0("n = ",comma(nrow(plt))))+
  scale_colour_continuous("Inflammatory\nMacrophage\nSignature Score")
save_plts(inflam_all, "inflammatory_umap_all", w=7,h=5)

inflam_all_agesplit<-ggplot(plt, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=inflammatory_macs_score1),size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~AgeGroup)+
  geom_text(aes(x=-10, y=-12, label=paste0("n = ",comma(CellCount))), cell_num_all)+
  scale_colour_continuous("Inflammatory\nMacrophage\nSignature Score")
save_plts(inflam_all_agesplit, "inflammatory_umap_all_agesplit", w=12,h=5)

inflam_myeloid<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2, color=inflammatory_macs_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  annotate("text", x = -6, y = -12, label = paste0("n = ",comma(nrow(plt_myeloid))))+
  scale_colour_continuous("Inflammatory\nMacrophage\nSignature Score")
save_plts(inflam_myeloid, "inflammatory_umap_myeloid", w=7,h=5)

inflam_myeloid_agesplit<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=inflammatory_macs_score1), size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~AgeGroup)+
  geom_text(aes(x=-6, y=-12, label=paste0("n = ",comma(CellCount))), cell_num_myeloid)+
  scale_colour_continuous("Inflammatory\nMacrophage\nSignature Score")
save_plts(inflam_myeloid_agesplit, "inflammatory_umap_myeloid_agesplit", w=12,h=5)


# Group the data by cell type do another anova
stats_plt_age<-lapply(cell_types, function(x){
  print("##################")
  cell_meta<-meta[which(meta$CellType_rough==x),]
  print(unique(cell_meta$CellType_rough))
  print(summary(aov(inflammatory_macs_score1 ~ AgeGroup, data = cell_meta)))
  pairwise<-pairwise.t.test(cell_meta$inflammatory_macs_score1, cell_meta$AgeGroup,p.adjust.method = "BH", pool.sd = FALSE)
  print(pairwise)
  pairwise<-melt(pairwise$p.value)
  pairwise[!is.na(pairwise$value),]
  pairwise$AgeGroup<-x
  pairwise
})

stats_plt_age<-do.call(rbind, stats_plt_age)
stats_plt_age$Var1<-as.character(stats_plt_age$Var1)
stats_plt_age$Var2<-as.character(stats_plt_age$Var2)

comp_simple<-list(c("Ped","Adult"))
### CONFIRM P vale match the pairwise above
all_inflam_box<-ggplot(plt, aes(AgeGroup,inflammatory_macs_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=AgeGroup))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~CellType_rough)+fillscale_age+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")#
save_plts(all_inflam_box, "inflammatory_box_all_cell", w=12,h=3)

inflam_box_myeloid<-ggplot(plt_myeloid, aes(AgeGroup,inflammatory_macs_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=AgeGroup))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~CellType_rough)+fillscale_age+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")+
  geom_text(aes(y=-1, x=AgeGroup,label=paste0("n = ",comma(CellCount))),cell_num_myeloid, hjust=-0.1, size=3)+
  xlab("Age Group")+ylab("Inflammatory Macrophage Signature Score")
save_plts(inflam_box_myeloid, "inflammatory__box_myeloid", w=4,h=4)



###############exhausted_tcells_score1
tcell_exhaust_all<-ggplot(plt, aes(UMAP_1,UMAP_2, color=exhausted_tcells_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  annotate("text", x = -10, y = -12, label = paste0("n = ",comma(nrow(plt))))+
  scale_colour_continuous("Exhausted\nT Cell\nSignature Score")
save_plts(tcell_exhaust_all, "tcell_exhaustumap_all", w=7,h=5)


tcell_exhaust_all_agesplit<-ggplot(plt, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=exhausted_tcells_score1),size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~AgeGroup)+
  geom_text(aes(x=-10, y=-12, label=paste0("n = ",comma(CellCount))), cell_num_all)+
  scale_colour_continuous("Exhausted\nT Cell\nSignature Score")
save_plts(tcell_exhaust_all_agesplit, "tcell_exhaustumap_all_agesplit", w=12,h=5)


tcell_exhaust_myeloid<-ggplot(plt_Tcell, aes(UMAP_1,UMAP_2, color=exhausted_tcells_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  annotate("text", x = -6, y = -12, label = paste0("n = ",comma(nrow(plt_Tcell))))+
  scale_colour_continuous("Exhausted\nT Cell\nSignature Score")
save_plts(tcell_exhaust_myeloid, "tcell_exhaustumap_tcell", w=7,h=5)


tcell_exhaust_myeloid_agesplit<-ggplot(plt_Tcell, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=exhausted_tcells_score1), size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~AgeGroup)+
  geom_text(aes(x=-6, y=-12, label=paste0("n = ",comma(CellCount))), cell_num_tcell)+
  scale_colour_continuous("Exhausted\nT Cell\nSignature Score")
save_plts(tcell_exhaust_myeloid_agesplit, "tcell_exhaustumap_tcell_agesplit", w=12,h=5)


# Group the data by cell type do another anova
stats_plt_age<-lapply(cell_types, function(x){
  print("##################")
  cell_meta<-meta[which(meta$CellType_rough==x),]
  print(unique(cell_meta$CellType_rough))
  print(summary(aov(exhausted_tcells_score1 ~ AgeGroup, data = cell_meta)))
  pairwise<-pairwise.t.test(cell_meta$exhausted_tcells_score1, cell_meta$AgeGroup,p.adjust.method = "BH", pool.sd = FALSE)
  print(pairwise)
  pairwise<-melt(pairwise$p.value)
  pairwise[!is.na(pairwise$value),]
  pairwise$AgeGroup<-x
  pairwise
})

stats_plt_age<-do.call(rbind, stats_plt_age)
stats_plt_age$Var1<-as.character(stats_plt_age$Var1)
stats_plt_age$Var2<-as.character(stats_plt_age$Var2)

comp_simple<-list(c("Ped","Adult"))
### CONFIRM P vale match the pairwise above
all_tcell_exhaust_box<-ggplot(plt, aes(AgeGroup,exhausted_tcells_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=AgeGroup))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~CellType_rough)+fillscale_age+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")+
  xlab("Age Group")+ylab("Exhausted T Cell Signature Score")
save_plts(all_tcell_exhaust_box, "tcell_exhaust_box_all_cell", w=12,h=3)


tcellexhaust_box_myeloid<-ggplot(plt_Tcell, aes(AgeGroup,inflammatory_macs_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=AgeGroup))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~CellType_rough)+fillscale_age+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")+
  geom_text(aes(y=-1, x=AgeGroup,label=paste0("n = ",comma(CellCount))),cell_num_myeloid, hjust=-0.1, size=3)+
  xlab("Age Group")+ylab("Exhausted T Cell Signature Score")
save_plts(tcellexhaust_box_myeloid, "tcell_exhaust_box_myeloid", w=8,h=4)
