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
library(ggsignif)


source("scripts/00_pretty_plots.R")
source("scripts/00_fanciest_UMAP.R")
source("scripts/00_plot_gene_exp.R")
source("scripts/00_entropy_d10x.R")




## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))


######
## add cell type labels
######
load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
d10x_myeloid<-subset(d10x, subset = CellType_refined %in% c("RR Myeloid","Macrophage\n(MHCII high)","KC Like","Macrophage\n(CLEC9A high)","Cycling Myeloid","Myeloid Erythrocytes\n(phagocytosis)"))




##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x_myeloid <- NormalizeData(d10x_myeloid,scale.factor = 10000, normalization.method = "LogNormalize")
Idents(d10x_myeloid) <- "age_condition"
cell_types<-unique(d10x_myeloid$CellType_refined)


#########
## Signatures
#########
recent_recruit_myeloid<-c("S100A8","S100A9","CD68","LYZ")
kuffer_signature<-c("VSIG4","MARCO","CD5L","HMOX1")


######
## plot genes
######
plot_key_genes<-function(keygenes, label){
  all_plots<-lapply(1:length(keygenes), function(y){
    plots<-lapply(1:length(cell_types),function(x){
      p<-VlnPlot(subset(d10x_myeloid, subset = CellType_refined == cell_types[x]) , features = keygenes[y], pt.size = 0, log=T)+
        fillscale_agecondition +xlab("") + ylab("")+ theme(legend.position="none")
      p})
    plot_grid(plotlist = plots, ncol=1)})
  
  label_blank<-lapply(1:length(cell_types), function(x){
    ggplot()+geom_blank()+theme_void()+ggtitle(cell_types[x])+ theme(plot.title = element_text(hjust = 0.5,vjust = -30))  })
  label_blank<-plot_grid(plotlist = label_blank, ncol=1)
  
  plot_grid(label_blank, plot_grid(plotlist=all_plots, ncol=length(keygenes)), rel_widths=c(0.1,1))
  ggsave2(paste0(here("figures/"), label, "_adult_ped.pdf"), w=20,h=20)
  ggsave2(paste0(here("figures/jpeg/"),label, "_adult_ped.jpeg"), w=20,h=20,bg="white")}

plot_key_genes(recent_recruit_myeloid, "recent_recruit_myeloid_IFALD")
plot_key_genes(kuffer_signature, "kuffer_signature_IFALD")




######
## Score Signatures
######
d10x_myeloid <- AddModuleScore(
  object = d10x_myeloid,
  features = list(recent_recruit_myeloid),
  ctrl = 5,
  name = 'recent_recruit_myeloid'
)

d10x_myeloid <- AddModuleScore(
  object = d10x_myeloid,
  features = list(kuffer_signature),
  ctrl = 5,
  name = 'kuffer_signature'
)


score_data<-d10x_myeloid@meta.data[,c("recent_recruit_myeloid1","kuffer_signature1")]


### load integrate for UMAP etc
load(here("data","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

d10x.combined_myeloid<-subset(d10x.combined, subset = CellType_refined %in% c("RR Myeloid","Macrophage\n(MHCII high)","KC Like","Macrophage\n(CLEC9A high)","Cycling Myeloid","Myeloid Erythrocytes\n(phagocytosis)"))
rm(d10x.combined)
gc()
d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)
d10x.combined_myeloid <- RunUMAP(d10x.combined_myeloid, reduction = "pca", dims = 1:30)

score_data<-score_data[match(rownames(d10x.combined_myeloid@meta.data),rownames(score_data)),]
identical(rownames(d10x.combined_myeloid@meta.data),rownames(score_data))
d10x.combined_myeloid <- AddMetaData(d10x.combined_myeloid, metadata = score_data)

######
## plot scores
######
umap_mat<-as.data.frame(Embeddings(object = d10x.combined_myeloid, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x.combined_myeloid@meta.data
meta$cell<-rownames(meta)

plt_myeloid<-merge(meta, umap_mat, by="cell")

cell_num_myeloid<-as.data.frame(table(plt_myeloid$age_condition))
colnames(cell_num_myeloid)<-c("age_condition","CellCount")





#######
### KC like signature
#######

KC_myeloid<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=kuffer_signature1), size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  annotate("text", x = -6, y = -12, label = paste0("n = ",comma(nrow(plt_myeloid))))+
  scale_color_continuous_sequential(palette = "Viridis", rev=F,name="KC Signature Score")
save_plts(KC_myeloid, "KC_Signature_umap_myeloid_IFALD", w=7,h=5)

KC_myeloid_agesplit<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=kuffer_signature1), size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~age_condition)+
  geom_text(aes(x=-6, y=-12, label=paste0("n = ",comma(CellCount))), cell_num_myeloid)+
  scale_color_continuous_sequential(palette = "Viridis", rev=F,name="KC Signature Score")

save_plts(KC_myeloid_agesplit, "KC_Signature_umap_myeloid_agesplit_IFALD", w=12,h=5)

print(summary(aov(kuffer_signature1 ~ age_condition, data = meta)))
pairwise<-pairwise.t.test(meta$kuffer_signature1, meta$age_condition,p.adjust.method = "BH", pool.sd = FALSE)
print(pairwise)

comp_simple<-list(c("Ped IFALD","Adult Healthy"),c("Ped Healthy","Adult Healthy"),c("Ped IFALD","Ped Healthy"))

KC_box<-ggplot(plt_myeloid, aes(age_condition,kuffer_signature1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=age_condition))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  fillscale_agecondition+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")#
save_plts(KC_box, "KC_Signature_box_all_cell_IFALD", w=5,h=4)


  

#######
### RR signature
#######

RR_myeloid<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2, color=recent_recruit_myeloid1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  annotate("text", x = -6, y = -12, label = paste0("n = ",comma(nrow(plt_myeloid))))+
  scale_color_continuous_sequential(palette = "Viridis", rev=F,name="RR Signature Score")
save_plts(RR_myeloid, "RR_Signature_umap_myeloid_IFALD", w=7,h=5)

RR_myeloid_agesplit<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=recent_recruit_myeloid1), size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~age_condition)+
  geom_text(aes(x=-6, y=-12, label=paste0("n = ",comma(CellCount))), cell_num_myeloid)+
  scale_color_continuous_sequential(palette = "Viridis", rev=F,name="RR Signature Score")

save_plts(RR_myeloid_agesplit, "RR_Signature_umap_myeloid_agesplit_IFALD", w=12,h=5)




print(summary(aov(recent_recruit_myeloid1 ~ age_condition, data = meta)))
pairwise<-pairwise.t.test(meta$recent_recruit_myeloid1, meta$age_condition,p.adjust.method = "BH", pool.sd = FALSE)
print(pairwise)

comp_simple<-list(c("Ped IFALD","Adult Healthy"),c("Ped Healthy","Adult Healthy"),c("Ped IFALD","Ped Healthy"))

RR_box<-ggplot(plt_myeloid, aes(age_condition,recent_recruit_myeloid1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=age_condition))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  fillscale_agecondition+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")#
save_plts(RR_box, "RR_Signature_box_all_cell_IFALD", w=5,h=4)


