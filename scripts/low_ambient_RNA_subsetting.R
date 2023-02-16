

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


####################
## Are there adults with less ambient RNA?
####################
load(here("data","adult_ped_integrated_refinedlabels_withDropletQC.rds"))

d10x_adult<-subset(d10x.combined, subset = AgeGroup %in% c("Adult"))
ALB<- FetchData(d10x_adult, vars = "ALB") 
ALB$cell<-rownames(ALB)



umap_mat<-as.data.frame(Embeddings(object = d10x.combined, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x.combined@meta.data
meta$cell<-rownames(meta)

plt<-merge(meta, umap_mat, by="cell")
plt<-merge(plt, ALB, by="cell")

plt<-plt[order(plt$cell_status),]

soup_change<-plot_grid(ggplot(plt, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = CellType_refined), size=0.5)+theme_bw()+colscale_cellType+ guides(colour = guide_legend(override.aes = list(size=2))) + facet_wrap(~individual, ncol=1),
                       ggplot(plt, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = relALBChange), size=0.5)+theme_bw()+scale_colour_gradientn(colours=c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704'),limits=c(NA,NA), name="Relative\nALB\nChange")+ facet_wrap(~individual, ncol=1),
                       ggplot(plt, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = rna_ALB), size=0.5)+theme_bw()+ scale_color_continuous_sequential(palette = "Mako")+ facet_wrap(~individual, ncol=1),
                       ncol=3, align="v")

soup_change
save_plts(soup_change, "adults_sample_soup_change", w=12,h=12)


#################
## low ambient adults and peds
#################
d10x_adultlow_ped<-subset(d10x.combined, subset = individual %in% c("C70_caud5pr","C64_caud5pr","C85_caud3pr","C93_caud3pr", "C96_caud3pr", "C82_caud3pr"))


cell_count<-as.data.frame(tapply(d10x_adultlow_ped$orig.ident, d10x_adultlow_ped$AgeGroup, length))
colnames(cell_count)<-"cell_count"
cell_count$AgeGroup<-rownames(cell_count)

lowadult_ped_refined_cluster_umap_nolab<-DimPlot(d10x_adultlow_ped, reduction = "umap", pt.size=0.25, group.by = "CellType_refined", split.by = "AgeGroup")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  geom_text(aes(x = -11, y = -12, label=paste0("n = ",comma(cell_count))), cell_count)
lowadult_ped_refined_cluster_umap_nolab
save_plts(lowadult_ped_refined_cluster_umap_nolab, "low_ambient_adult_plus_all_ped", w=9,h=5)


##################
## GSEA
##################
d10x<-readRDS(file = here("data","d10x_adult_ped_raw.rds"))
d10x_adultlow_ped<-subset(d10x, subset = individual %in% c("C70_caud5pr","C64_caud5pr","C85_caud3pr","C93_caud3pr", "C96_caud3pr", "C82_caud3pr"))


######
## add cell type labels
######
load(here("data","adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x_adultlow_ped), cell_label$index),]
identical(colnames(d10x_adultlow_ped), cell_label$index)

d10x_adultlow_ped <- AddMetaData(d10x_adultlow_ped, metadata = cell_label)




##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x_adultlow_ped <- NormalizeData(d10x_adultlow_ped,scale.factor = 10000, normalization.method = "LogNormalize")


## testing factor
levels(d10x_adultlow_ped$CellType_refined)[which(levels(d10x_adultlow_ped$CellType_refined)=="LSEC\n(Hepatocyte Like)")]<-"LSEC_hep"
levels(d10x_adultlow_ped$CellType_refined)[which(levels(d10x_adultlow_ped$CellType_refined)=="LSEC")]<-"LSEC_nothep"
levels(d10x_adultlow_ped$CellType_refined)[which(levels(d10x_adultlow_ped$CellType_refined)=="Neutrophil\n(DEFA+)")]<-"Neutrophil_DEFA"
levels(d10x_adultlow_ped$CellType_refined)[which(levels(d10x_adultlow_ped$CellType_refined)=="Neutrophil")]<-"Neutrophil_notDEFA"

d10x_adultlow_ped$cell_age<-paste(d10x_adultlow_ped$CellType_refined, d10x_adultlow_ped$AgeGroup, sep = "_")
Idents(d10x_adultlow_ped) <- "cell_age"

table(d10x_adultlow_ped$CellType_refined, d10x_adultlow_ped$AgeGroup)


#MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene,
#consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
#and linear regression for the continuous process (i.e., the expression level). 

cell_types<-unique(as.character(d10x_adultlow_ped$CellType_refined))
#no neutrophils in peds
cell_types<-cell_types[-grep("Neutrophil",cell_types)]
cell_types<-cell_types[-grep("Low Quality",cell_types)]
cell_types[grep("CD3",cell_types)]<-"CD3"

contrasts_celltype_age<-do.call(rbind,lapply(1:length(cell_types), function(x){
  combinations(n = 2, r = 2, v = d10x_adultlow_ped$cell_age[grep(cell_types[x],d10x_adultlow_ped$cell_age)], repeats.allowed = FALSE)}))

contrasts_celltype_age

nrow(contrasts_celltype_age)

tapply(d10x_adultlow_ped$cell, list(d10x_adultlow_ped$CellType_refined,d10x_adultlow_ped$AgeGroup), length)

### get fold change in whole cohort for pathway analysis
## run DE 

diff_exp_all<-lapply(1:nrow(contrasts_celltype_age), function(x){
  de<-FindMarkers(d10x_adultlow_ped, ident.1 = contrasts_celltype_age[x,1], ident.2 = contrasts_celltype_age[x,2], test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
  print(paste(contrasts_celltype_age[x,1],"vs", contrasts_celltype_age[x,2],":", nrow(de), sep=" "))
  de$gene<-rownames(de)
  rownames(de)<-NULL
  de<-de[,c(6,1:5)]
  de$cell.1<-contrasts_celltype_age[x,1]
  de$cell.2<-contrasts_celltype_age[x,2]
  de})


diff_exp_all_adultlow_subset<-do.call(rbind, diff_exp_all)
save(file=here("data","diff_exp_all_adultlow_subset.RData"),diff_exp_all_adultlow_subset)


# DE_monte_carlo_sig<-diff_exp_all_adultlow_subset[which(diff_exp_all_adultlow_subset$p_val_adj<0.001),]
# 
# table(DE_monte_carlo_sig$cell.1)
# table(DE_monte_carlo_sig$gene)[which(table(DE_monte_carlo_sig$gene)>0)]
# length(unique(DE_monte_carlo_sig$gene))
# sig_genes<-table(DE_monte_carlo_sig$gene)[which(table(DE_monte_carlo_sig$gene)>0)]
# ord_name<-names(sig_genes[rev(order(sig_genes))])
# DE_monte_carlo_sig$cell<-gsub("_Adult","",DE_monte_carlo_sig$cell.1)
# 
# summary_tbl<-as.data.frame(DE_monte_carlo_sig %>%
#                              select(gene, cell) %>% 
#                              group_by(gene) %>%
#                              mutate(all_cells = paste(cell, collapse = " | "))%>%
#                              select(gene, all_cells))
# summary_tbl<-summary_tbl[!duplicated(summary_tbl),]
# summary_tbl$all_cells<-gsub("\n"," ", summary_tbl$all_cells)
# 
# write.csv(file=here("data","Significant_genes_adult_ped.csv"),summary_tbl[match(ord_name,summary_tbl$gene),])
# 

#########
## pathway adult versus ped in KC and RR
#########
source("scripts/00_GSEA_function.R")
GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

pathway_plt<-function(de){
  gene_list = de$avg_log2FC
  names(gene_list) = de$gene
  gene_list = sort(gene_list, decreasing = TRUE)
  gene_list = gene_list[!duplicated(names(gene_list))]
  
  res = GSEA(gene_list, GO_file, pval = 0.05)
  
  plt_path<-res$Results
  if(nrow(plt_path)>0){
  plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
  plt_path$Enrichment_Cell<-"Up-regulated in \nAdult"
  plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \nPediatric"
  
  plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))
  
  plt_path$direction_label<-as.factor(plt_path$Enrichment)
  levels(plt_path$direction_label)<-c(0.1,-0.1)
  plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))
  
  # top and bottom 15
  if(nrow(plt_path)>30){plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])}
  
  ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
    theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
    geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
    geom_vline(xintercept = 0, color="grey40")+scale_fill_manual(values=c("#fd8d3c","#6baed6"))+ 
    guides(fill = guide_legend(override.aes = list(size=5)))}}

RR_GSEA<-pathway_plt(diff_exp_all_adultlow_subset[which(diff_exp_all_adultlow_subset$cell.1=="RR Myeloid_Adult"),])
save_plts(RR_GSEA, "GSEA_adult_ped_recently_recruited_lowambient", w=20,h=10)

KC_GSEA<-pathway_plt(diff_exp_all_adultlow_subset[which(diff_exp_all_adultlow_subset$cell.1=="KC Like_Adult"),])
save_plts(KC_GSEA, "GSEA_adult_ped_KClike_lowambient", w=20,h=10)

NK_GSEA<-pathway_plt(diff_exp_all_adultlow_subset[which(diff_exp_all_adultlow_subset$cell.1=="NK-like cells_Adult"),])
save_plts(NK_GSEA, "GSEA_adult_ped_NKlike_lowambient", w=20,h=10)

## only 3 adult cells!
# HSC_GSEA<-pathway_plt(diff_exp_all_adultlow_subset[which(diff_exp_all_adultlow_subset$cell.1=="HSC_Adult"),])
# save_plts(HSC_GSEA, "GSEA_adult_ped_HSClike_lowambient", w=20,h=10)



###################
## Correlation in FC in full cohort and subset
###################
load(here("data/diff_exp_all_adultlow_subset.RData"))
load(here("data/diff_exp_all_fullcohort.RData"))

diff_exp_all_celltype<-merge(diff_exp_all[grep("RR", diff_exp_all$cell.1),], diff_exp_all_adultlow_subset[grep("RR", diff_exp_all_adultlow_subset$cell.1),], by="gene",  suffixes = c("_fullcohort","_adultsubset"))

load(here("data",paste("RR Myeloid","adult_ped_diff_motecarlo_1000.RData",sep="_")))
DE_monte_carlo_sig<-DE_monte_carlo[which(DE_monte_carlo$monte_carlo_sig<0.001),]

sig_MC_RR<-DE_monte_carlo_sig[grep("RR",DE_monte_carlo_sig$cell),]

diff_exp_all_celltype$sig<-sapply(1:nrow(diff_exp_all_celltype), function(x){
 if(diff_exp_all_celltype$gene[x]%in%sig_MC_RR$gene){"Significant"}else{"NS"}})
diff_exp_all_celltype_label<-diff_exp_all_celltype[which(diff_exp_all_celltype$gene%in%c(sig_MC_RR$gene)),]

ggplot(diff_exp_all_celltype, aes(avg_log2FC_fullcohort, avg_log2FC_adultsubset, color=sig))+
  geom_point()+th_present+theme_bw()+
  xlab("Full Cohort Differential Expression\n(Fold change)")+
  ylab("Adult Subset Differential Expression\n(log2 Fold change)")+
  scale_color_manual(values = c("grey60","#b80783"))+
  geom_text(aes(label=gene),diff_exp_all_celltype_label,color="black",vjust=-0.75, hjust=1,size=2)+
  geom_vline(xintercept = c(-1,1), color="grey")+  geom_hline(yintercept = c(-1,1), color="grey")+
  ylim(c(min(diff_exp_all_celltype$avg_log2FC_fullcohort),max(diff_exp_all_celltype$avg_log2FC_fullcohort)))+
  xlim(c(min(diff_exp_all_celltype$avg_log2FC_fullcohort),max(diff_exp_all_celltype$avg_log2FC_fullcohort)))+
  annotate("text", x=2.7, y=4.5, label="Higher Expression in Adults")+
  annotate("text", x=-1, y=-2.5, label="Higher Expression in Peds")+
  geom_segment(aes(xend = avg_log2FC_fullcohort, yend = avg_log2FC_fullcohort),alpha=0.25) +
  geom_abline(slope=1, intercept=0, color="grey")
ggsave(file=here("figures/FC_correlation_RR_fullcohort_adultsubset_differential_Ped_adult.pdf"), h=7, w=8)
ggsave(file=here("figures/jpeg/FC_correlation_RR_fullcohort_adultsubset_differential_Ped_adult.jpeg"), h=7, w=8)


diff_exp_all_celltype<-merge(diff_exp_all[grep("KC", diff_exp_all$cell.1),], diff_exp_all_adultlow_subset[grep("KC", diff_exp_all_adultlow_subset$cell.1),], by="gene",  suffixes = c("_fullcohort","_adultsubset"))

load(here("data",paste("KC Like","adult_ped_diff_motecarlo_1000.RData",sep="_")))
DE_monte_carlo_sig<-DE_monte_carlo[which(DE_monte_carlo$monte_carlo_sig<0.001),]

sig_MC_KC<-DE_monte_carlo_sig[grep("KC",DE_monte_carlo_sig$cell),]

diff_exp_all_celltype$sig<-sapply(1:nrow(diff_exp_all_celltype), function(x){
  if(diff_exp_all_celltype$gene[x]%in%sig_MC_KC$gene){"Significant"}else{"NS"}})
diff_exp_all_celltype_label<-diff_exp_all_celltype[which(diff_exp_all_celltype$gene%in%c(sig_MC_KC$gene)),]

ggplot(diff_exp_all_celltype, aes(avg_log2FC_fullcohort, avg_log2FC_adultsubset, color=sig))+
  geom_point()+th_present+theme_bw()+
  xlab("Full Cohort Differential Expression\n(Fold change)")+
  ylab("Adult Subset Differential Expression\n(log2 Fold change)")+
  scale_color_manual(values = c("grey60","#b80783"))+
  geom_text(aes(label=gene),diff_exp_all_celltype_label,color="black",vjust=-0.75, hjust=1,size=2)+
  geom_vline(xintercept = c(-1,1), color="grey")+  geom_hline(yintercept = c(-1,1), color="grey")+
  ylim(c(min(diff_exp_all_celltype$avg_log2FC_fullcohort),max(diff_exp_all_celltype$avg_log2FC_fullcohort)))+
  xlim(c(min(diff_exp_all_celltype$avg_log2FC_fullcohort),max(diff_exp_all_celltype$avg_log2FC_fullcohort)))+
  annotate("text", x=1.5, y=3.75, label="Higher Expression in Adults")+
  annotate("text", x=-1.9, y=-3.99, label="Higher Expression in Peds")+
  geom_segment(aes(xend = avg_log2FC_fullcohort, yend = avg_log2FC_fullcohort),alpha=0.25) +
  geom_abline(slope=1, intercept=0, color="grey")
ggsave(file=here("figures/FC_correlation_KC_fullcohort_adultsubset_differential_Ped_adult.pdf"), h=7, w=8)
ggsave(file=here("figures/jpeg/FC_correlation_KC_fullcohort_adultsubset_differential_Ped_adult.jpeg"), h=7, w=8)



  
