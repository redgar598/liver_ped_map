### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(reshape2)
library(gtools)
library(SoupX)
library(colorspace)
library(cowplot)
library(DropletQC)
library(anndata)

source("scripts/00_pretty_plots.R")
source("scripts/00_entropy_d10x.R")

#################
## Load raw QC'ed data
#################
d10x_fetal<-readRDS(here("data","d10x_fetal_raw.rds"))
#d10x_fetal<-readRDS("/media/redgar/Seagate Portable Drive/fetal_liver/d10x_fetal_raw.rds")

d10x_ped_IFALD<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))
## add cell type labels
load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x_ped_IFALD), cell_label$index),]
identical(colnames(d10x_ped_IFALD), cell_label$index)
d10x_ped_IFALD <- AddMetaData(d10x_ped_IFALD, metadata = cell_label)


#################
## Match meta data
#################
d10x_fetal@meta.data$Barcodes<-NULL
colnames(d10x_fetal@meta.data)[which(colnames(d10x_fetal@meta.data)=="Extract.Name")]<-"Sample"
colnames(d10x_fetal@meta.data)[which(colnames(d10x_fetal@meta.data)=="Characteristics.individual.")]<-"individual"
d10x_fetal@meta.data$Treatment<-"Healthy"
d10x_fetal@meta.data$Tissue<-"TLH"
d10x_fetal@meta.data$chemistry<-"3pr"
colnames(d10x_fetal@meta.data)[which(colnames(d10x_fetal@meta.data)=="Characteristics.sex.")]<-"Sex"
colnames(d10x_fetal@meta.data)[which(colnames(d10x_fetal@meta.data)=="Characteristics.age.")]<-"Age"
d10x_fetal@meta.data$AgeGroup<-"Fetal"
d10x_fetal@meta.data$FreshorFrozen<-"fresh"
d10x_fetal@meta.data$BMI<-NA
d10x_fetal@meta.data$relALBChange<-NA
d10x_fetal@meta.data$nuclear_fraction<-NA
d10x_fetal@meta.data$cell_status<-NA
colnames(d10x_fetal@meta.data)[which(colnames(d10x_fetal@meta.data)=="Cell.Labels")]<-"CellType_refined"
d10x_fetal@meta.data$age_condition<-paste(d10x_fetal$AgeGroup, d10x_fetal$Treatment, sep=" ")


d10x_ped_IFALD@meta.data$file<-NULL
d10x_ped_IFALD@meta.data$Approx_bam_GB<-NULL
d10x_ped_IFALD@meta.data$Characteristics.facs.sorting.<-NA
d10x_ped_IFALD@meta.data$Sample<-d10x_ped_IFALD@meta.data$individual
d10x_ped_IFALD@meta.data$individual<-sapply(1:nrow(d10x_ped_IFALD@meta.data), function(x) strsplit(d10x_ped_IFALD@meta.data$individual[x],"_")[[1]][1])
d10x_ped_IFALD@meta.data$CellType_rough <-NULL
d10x_ped_IFALD@meta.data$second_best_cell  <-NULL
d10x_ped_IFALD@meta.data$S.Score   <-NULL
d10x_ped_IFALD@meta.data$G2M.Score <-NULL
d10x_ped_IFALD@meta.data$Phase     <-NULL
d10x_ped_IFALD@meta.data$old.ident<-NULL
d10x_ped_IFALD@meta.data$integrated_snn_res.0.5 <-NULL
d10x_ped_IFALD@meta.data$seurat_clusters<-NULL
d10x_ped_IFALD@meta.data$age_id<-NULL
d10x_ped_IFALD@meta.data$index<-NULL


head(d10x_ped_IFALD@meta.data)
head(d10x_fetal@meta.data)

d10x_fetal@meta.data<-d10x_fetal@meta.data[,colnames(d10x_ped_IFALD@meta.data)]

head(d10x_fetal@meta.data)
head(d10x_ped_IFALD@meta.data)






#################
## Merge
#################
d10x <- merge(d10x_ped_IFALD,d10x_fetal, merge.data=TRUE, project = "IFALD_fetal_adult_ped_map")
saveRDS(d10x, file = here("data","Fetal_IFALD_d10x_adult_ped_raw.rds"))
# 
# #########
# ## Signature Scores on unintegrated data
# #########
# 
# ##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# # This is log(TP10K+1)
# d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")
# 
# myeloid_immune_supressive<-c("CTSB","CD163","MS4A7","FOLR2","GPNMB","VSIG4","HMOX1","MSR1")
# inflammatory_macs<-c("CD74","HLA-DRA","C1QC","HLA-DPA1","HLA-DPB1","LYZ","S100A6")
# exhausted_tcells<-c("TOX","PDCD1","LAG3","TNFRSF9","CXCL13","ENTPD1","HAVCR2","CD38")
# 
# recent_recruit_myeloid<-c("S100A8","S100A9","CD68","LYZ")
# kuffer_signature<-c("VSIG4","MARCO","CD5L","HMOX1")
# 
# ######
# ## Score Signatures
# ######
# d10x <- AddModuleScore(
#   object = d10x,
#   features = list(myeloid_immune_supressive),
#   ctrl = 5,
#   name = 'myeloid_immune_supressive_score'
# )
# 
# d10x <- AddModuleScore(
#   object = d10x,
#   features = list(inflammatory_macs),
#   ctrl = 5,
#   name = 'inflammatory_macs_score'
# )
# 
# d10x <- AddModuleScore(
#   object = d10x,
#   features = list(exhausted_tcells),
#   ctrl = 5,
#   name = 'exhausted_tcells_score'
# )
# 
# d10x <- AddModuleScore(
#   object = d10x,
#   features = list(recent_recruit_myeloid),
#   ctrl = 5,
#   name = 'recently_recruited_myeloid'
# )
# 
# d10x <- AddModuleScore(
#   object = d10x,
#   features = list(kuffer_signature),
#   ctrl = 5,
#   name = 'kuffer_like_score'
# )
# 
# score_data<-d10x@meta.data[,c("myeloid_immune_supressive_score1","inflammatory_macs_score1","exhausted_tcells_score1","recently_recruited_myeloid1","kuffer_like_score1")]
# score_data$cell<-rownames(score_data)
# 
# ###############
# ## Integrate by donor
# ###############
# #https://satijalab.org/seurat/articles/integration_rpca.html
# print("RUNNING INTEGRATION")
# 
# ## run integration across donor and hopefully that will also smooth out differences with chemistry and batch?
# d10x.list<- SplitObject(d10x, split.by = "Sample")
# 
# # normalize, identify variable features and score cell cycle for each dataset independently
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# 
# d10x.list <- lapply(X = d10x.list, FUN = function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
#   x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# })
# 
# # select features that are repeatedly variable across datasets for integration run PCA on each
# # dataset using these features
# features <- SelectIntegrationFeatures(object.list = d10x.list)
# d10x.list <- lapply(X = d10x.list, FUN = function(x) {
#   #x <- ScaleData(x, features = features, verbose = FALSE)
#   x <- ScaleData(x, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), features = features, verbose = FALSE)
#   x <- RunPCA(x, features = features, verbose = FALSE)
# })
# 
# ## Identify anchors
# chem.anchors <- FindIntegrationAnchors(object.list = d10x.list, anchor.features = features, reduction = "rpca")
# d10x.fetal_ped_IFALD <- IntegrateData(anchorset = chem.anchors)
# 
# DefaultAssay(d10x.fetal_ped_IFALD) <- "integrated"
# 
# print("INTEGRATED")
# 
# 
# # Run the standard workflow for visualization and clustering
# d10x.fetal_ped_IFALD <- ScaleData(d10x.fetal_ped_IFALD, verbose = FALSE)
# d10x.fetal_ped_IFALD <- RunPCA(d10x.fetal_ped_IFALD, npcs = 30, verbose = FALSE)
# d10x.fetal_ped_IFALD <- RunUMAP(d10x.fetal_ped_IFALD, reduction = "pca", dims = 1:30)
# d10x.fetal_ped_IFALD <- RunTSNE(d10x.fetal_ped_IFALD, dims = 1:30)
# 
# d10x.fetal_ped_IFALD <- FindNeighbors(d10x.fetal_ped_IFALD, reduction = "pca", dims = 1:30)
# d10x.fetal_ped_IFALD <- FindClusters(d10x.fetal_ped_IFALD, resolution = 0.5)
# 
# d10x.fetal_ped_IFALD
# 
# 
# ##############
# ## Save integrated to look local
# ##############
# save(d10x.fetal_ped_IFALD, file=paste(here("data/"),"Fetal_IFALD_adult_ped_integrated.rds", sep=""))
# cell_label<-d10x.fetal_ped_IFALD@meta.data
# save(cell_label, file=paste(here("data/"),"Fetal_IFALD_adult_ped_cellRough.rds", sep=""))
# 
# 
# 
# 
# 
# 
# ###################
# ### save for local plotting
# ###################
# load(here("data","Fetal_IFALD_adult_ped_integrated.rds"))
# 
Macrophage_genes<-c( "PTPRC", "MARCO","CD74")
LEC_genes<-c("CALCRL","RAMP2")
Hepatocyte_genes<-c("ALB", "CYP3A4")
Cholangiocytes_genes<-c( "EPCAM", "KRT7")
HSCs_genes<-c( "IGFBP7",  "SPARC")
T_genes<-c("CD3D","CD8A")
NK_genes<-c("NKG7","CD7")
gd_genes<-c("GNLY")
RBC<-c("HBB","HBA2","HBA1","FCGR3A")
MAST<-c("TPSAB1", "AREG")
recent_recruit_myeloid<-c("S100A8","S100A9","CD68","LYZ")
kuffer_signature<-c("VSIG4","CD5L")
neutro_gene<-c("CSF3R","FCGR3B")
MHCII<-c("HLA-DRA","HLA-DPB1")
b_genes_noIG<-c("MS4A1", "CD79B")
immunoglobins<-c("IGKC","IGHG1")
# 
# head(d10x.fetal_ped_IFALD)
# 
# umap_mat<-as.data.frame(Embeddings(object = d10x.fetal_ped_IFALD, reduction = "umap"))#
# umap_mat$cell<-rownames(umap_mat)
# 
# meta<-d10x.fetal_ped_IFALD@meta.data
# meta$cell<-rownames(meta)
# 
# plt<-merge(meta, umap_mat, by="cell")
# 
# ## raw expression values
# gene_exp<-FetchData(d10x, vars=c(Macrophage_genes,LEC_genes,Hepatocyte_genes,Cholangiocytes_genes,HSCs_genes,T_genes,NK_genes,gd_genes,RBC,
#                                  MAST, recent_recruit_myeloid, kuffer_signature, neutro_gene, MHCII, b_genes_noIG, immunoglobins))
# gene_exp$cell<-rownames(gene_exp)
# gene_exp<-melt(gene_exp)
# 
# plt<-merge(gene_exp,plt, by="cell")
# 
# plt<-merge(score_data,plt, by="cell")
# 
# head(plt)
# 
# 
# 
# 
# save(plt, file=here("data","Fetal_IFALD_adult_ped_pltData.RData"))
# 
# 

#############
## PCAs for plotting
#############
d10x$CellType_harmonized<-d10x$CellType_refined

levels(d10x$CellType_harmonized)[which(levels(d10x$CellType_harmonized)=="Kupffer Cell")]<-"KC Like"
levels(d10x$CellType_harmonized)[which(levels(d10x$CellType_harmonized)%in%c("Mono-Mac","Monocyte-DC precursor","Monocyte" ))]<-"RR Myeloid"
levels(d10x$CellType_harmonized)[which(levels(d10x$CellType_harmonized)%in%c("DC2","DC1" ))]<-"Macrophage\n(MHCII high)"
levels(d10x$CellType_harmonized)[which(levels(d10x$CellType_harmonized)%in%c("NK-like cells","NK" ))]<-"NK cell"

d10x.combined_myeloid<-subset(d10x, subset = CellType_harmonized %in% c("RR Myeloid","KC Like","Macrophage\n(MHCII high)"))
rm(d10x)
gc

d10x.combined_myeloid <- NormalizeData(d10x.combined_myeloid)
d10x.combined_myeloid <- FindVariableFeatures(d10x.combined_myeloid, selection.method = "vst", nfeatures = 2000)
d10x.combined_myeloid <- ScaleData(d10x.combined_myeloid) 
d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)

d10x.combined_myeloid[["pca"]]

head(Embeddings(d10x.combined_myeloid, reduction = "pca")[, 1:5])
head(Loadings(d10x.combined_myeloid, reduction = "pca")[, 1:5])
head(Stdev(d10x.combined_myeloid, reduction = "pca"))

#' ## PCA for batch effect
Loadings<-as.data.frame(Loadings(d10x.combined_myeloid, reduction = "pca"))
embed<-as.data.frame(Embeddings(d10x.combined_myeloid, reduction = "pca"))
vars <- Stdev(d10x.combined_myeloid, reduction = "pca")^2
Importance<-vars/sum(vars)
print(Importance[1:10])

meta_categorical <- d10x.combined_myeloid@meta.data[, c("CellType_refined","CellType_harmonized","age_condition")]  # input column numbers in meta that contain categorical variables
meta_continuous <- d10x.combined_myeloid@meta.data[, c("percent.mt","nFeature_RNA","Age")]  # input column numbers in meta that contain continuous variables

save(embed,vars, Importance, meta_categorical, meta_continuous,Loadings, file=here("data","Fetal_ped_IFALD_adult_PCA_myeloid.RData"))

# 
# #############
# ### Plots
# #############
# #load(here("data","Fetal_IFALD_adult_ped_pltData.RData"))
# load(here("/media/redgar/Seagate Portable Drive/fetal_liver/","Fetal_IFALD_adult_ped_pltData.RData"))
# 
# plt_not_gene<-plt[which(plt$variable=="MARCO"),]
# plt_not_gene[which(is.na(plt_not_gene$CellType_refined)),]
# plt[which(is.na(plt$variable)),]
# plt_not_gene$CellType_refined<-factor(plt_not_gene$CellType_refined, levels=names(combo_colors))
# 
# 
# all_UMAP<-ggplot(plt_not_gene, aes(UMAP_1,UMAP_2, color=CellType_refined))+
#   geom_point(size=0.5)+
#   theme_classic()+th+theme(legend.text=element_text(size=10),
#                            legend.title=element_text(size=12),
#                            plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#   guides(colour = guide_legend(override.aes = list(size=5)))+
#   annotate("text", x = -10, y = -15, label = paste0("n = ",comma(nrow(plt_not_gene))))+
#   scale_color_manual(name="Cell Type",values = combo_colors, drop = T, limits=force)
# save_plts(all_UMAP, "UMAP_Fetal_ped_adult_IFALD", w=20,h=12)
# 
# 
# countcells<-data.frame(tapply(plt_not_gene$cell, plt_not_gene$age_condition, function(x) length(unique(x))))
# colnames(countcells)<-c("Count")
# countcells$age_condition<-rownames(countcells)
# 
# 
# condition_split_UMAP<-ggplot(plt_not_gene, aes(UMAP_1,UMAP_2))+
#   geom_point(aes( color=CellType_refined),size=0.5)+
#   theme_classic()+th+theme(legend.text=element_text(size=10),
#                            legend.title=element_text(size=12),
#                            plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#   guides(colour = guide_legend(override.aes = list(size=5)))+
#   geom_text(aes(x = -10, y = -15, label=paste0("n = ",comma(Count))), countcells)+
#   scale_color_manual(name="Cell Type",values = combo_colors, drop = T, limits=force)+
#   facet_wrap(~age_condition)
# save_plts(condition_split_UMAP, "condition_split_UMAP_Fetal_ped_adult_IFALD", w=20,h=12)
# 
# 
# all_UMAP<-ggplot(plt_not_gene, aes(UMAP_1,UMAP_2, color=seurat_clusters))+
#   geom_point(size=0.5)+
#   theme_classic()+th+theme(legend.text=element_text(size=10),
#                            legend.title=element_text(size=12),
#                            plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#   guides(colour = guide_legend(override.aes = list(size=5)))+
#   annotate("text", x = -10, y = -15, label = paste0("n = ",comma(nrow(plt_not_gene))))
# save_plts(all_UMAP, "UMAP_Fetal_ped_adult_IFALD_seuratclusters", w=16,h=12)
# 
# 
# ##############
# ## entropy 
# ##############
# entropy_plt<-function(covariate, d10x){
#   entrophy_cluster_df<-do.call(rbind, lapply(as.character(unique(d10x$seurat_clusters)), function(cluster){
#     data.frame(seurat_clusters=cluster, entropy=entropy(d10x[,covariate][which(d10x$seurat_clusters==cluster)]))
#   }))
# 
#   cell_cluster_count<-plt_not_gene %>%  group_by(seurat_clusters,rlang::parse_expr(covariate) %>% eval) %>%
#     summarise(n = n()) %>%
#     mutate(freq = n / sum(n))
# 
#   cell_cluster_count<-as.data.frame(cell_cluster_count)
#   cell_cluster_count<-merge(cell_cluster_count, entrophy_cluster_df, by="seurat_clusters")
#   colnames(cell_cluster_count)[2]<-covariate
#   plt_entropy<-cell_cluster_count
# 
#   max_count<-as.data.frame(plt_entropy %>%
#                              group_by(seurat_clusters) %>%
#                              summarise(total = sum(n)))
# 
#   entrophy_label<-plt_entropy[!duplicated(plt_entropy[,c("seurat_clusters",  "entropy")]), c("seurat_clusters",  "entropy")]
#   entrophy_label<-merge(entrophy_label, max_count, by="seurat_clusters")
# 
#   bar_individual<-ggplot() +
#     geom_bar(aes(fill=rlang::parse_expr(covariate) %>% eval, y=n, x=seurat_clusters),plt_entropy, position="stack", stat="identity", color="black")+
#     theme_bw()+th+ylab("Cell Count")+xlab("Seurat Cluster")+
#     geom_text(aes(label=round(entropy,2), y=(total+(0.1*max(total))), x=seurat_clusters), entrophy_label)+
#     scale_fill_manual(values=sequential_hcl(length(unique(d10x[,covariate])), "Batlow"), name=covariate)
# 
#   individual_UMAP<-ggplot(d10x, aes(UMAP_1,UMAP_2, color=get(covariate)))+geom_point(size=0.5)+theme_void()+
#     scale_color_manual(values=sequential_hcl(length(unique(d10x[,covariate])), "Batlow"))+
#     theme(legend.position = "none")
# 
#   label.df_2 <- d10x %>%
#     group_by(seurat_clusters) %>%
#     summarize(x = mean(UMAP_1), y = mean(UMAP_2))
# 
#   cluster_UMAP<-ggplot(d10x, aes(UMAP_1,UMAP_2, color=seurat_clusters))+geom_point(size=0.5)+theme_void()+ theme(legend.position = "none") +
#     geom_text(aes(x,y,label=seurat_clusters),label.df_2, color="black")
# 
#   plot_grid(plot_grid(cluster_UMAP, individual_UMAP), bar_individual, ncol=1, rel_heights = c(1.5,1))}
# 
# 
# entropy_plt("age_condition", plt_not_gene)
# save_plts(entropy_plt("age_condition", plt_not_gene), "IFALD_Fetal_entropy_condition", w=15,h=10)
# 
# 
# ### with color pallette
# covariate<-"CellType_refined"
# d10x<-plt_not_gene
# 
# entrophy_cluster_df<-do.call(rbind, lapply(as.character(unique(d10x$seurat_clusters)), function(cluster){
#   data.frame(seurat_clusters=cluster, entropy=entropy(d10x[,covariate][which(d10x$seurat_clusters==cluster)]))
# }))
# 
# cell_cluster_count<-plt_not_gene %>%  group_by(seurat_clusters,rlang::parse_expr(covariate) %>% eval) %>%
#   summarise(n = n()) %>%
#   mutate(freq = n / sum(n))
# 
# cell_cluster_count<-as.data.frame(cell_cluster_count)
# cell_cluster_count<-merge(cell_cluster_count, entrophy_cluster_df, by="seurat_clusters")
# colnames(cell_cluster_count)[2]<-covariate
# plt_entropy<-cell_cluster_count
# 
# max_count<-as.data.frame(plt_entropy %>%
#                            group_by(seurat_clusters) %>%
#                            summarise(total = sum(n)))
# 
# entrophy_label<-plt_entropy[!duplicated(plt_entropy[,c("seurat_clusters",  "entropy")]), c("seurat_clusters",  "entropy")]
# entrophy_label<-merge(entrophy_label, max_count, by="seurat_clusters")
# 
# bar_individual<-ggplot() +
#   geom_bar(aes(fill=rlang::parse_expr(covariate) %>% eval, y=n, x=seurat_clusters),plt_entropy, position="stack", stat="identity", color="black")+
#   theme_bw()+th+ylab("Cell Count")+xlab("Seurat Cluster")+
#   geom_text(aes(label=round(entropy,2), y=(total+(0.1*max(total))), x=seurat_clusters), entrophy_label)+
#   scale_fill_manual(values=combo_colors, name=covariate)
# 
# individual_UMAP<-ggplot(d10x, aes(UMAP_1,UMAP_2, color=get(covariate)))+geom_point(size=0.5)+theme_void()+
#   scale_color_manual(values=combo_colors)+
#   theme(legend.position = "none")
# 
# label.df_2 <- d10x %>%
#   group_by(seurat_clusters) %>%
#   summarize(x = mean(UMAP_1), y = mean(UMAP_2))
# 
# cluster_UMAP<-ggplot(d10x, aes(UMAP_1,UMAP_2, color=seurat_clusters))+geom_point(size=0.5)+theme_void()+ theme(legend.position = "none") +
#   geom_text(aes(x,y,label=seurat_clusters),label.df_2, color="black")
# 
# entropy_celltype<-plot_grid(plot_grid(cluster_UMAP, individual_UMAP), bar_individual, ncol=1, rel_heights = c(1.5,1))
# 
# save_plts(entropy_celltype, "IFALD_Fetal_entropy_celltype", w=20,h=15)
# 
# 
# 
# 
# 
# ###############
# ## summarize dotplot
# ###############
# scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}
# 
# plt_summary<-plt %>% group_by(CellType_refined, variable) %>%
#   summarise(mn=mean(value), count=length(which(value>0)), percent_exp=(length(which(value>0))/length(value))*100)
# plt_summary <- plt_summary %>% group_by(variable) %>%
#   dplyr::mutate(scaled = scale_this(mn))
# plt_summary<-as.data.frame(plt_summary)
# 
# # remove dots where 0 cell expressing marker
# plt_summary<-plt_summary[(which(plt_summary$count>0)),]
# 
# plt_summary$variable<-gsub("rna_","",plt_summary$variable)
# plt_summary$variable<-factor(plt_summary$variable, levels=rev(c(b_genes_noIG, immunoglobins, T_genes,gd_genes,NK_genes, Cholangiocytes_genes,LEC_genes,
#                                                                 Macrophage_genes,recent_recruit_myeloid,MHCII, kuffer_signature,RBC, neutro_gene, MAST,
#                                                                 HSCs_genes,Hepatocyte_genes)))
# 
# plt_summary$CellType_refined<-factor(plt_summary$CellType_refined, levels=c("B-cells","Mature B-cells","Plasma cells","pre B cell","pro B cell","B cell","pre pro B cell ","Platelets",
#                                                                             "CD3+ T-cells","gd T-cells","NK-like cells", "NK and T cells","Early lymphoid/T lymphocyte","CLNK T-cells","Cycling T-cells","NK","ILC precursor","NKT cells",
#                                                                             "Cholangiocytes", "LSEC","Endothelial cell","Fibroblast",
#                                                                             "Myeloid cells","RR Myeloid","Mono-Mac","Kupffer Cell","KC Like","Macrophage\n(MHCII high)","Macrophage\n(CLEC9A high)",
#                                                                             "Monocyte","pDC precursor","Megakaryocyte",
#                                                                             "Cycling Myeloid","DC1","DC2","Monocyte-DC precursor","Neutrophil-myeloid progenitor","Mono-NK",
#                                                                             "Myeloid Erythrocytes\n(phagocytosis)","VCAM1+ Erythroblastic Island Macrophage","Erythroblastic Island Macrophage",
#                                                                             "Early Erythroid","Mid  Erythroid","Late Erythroid","Erythrocytes","MEMP",
#                                                                             "Neutrophil","Neutrophil\n(DEFA+)",
#                                                                             "Mast cell",
#                                                                             "HSC","HSC/MPP",
#                                                                             "Hepatocytes","Hepatocyte",
#                                                                             "Doublet","Low Quality"))
# 
# 
# fancy_dotplot<-plot_grid(
#   ggplot(plt_summary, aes(CellType_refined, variable, color=scaled, size=percent_exp))+geom_point()+
#     th+theme_classic()+
#     scale_color_continuous_sequential(palette = "Oslo", rev=F, name="Scaled\nMean\nExpression")+
#     scale_size(name="Percent\nCells\nExpressing")+
#     theme(axis.text.x = element_blank(),axis.title = element_blank(),axis.ticks.x = element_blank(),
#           plot.margin = margin(0.25,0.25,0,0.25,"cm"))+
#     geom_hline(yintercept = 32.5)+ geom_hline(yintercept = 27.5)+ geom_hline(yintercept = 23.5)+
#     geom_hline(yintercept = 12.5)+ geom_hline(yintercept = 8.5)+ geom_hline(yintercept = 6.5)+
#     geom_hline(yintercept = 4.5)+ geom_hline(yintercept = 2.5),
#   ggplot(plt_summary, aes(CellType_refined, y=1, fill=CellType_refined))+geom_tile(color="black")+
#     th+theme_classic()+  scale_fill_manual(name="Cell Type",
#                                             values = c(myColors_celltype,myColors_celltype_fetal), drop = T, limits=force)+
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#           axis.text.y = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
#           legend.position = "none",axis.line  = element_blank(),
#           plot.margin = margin(0,0,1,1,"cm")),
#   ncol=1, rel_heights = c(6,2), align = "v", axis="lr")
# fancy_dotplot
# save_plts(fancy_dotplot, "dot_plot_celltype_Fetal_ped_adult_IFALD", w=12,h=13)
# 
# 
# 
# 
#################
## individual gene expression
#################
# countcells<-data.frame(tapply(plt_not_gene$cell, plt_not_gene$age_condition, function(x) length(unique(x))))
# colnames(countcells)<-c("Count")
# countcells$age_condition<-rownames(countcells)
# 
# 
# gene_UMAP<-function(gene,percentile, splt=F) {
#   plt_not_gene<-plt[which(plt$variable==gene),]
# 
#   exp_limit<-quantile(plt_not_gene$value, percentile)
#   plt_not_gene$gene_exp_limited<-NA
#   plt_not_gene$gene_exp_limited[which(plt_not_gene$value>exp_limit)]<-plt_not_gene$value[which(plt_not_gene$value>exp_limit)]
#   #plt_not_gene<-plt_not_gene[rev(order(plt_not_gene$gene_exp_limited)),]
# 
#   plt_not_gene<-rbind(plt_not_gene[which(is.na(plt_not_gene$gene_exp_limited)),],
#                       plt_not_gene[which(!(is.na(plt_not_gene$gene_exp_limited))),][(order(plt_not_gene[which(!(is.na(plt_not_gene$gene_exp_limited))),]$gene_exp_limited)),])
# 
# 
#   if(splt==F){
#     umap_plt<-ggplot(plt_not_gene, aes(UMAP_1,UMAP_2, color=log(gene_exp_limited)))+
#       geom_point(size=0.15)+
#       theme_classic()+th+theme(legend.text=element_text(size=10),
#                                legend.title=element_text(size=12),
#                                plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#       annotate("text", x = -10, y = -15, label = paste0("n = ",comma(nrow(plt_not_gene))))+
#       scale_color_continuous_sequential(palette = "Viridis", rev=F, nam=paste(gene, "\nExpression\n(log)"),na.value = "grey80")
#     }
#   if(splt==T){
#     umap_plt<-ggplot(plt_not_gene, aes(UMAP_1,UMAP_2))+
#       geom_point(size=0.15, aes(color=log(gene_exp_limited)))+  facet_wrap(~age_condition)+
#       theme_classic()+th+theme(legend.text=element_text(size=10),
#                                legend.title=element_text(size=12),
#                                plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#       geom_text(aes(x = -10, y = -15, label=paste0("n = ",comma(Count))), countcells)+
#       scale_color_continuous_sequential(palette = "Viridis", rev=F, nam=paste(gene, "\nExpression\n(log)"),na.value = "grey80")
#   }
#   umap_plt
#   }
# 
# gene_UMAP("ALB", 0)
# gene_UMAP("ALB", 0.1)
# gene_UMAP("ALB", 0.9)
# 
# gene_UMAP("VSIG4", 0.9)
# gene_UMAP("CD5L", 0.9)
# gene_UMAP("VSIG4", 0.9,T)
# gene_UMAP("CD5L", 0.9,T)
# 
# gene_UMAP("S100A8", 0.9,T)
# gene_UMAP("S100A9", 0.9,T)
# gene_UMAP("CD68", 0.9,T)
# gene_UMAP("LYZ", 0.9,T)
# 
# gene_UMAP("HLA-DRA", 0.9,T)
# gene_UMAP("HLA-DPB1", 0.9,T)
# gene_UMAP("HLA-DRA", 0.9)
# gene_UMAP("HLA-DPB1", 0.9)
# 
# # 
# # save_plts(all_UMAP, "UMAP_Fetal_ped_adult_IFALD", w=20,h=12)
# # 
# # 
# #############
# ## Cell Label Harmonization
# #############
# #load(here("/media/redgar/Seagate Portable Drive/fetal_liver/","Fetal_IFALD_adult_ped_pltData.RData"))
# 
# plt_not_gene<-plt[which(plt$variable=="MARCO"),]
# plt_not_gene$CellType_refined<-factor(plt_not_gene$CellType_refined, levels=names(combo_colors))
# 
# plt_not_gene$CellType_harmonized<-plt_not_gene$CellType_refined
# 
# levels(plt_not_gene$CellType_harmonized)[which(levels(plt_not_gene$CellType_harmonized)=="Kupffer Cell")]<-"KC Like"
# levels(plt_not_gene$CellType_harmonized)[which(levels(plt_not_gene$CellType_harmonized)%in%c("Mono-Mac","Monocyte-DC precursor","Monocyte" ))]<-"RR Myeloid"
# levels(plt_not_gene$CellType_harmonized)[which(levels(plt_not_gene$CellType_harmonized)%in%c("DC2","DC1" ))]<-"Macrophage\n(MHCII high)"
# levels(plt_not_gene$CellType_harmonized)[which(levels(plt_not_gene$CellType_harmonized)%in%c("NK-like cells","NK" ))]<-"NK cell"
# 
# levels(plt_not_gene$CellType_harmonized)
# 
# table(plt_not_gene$CellType_refined,plt_not_gene$seurat_clusters)
# table(plt_not_gene$CellType_refined[which(plt_not_gene$seurat_clusters=="3")])[order(table(plt_not_gene$CellType_refined[which(plt_not_gene$seurat_clusters=="3")]))]
# table(plt_not_gene$CellType_harmonized)[order(table(plt_not_gene$CellType_harmonized))]
# 
# 
# ggplot(plt_not_gene, aes(UMAP_1,UMAP_2, color=CellType_harmonized))+
#     geom_point(size=0.5)+
#     theme_classic()+th+theme(legend.text=element_text(size=10),
#                              legend.title=element_text(size=12),
#                              plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#     guides(colour = guide_legend(override.aes = list(size=5)))+
#     annotate("text", x = -10, y = -15, label = paste0("n = ",comma(nrow(plt_not_gene))))+
#     scale_color_manual(name="Cell Type",values = combo_colors, drop = T, limits=force)
# 
# 
# #############
# ## Myeloid Score plots
# #############
# plt_not_gene_myeloid<-plt_not_gene[which(plt_not_gene$CellType_harmonized%in%c("Cycling Myeloid","RR Myeloid","KC Like","Macrophage\n(MHCII high)")),]
# 
# myeloid_UMAP<-ggplot(plt_not_gene_myeloid, aes(UMAP_1,UMAP_2, color=CellType_harmonized))+
#   geom_point(size=0.5)+
#   theme_classic()+th+theme(legend.text=element_text(size=10),
#                            legend.title=element_text(size=12),
#                            plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#   guides(colour = guide_legend(override.aes = list(size=5)))+
#   annotate("text", x = -10, y = -15, label = paste0("n = ",comma(nrow(plt_not_gene_myeloid))))+
#   scale_color_manual(name="Cell Type",values = combo_colors, drop = T, limits=force)
# myeloid_UMAP
# save_plts(myeloid_UMAP, "UMAP_Fetal_ped_adult_IFALD_myeloid", w=20,h=12)
# 
# recruit_all_myeloid<-ggplot(plt_not_gene_myeloid, aes(UMAP_1,UMAP_2, color=recently_recruited_myeloid1))+
#   geom_point(size=0.5)+
#   theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#   annotate("text", x = -5, y = -6, label = paste0("n = ",comma(nrow(plt_not_gene_myeloid))))+
#   scale_color_continuous_sequential(palette = "Viridis", rev=F, name="Recently\nRecruited\nMyeloid\nSignature Score")
# recruit_all_myeloid
# save_plts(recruit_all_myeloid, "Fetal_adult_ped_IFALD_recruit_myeloid", w=7,h=5)
# 
# 
# countcells<-data.frame(tapply(plt_not_gene_myeloid$cell, plt_not_gene_myeloid$age_condition, function(x) length(unique(x))))
# colnames(countcells)<-c("Count")
# countcells$age_condition<-rownames(countcells)
# 
# condition_split_UMAP_myeloid<-ggplot(plt_not_gene_myeloid, aes(UMAP_1,UMAP_2))+
#   geom_point(aes( color=CellType_harmonized),size=0.5)+
#   theme_classic()+th+theme(legend.text=element_text(size=10),
#                            legend.title=element_text(size=12),
#                            plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#   guides(colour = guide_legend(override.aes = list(size=5)))+
#   geom_text(aes(x = -10, y = -15, label=paste0("n = ",comma(Count))), countcells)+
#   scale_color_manual(name="Cell Type",values = combo_colors, drop = T, limits=force)+
#   facet_wrap(~age_condition)
# condition_split_UMAP_myeloid
# save_plts(condition_split_UMAP_myeloid, "condition_split_UMAP_Fetal_ped_adult_IFALD_myeloid", w=20,h=12)
# 
# 
# condition_split_UMAP<-ggplot(plt_not_gene_myeloid, aes(UMAP_1,UMAP_2))+
#   geom_point(aes( color=recently_recruited_myeloid1),size=0.5)+
#   theme_classic()+th+theme(legend.text=element_text(size=10),
#                            legend.title=element_text(size=12),
#                            plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#   scale_color_continuous_sequential(palette = "Viridis", rev=F, name="Recently\nRecruited\nMyeloid\nSignature Score")+
#   facet_wrap(~age_condition)+
#   geom_text(aes(x = -10, y = -15, label=paste0("n = ",comma(Count))), countcells)
# condition_split_UMAP
# save_plts(condition_split_UMAP, "Fetal_adult_ped_IFALD_recruit_myeloid_condition_split", w=20,h=12)
# 
# condition_split_UMAP<-ggplot(plt_not_gene_myeloid, aes(UMAP_1,UMAP_2))+
#   geom_point(aes( color=kuffer_like_score1),size=0.5)+
#   theme_classic()+th+theme(legend.text=element_text(size=10),
#                            legend.title=element_text(size=12),
#                            plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#   scale_color_continuous_sequential(palette = "Viridis", rev=F, name="Kupffer-like\nSignature Score")+
#   facet_wrap(~age_condition)+
#   geom_text(aes(x = -10, y = -15, label=paste0("n = ",comma(Count))), countcells)
# condition_split_UMAP
# save_plts(condition_split_UMAP, "Fetal_adult_ped_IFALD_KC_myeloid_condition_split", w=20,h=12)
# 
# condition_split_UMAP<-ggplot(plt_not_gene_myeloid, aes(UMAP_1,UMAP_2))+
#   geom_point(aes( color=inflammatory_macs_score1),size=0.5)+
#   theme_classic()+th+theme(legend.text=element_text(size=10),
#                            legend.title=element_text(size=12),
#                            plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#   scale_color_continuous_sequential(palette = "Viridis", rev=F, name="MHCII\nMyeloid\nSignature Score")+
#   facet_wrap(~age_condition)+
#   geom_text(aes(x = -10, y = -15, label=paste0("n = ",comma(Count))), countcells)
# condition_split_UMAP
# save_plts(condition_split_UMAP, "Fetal_adult_ped_IFALD_MHCII_myeloid_condition_split", w=20,h=12)
# 
# 
# ################
# ## Statistics
# ################
# plt_myeloid_adult<-plt_not_gene_myeloid[which(plt_not_gene_myeloid$age_condition=="Adult Healthy"),]
# adult_cellcount<-nrow(plt_myeloid_adult)
# print(paste(adult_cellcount, " adult myeloid cells",sep=""))
# 
# plt_myeloid_ped<-plt_not_gene_myeloid[which(plt_not_gene_myeloid$age_condition=="Ped Healthy"),]
# ped_cellcount<-nrow(plt_myeloid_ped)
# print(paste(ped_cellcount, " ped myeloid cells",sep=""))
# 
# plt_myeloid_fetal<-plt_not_gene_myeloid[which(plt_not_gene_myeloid$age_condition=="Fetal Healthy"),]
# fetal_cellcount<-nrow(plt_myeloid_fetal)
# print(paste(fetal_cellcount, " fetal myeloid cells",sep=""))
# 
# 
# samp_num<-10000
# pval<-0.05
# myeloid_pval_montecarlo<-do.call(rbind, lapply(1:samp_num, function(x){
#   set.seed(x)
#   plt_myeloid_adult_random<-plt_myeloid_adult[sample(adult_cellcount,ped_cellcount),]
# 
#   supressive<-t.test(plt_myeloid_adult_random$myeloid_immune_supressive_score1, plt_myeloid_ped$myeloid_immune_supressive_score1)$p.value
#   inflammatory<-t.test(plt_myeloid_adult_random$inflammatory_macs_score1, plt_myeloid_ped$inflammatory_macs_score1)$p.value
#   recruit<-t.test(plt_myeloid_adult_random$recently_recruited_myeloid1, plt_myeloid_ped$recently_recruited_myeloid1)$p.value
#   kuffer<-t.test(plt_myeloid_adult_random$kuffer_like_score1, plt_myeloid_ped$kuffer_like_score1)$p.value
# 
#   data.frame(supressive=supressive, inflammatory=inflammatory, recruit=recruit ,kuffer=kuffer)
# }))
# 
# print(paste("Comparing the myeloid immune supressive score in myeloid cells, there is a sig difference between ped and adult (p <", pval,
#             ") in ", samp_num, " random samples at a p value of ",
#             round((length(myeloid_pval_montecarlo$supressive[which(myeloid_pval_montecarlo$supressive>pval)])+1)/(samp_num+1), 3), sep=""))
# print(paste("Comparing the inflammatory macs score in myeloid cells, there is a sig difference between ped and adult (p <", pval,
#             ") in ", samp_num, " random samples at a p value of ",
#             round((length(myeloid_pval_montecarlo$inflammatory[which(myeloid_pval_montecarlo$inflammatory>pval)])+1)/(samp_num+1), 3), sep=""))
# print(paste("Comparing the recently recruited myeloid score in myeloid cells, there is a sig difference between ped and adult (p <", pval,
#             ") in ", samp_num, " random samples at a p value of ",
#             round((length(myeloid_pval_montecarlo$recruit[which(myeloid_pval_montecarlo$recruit>pval)])+1)/(samp_num+1), 6), sep=""))
# print(paste("Comparing the kuffer-like score in myeloid cells, there is a sig difference between ped and adult (p <", pval,
#             ") in ", samp_num, " random samples at a p value of ",
#             round((length(myeloid_pval_montecarlo$kuffer[which(myeloid_pval_montecarlo$kuffer>pval)])+1)/(samp_num+1), 6), sep=""))
# 
# 
# samp_num<-10000
# pval<-0.05
# myeloid_pval_montecarlo<-do.call(rbind, lapply(1:samp_num, function(x){
#   set.seed(x)
#   plt_myeloid_fetal_random<-plt_myeloid_fetal[sample(fetal_cellcount,ped_cellcount),]
# 
#   supressive<-t.test(plt_myeloid_fetal_random$myeloid_immune_supressive_score1, plt_myeloid_ped$myeloid_immune_supressive_score1)$p.value
#   inflammatory<-t.test(plt_myeloid_fetal_random$inflammatory_macs_score1, plt_myeloid_ped$inflammatory_macs_score1)$p.value
#   recruit<-t.test(plt_myeloid_fetal_random$recently_recruited_myeloid1, plt_myeloid_ped$recently_recruited_myeloid1)$p.value
#   kuffer<-t.test(plt_myeloid_fetal_random$kuffer_like_score1, plt_myeloid_ped$kuffer_like_score1)$p.value
# 
#   data.frame(supressive=supressive, inflammatory=inflammatory, recruit=recruit ,kuffer=kuffer)
# }))
# 
# print(paste("Comparing the myeloid immune supressive score in myeloid cells, there is a sig difference between ped and adult (p <", pval,
#             ") in ", samp_num, " random samples at a p value of ",
#             round((length(myeloid_pval_montecarlo$supressive[which(myeloid_pval_montecarlo$supressive>pval)])+1)/(samp_num+1), 3), sep=""))
# print(paste("Comparing the inflammatory macs score in myeloid cells, there is a sig difference between ped and adult (p <", pval,
#             ") in ", samp_num, " random samples at a p value of ",
#             round((length(myeloid_pval_montecarlo$inflammatory[which(myeloid_pval_montecarlo$inflammatory>pval)])+1)/(samp_num+1), 3), sep=""))
# print(paste("Comparing the recently recruited myeloid score in myeloid cells, there is a sig difference between ped and adult (p <", pval,
#             ") in ", samp_num, " random samples at a p value of ",
#             round((length(myeloid_pval_montecarlo$recruit[which(myeloid_pval_montecarlo$recruit>pval)])+1)/(samp_num+1), 6), sep=""))
# print(paste("Comparing the kuffer-like score in myeloid cells, there is a sig difference between ped and adult (p <", pval,
#             ") in ", samp_num, " random samples at a p value of ",
#             round((length(myeloid_pval_montecarlo$kuffer[which(myeloid_pval_montecarlo$kuffer>pval)])+1)/(samp_num+1), 6), sep=""))
# 
# 
# 
# # BOX PLOT
# plt_not_gene_myeloid$age_condition<-factor(plt_not_gene_myeloid$age_condition, levels=c("Fetal Healthy","Ped Healthy","Ped IFALD","Adult Healthy" ))
# 
# cell_num_myeloid<-as.data.frame(table(plt_not_gene_myeloid$age_condition))
# colnames(cell_num_myeloid)<-c("age_condition","CellCount")
# 
# plt_max<-ceiling(max(plt_not_gene_myeloid$recently_recruited_myeloid1))
# plt_min<-floor(min(plt_not_gene_myeloid$recently_recruited_myeloid1))+0.5
# 
# myeloid_recruit_box<-
#   ggplot(plt_not_gene_myeloid, aes(age_condition,recently_recruited_myeloid1))+
#   geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=age_condition))+
#   theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#   geom_signif(stat="identity",
#               data=data.frame(x=c(1, 1, 2), xend=c(1, 2, 2),
#                               y=c(plt_max, plt_max+0.25, plt_max+0.25), yend=c(plt_max+0.25, plt_max+0.25, plt_max),
#                               annotation=c("*")),
#               aes(x=x,xend=xend, y=y, yend=yend, annotation=annotation), color="grey50")+ylim(plt_min, plt_max+1)+
#   xlab("Age Group")+ylab("Recently Recruited Myeloid Signature Score")+
#   geom_text(aes(y=-1, x=age_condition,label=paste0("n = ",comma(CellCount))),cell_num_myeloid, hjust=-0.1, size=3)
# myeloid_recruit_box
# 
# 
# plt_max<-ceiling(max(plt_not_gene_myeloid$kuffer_like_score1))
# plt_min<-floor(min(plt_not_gene_myeloid$kuffer_like_score1))+0.5
# 
# myeloid_kc_box<-
#   ggplot(plt_not_gene_myeloid, aes(age_condition,kuffer_like_score1))+
#   geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=age_condition))+
#   theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#   xlab("Age Group")+ylab("Kupffer-like Signature Score")+
#   geom_text(aes(y=-1, x=age_condition,label=paste0("n = ",comma(CellCount))),cell_num_myeloid, hjust=-0.1, size=3)
# myeloid_kc_box
# 
# 
# plt_max<-ceiling(max(plt_not_gene_myeloid$inflammatory_macs_score1))
# plt_min<-floor(min(plt_not_gene_myeloid$inflammatory_macs_score1))+0.5
# 
# myeloid_MHCII_box<-
#   ggplot(plt_not_gene_myeloid, aes(age_condition,inflammatory_macs_score1))+
#   geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=age_condition))+
#   theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#   xlab("Age Group")+ylab("MHCII Myeloid Signature Score")+
#   geom_text(aes(y=-1, x=age_condition,label=paste0("n = ",comma(CellCount))),cell_num_myeloid, hjust=-0.1, size=3)
# myeloid_MHCII_box
# 
# 
# plt_not_gene_myeloid$Age<-factor(plt_not_gene_myeloid$Age,
#                                  levels=c( "7 weeks gestation","8 weeks gestation", "9 weeks gestation",
#                                            "11 weeks gestation","12 weeks gestation", "13 weeks gestation",
#                                           "14 weeks gestation","16 weeks gestation","17 weeks gestation",
#                                           "2","3","11","12","13","17","26","48","57","61","65","67","69"))
# levels(plt_not_gene_myeloid$Age)<-c( "7\nweeks gestation","8\nweeks gestation", "9\nweeks gestation",
#                                      "11\nweeks gestation","12\nweeks gestation", "13\nweeks gestation",
#                                      "14\nweeks gestation","16\nweeks gestation","17\nweeks gestation",
#                                      "2","3","11","12","13","17","26","48","57","61","65","67","69")
# 
#   ggplot(plt_not_gene_myeloid, aes(Age,inflammatory_macs_score1))+
#   geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1)+
#   theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#   xlab("Age")
# 
#   ggplot(plt_not_gene_myeloid, aes(Age,kuffer_like_score1))+
#     geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1)+
#     theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#     xlab("Age")
# 
#   ggplot(plt_not_gene_myeloid, aes(Age,recently_recruited_myeloid1))+
#     geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1)+
#     theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#     xlab("Age")
# 

