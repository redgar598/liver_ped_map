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


source("scripts/00_pretty_plots.R")


load(here("data","adult_ped_integrated_refinedlabels_withDropletQC.rds"))

d10x.combined_myeloid<-subset(d10x.combined, subset = CellType_rough %in% c("Myeloid cells"))
rm(d10x.combined)
gc()
d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)
d10x.combined_myeloid <- RunUMAP(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindNeighbors(d10x.combined_myeloid, reduction = "pca", dims = 1:30)
d10x.combined_myeloid <- FindClusters(d10x.combined_myeloid, resolution = 0.1)
DimPlot(d10x.combined_myeloid, label=T)

myeloid_cluster_umap<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=T, group.by = "CellType_refined")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-14, y=-12, label=paste0("n = ",comma(ncol(d10x.combined_myeloid))))
myeloid_cluster_umap
save_plts(myeloid_cluster_umap, "myeloid_cluster_umap_labelled_refined", w=7,h=5)

myeloid_cluster_umap<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=F, group.by = "CellType_refined")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-14, y=-12, label=paste0("n = ",comma(ncol(d10x.combined_myeloid))))
myeloid_cluster_umap
save_plts(myeloid_cluster_umap, "myeloid_cluster_umap_nolab_refined", w=7,h=5)

myeloid_cluster_umap<-DimPlot(d10x.combined_myeloid, reduction = "umap", pt.size=0.25, label=T, group.by = "CellType_refined")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-14, y=-12, label=paste0("n = ",comma(ncol(d10x.combined_myeloid))))
myeloid_cluster_umap


###################
## Some silly plots
###################
umap_mat<-as.data.frame(Embeddings(object = d10x.combined_myeloid, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x.combined_myeloid@meta.data
meta$cell<-rownames(meta)

plt<-merge(meta, umap_mat, by="cell")

myeloid_density<-ggplot(plt, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(aes(color=CellType_refined),alpha=0.1, shape=16, size=2) + 
  geom_density_2d(bins=6, color="grey30")+
  theme_bw()+colscale_cellType
save_plts(myeloid_density, "myeloid_density", w=6,h=4)

##### DO this on raw data
d10x<-readRDS(file = here("data","d10x_adult_ped_raw.rds"))
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

myeloid_markers<-FetchData(object = d10x, vars = c("MARCO", "VSIG4","LYZ"))
myeloid_markers$cell<-rownames(myeloid_markers)
plt<-merge(plt, myeloid_markers, by="cell")

plt_nonzero<-plt[which(plt$MARCO>0 & plt$LYZ>0),]

myeloid_density_markers<-ggplot()+geom_point(aes(MARCO, LYZ, fill=CellType_refined),plt,  shape=21, colour = "transparent", size=1)+  theme_bw()+fillscale_cellType+
  geom_density_2d(aes(MARCO, LYZ, color=CellType_refined),plt_nonzero, bins=6, color="grey30")+
  guides(fill = guide_legend(override.aes = list(size=2)))
save_plts(myeloid_density_markers, "myeloid_density_markers", w=6,h=4)



##############
## fibrotic or regenerative macrophages
##############
FeaturePlot(d10x.combined_myeloid, features = c("CD86","MRC1"), min.cutoff = "q9", pt.size=1)


            # myeliod_clusters<-d10x.combined_myeloid@meta.data
            # 
            # ## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
            # # but not normalized or scaled
            # d10x<-readRDS(file = here("data","d10x_adult_ped_raw.rds"))
            # d10x_raw_myeloid<-subset(d10x, cells = rownames(myeliod_clusters))
            # rm(d10x)
            # gc()
            # 
            # myeliod_clusters$index<-rownames(myeliod_clusters)
            # identical(colnames(d10x_raw_myeloid), myeliod_clusters$index)
            # 
            # d10x_raw_myeloid <- AddMetaData(d10x_raw_myeloid, metadata = myeliod_clusters)
            # 
            # d10x_raw_myeloid <- NormalizeData(d10x_raw_myeloid,scale.factor = 10000, normalization.method = "LogNormalize")
            # 
            # 
            # 
            # ## testing factor
            # Idents(d10x_raw_myeloid) <- "seurat_clusters"
            # table(d10x_raw_myeloid$seurat_clusters)
            # 
            # 
            # #MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene,
            # #consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
            # #and linear regression for the continuous process (i.e., the expression level). 
            # 
            # myeloid_clusters<-unique(d10x_raw_myeloid$seurat_clusters)
            # 
            # diff_exp_all<-lapply(1:length(myeloid_clusters), function(x){
            #   de<-FindMarkers(d10x_raw_myeloid, ident.1 = myeloid_clusters[x], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
            #   de$gene<-rownames(de)
            #   rownames(de)<-NULL
            #   de<-de[,c(6,1:5)]
            #   de$cluster<-myeloid_clusters[x]
            #   de})
            # 
            # 
            # diff_exp_all<-do.call(rbind, diff_exp_all)
            # 
            # diff_exp_sig<-diff_exp_all[which(diff_exp_all$p_val_adj<0.005),]
            # 
            # 
            # 
            # # ### Top DE genes
            # diff_exp_sig %>%
            #   group_by(cluster) %>%
            #   top_n(n = 10, wt = abs(avg_log2FC)) -> top10
            # 
            # top_DE<-as.data.frame(top10)
            # 
            # diff_exp_sig %>%
            #   group_by(cluster) %>%
            #   top_n(n = 10, wt = avg_log2FC) -> top_DE_up
            # 
            # top_DE_up<-as.data.frame(top_DE_up)
            # 
            # ### plot markers
            # FeaturePlot(d10x.combined_myeloid, features = "DEFA3", min.cutoff = "q9", pt.size=1)
            # 
            # VlnPlot(d10x.combined_myeloid, features = "DEFA3", pt.size = 0, log=T)
            # VlnPlot(d10x.combined_myeloid, features = "DEFA3", pt.size = 0, log=F)
            # 
            # recent_recruit_myeloid<-c("S100A8","S100A9","CD68","LYZ")
            # kuffer_signature<-c("VSIG4","MARCO","CD5L","HMOX1")
            # 
            # diff_exp_sig[which(diff_exp_sig$gene%in%recent_recruit_myeloid),]
            # diff_exp_sig[which(diff_exp_sig$gene%in%kuffer_signature),]
            # 
            # FeaturePlot(d10x.combined_myeloid, features = kuffer_signature, min.cutoff = "q9", pt.size=1)
            # FeaturePlot(d10x.combined_myeloid, features = recent_recruit_myeloid, min.cutoff = "q9", pt.size=1)
            # 
            # 
            # ###################################################
            # #0: recently recruited myeloid (monocyte)
            # top_DE_up[which(top_DE_up$cluster=="0"),]
            # top_DE[which(top_DE$cluster=="0"),]
            # myeloid_0<-FeaturePlot(d10x.combined_myeloid, features = c("LYZ","S100A10","AHNAK","VCAN"), min.cutoff = "q9", pt.size=1)
            # myeloid_0
            # save_plts(myeloid_0, "myeloid_cluster_0_markers", w=8,h=6)
            # 
            # 
            # #1: recently recruited myeloid (monocyte)
            # top_DE[which(top_DE$cluster=="1"),]
            # top_DE_up[which(top_DE_up$cluster=="1"),]
            # myeloid_1<-FeaturePlot(d10x.combined_myeloid, features = c("AREG","GPR183","IFI30","RPL17"), min.cutoff = "q9", pt.size=1)
            # myeloid_1
            # save_plts(myeloid_1, "myeloid_cluster_1_markers", w=8,h=6)
            # 
            # 
            # # 2,4,5 kupffer-like
            # top_DE_up[which(top_DE_up$gene%in%kuffer_signature),]
            # myeloid_245<-FeaturePlot(d10x.combined_myeloid, features = kuffer_signature, min.cutoff = "q9", pt.size=1)
            # myeloid_245
            # save_plts(myeloid_245, "myeloid_cluster_245_markers", w=8,h=6)
            # 
            # 
            # #7: unclear soupy?
            # top_DE[which(top_DE$cluster=="7"),]
            # top_DE_up[which(top_DE_up$cluster=="7"),]
            # myeloid_7<-FeaturePlot(d10x.combined_myeloid, features = c("MT1M","APOC1","MT1G","SAA2"), min.cutoff = "q9", pt.size=1)
            # myeloid_7
            # save_plts(myeloid_7, "myeloid_cluster_7_markers", w=8,h=6)
            # 
            # 
            # #6: Neutrophils
            # #neutrophil makers from literature
            # neutro_gene<-c("CSF3R","FCGR3B","NAMPT","CXCR2")
            # 
            # top_DE[which(top_DE$cluster=="6"),]
            # top_DE_up[which(top_DE_up$cluster=="6"),]
            # myeloid_6<-FeaturePlot(d10x.combined_myeloid, features = c("AZU1","DEFA3","DEFA4","ELANE"), min.cutoff = "q9", pt.size=1)
            # myeloid_6
            # save_plts(myeloid_6, "myeloid_cluster_6_markers", w=8,h=6)
            # FeaturePlot(d10x.combined_myeloid, features = neutro_gene, min.cutoff = "q9", pt.size=1)
            # diff_exp_sig[which(diff_exp_sig$gene%in%neutro_gene),]
            # 
            # 
            # #3: Neutrophils
            # top_DE[which(top_DE$cluster=="3"),]
            # top_DE_up[which(top_DE_up$cluster=="3"),]
            # FeaturePlot(d10x.combined_myeloid, features = c("DEFA3","DEFA4","OLFM4","LTF"), min.cutoff = "q9", pt.size=1)
            # 
            # diff_exp_sig[which(diff_exp_sig$gene%in%neutro_gene),]
            # myeloid_3<-FeaturePlot(d10x.combined_myeloid, features = neutro_gene, min.cutoff = "q9", pt.size=1)
            # myeloid_3
            # save_plts(myeloid_3, "myeloid_cluster_3_markers", w=8,h=6)
            # 

######
## Neutrophils have few genes expressed?
######
umap_nfeaturemyeloid<-FeaturePlot(d10x.combined_myeloid, features = "nFeature_RNA", min.cutoff = "q9", pt.size=1)
umap_nfeaturemyeloid + scale_color_continuous_sequential(palette = "ag_GrnYl") 

umap_ncountmyeloid<-FeaturePlot(d10x.combined_myeloid, features = "nCount_RNA", min.cutoff = "q9", pt.size=1)
umap_ncountmyeloid + scale_color_continuous_sequential(palette = "ag_GrnYl") 

#### Add back in cell cycle cause missing
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
d10x_raw_myeloid <- CellCycleScoring(d10x_raw_myeloid, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
identical(colnames(d10x_raw_myeloid), colnames(d10x.combined_myeloid))
d10x.combined_myeloid <- AddMetaData(d10x.combined_myeloid, metadata = d10x_raw_myeloid$Phase, col.name = 'Phase')


myeloid_cycle<-DimPlot(d10x.combined_myeloid, group.by = "Phase",pt.size=1) + scale_color_manual(values=c("#e6ab02","#386cb0","#1b9e77"))
myeloid_cycle
save_plts(myeloid_cycle, "myeloid_cell_cycle", w=8,h=6)



## RBC stuff
# Among the cognate receptors for these senescence signals that are expressed on RPMs and KCs, 
# TAM receptors AXL and MERTK, TIM4 and Fc receptor CD16 are the most abundant. A “don’t eat me” 
#signal that prevents the clearance of young and intact erythrocytes is provided by the interaction 
# between SIRPα at the surface of the macrophage with CD47 expressed by erythrocytes.
RBC_receptors<-FeaturePlot(d10x.combined_myeloid, features = c("TREM2","CD9","AXL","CD47","MERTK","TIMD4","FCGR3A"), min.cutoff = "q9", pt.size=1)
RBC_receptors
save_plts(RBC_receptors, "RBC_receptors", w=8,h=6)

FeaturePlot(d10x.combined_myeloid, features = c("HBA2", "AXL"), blend = TRUE)
RBC_Overlap<-FeaturePlot(d10x.combined_myeloid, features = c("HBA2", "FCGR3A"), blend = TRUE)
RBC_Overlap
save_plts(RBC_Overlap, "RBC_hem_receptor_Overlap", w=14,h=3)

                            # 
                            # d10x.combined_kupffer<-subset(d10x.combined_myeloid, subset = seurat_clusters %in% c(3,4,7))
                            # rm(d10x.combined)
                            # gc()
                            # d10x.combined_kupffer <- RunPCA(d10x.combined_kupffer, npcs = 30, verbose = FALSE)
                            # d10x.combined_kupffer <- RunUMAP(d10x.combined_kupffer, reduction = "pca", dims = 1:30)
                            # d10x.combined_kupffer <- RunTSNE(d10x.combined_kupffer, reduction = "pca", dims = 1:30)
                            # d10x.combined_kupffer <- FindNeighbors(d10x.combined_kupffer, reduction = "pca", dims = 1:30)
                            # d10x.combined_kupffer <- FindClusters(d10x.combined_kupffer, resolution = 0.1)
                            # DimPlot(d10x.combined_kupffer, label=T)
                            # DimPlot(d10x.combined_kupffer, label=T, reduction="tsne")
                            # 
                            # RBC_receptors_KC<-FeaturePlot(d10x.combined_kupffer, features = c("TREM2","CD9","AXL","CD47","MERTK","TIMD4","FCGR3A"), min.cutoff = "q9", pt.size=1)
                            # RBC_receptors_KC
                            # RBC_hem_KC<-FeaturePlot(d10x.combined_kupffer, features = c("HBA2","HBA1","HBB","MARCO"), min.cutoff = "q9", pt.size=1)
                            # RBC_hem_KC
                            # RBC_Overlap<-FeaturePlot(d10x.combined_kupffer, features = c("HBA2", "FCGR3A"), blend = TRUE)
                            # RBC_Overlap



######
## Kupffer-like cells
######
# Guilliams KC markers
KC_markers<-FeaturePlot(d10x.combined_myeloid, features = c("SLC16A9","CD5L","SLC40A1","VSIG4"), min.cutoff = "q9", pt.size=1)
KC_markers

# Guilliams Bile-duct associated lipid associated macropahges (BC-LAM) used human orthologos
LAM_markers<-FeaturePlot(d10x.combined_myeloid, features = c("TREM2","CLEC4A","GPNMB","SPP1","TREM2","LGALS1","FTL"), min.cutoff = "q9", pt.size=1)
LAM_markers



###########
## DGE KC-like and recent recruited mono
###########
myeliod_clusters<-d10x.combined_myeloid@meta.data

## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","d10x_adult_ped_raw.rds"))
d10x_raw_myeloid<-subset(d10x, cells = rownames(myeliod_clusters))
rm(d10x)
gc()

myeliod_clusters$index<-rownames(myeliod_clusters)
identical(colnames(d10x_raw_myeloid), myeliod_clusters$index)

d10x_raw_myeloid <- AddMetaData(d10x_raw_myeloid, metadata = myeliod_clusters)

d10x_raw_myeloid <- NormalizeData(d10x_raw_myeloid,scale.factor = 10000, normalization.method = "LogNormalize")


Idents(d10x_raw_myeloid) <- "CellType_refined"

de<-FindMarkers(d10x_raw_myeloid, ident.1 = "RR Myeloid", ident.2 = "KC Like", test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]


#########
## pathway KC versus recently recruited
#########
source("scripts/00_GSEA_function.R")

GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

gene_list = de$avg_log2FC
names(gene_list) = rownames(de)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

FeaturePlot(d10x.combined_myeloid, features = c("VCAN","LYZ","CD5L","MARCO"), min.cutoff = "q9", pt.size=1)

res = GSEA(gene_list, GO_file, pval = 0.05)
top10_pathways<-data.frame(pathway=sapply(1:10, function(x) strsplit(res$Results$pathway[x], "%")[[1]][1]), direction=sapply(1:10, function(x) res$Results$Enrichment[x]))

plt_path<-res$Results
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Enrichment_Cell<-"Up-regulated in \nrecently recruited myeloid"
plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \nKupffer-like cells"

plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))

plt_path$direction_label<-as.factor(plt_path$Enrichment)
levels(plt_path$direction_label)<-c(0.1,-0.1)
plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))

GO_KC_vs_recentrcruti<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_vline(xintercept = 0, color="grey40")+scale_fill_manual(values=c("#fd8d3c","#6baed6"))+ 
  guides(fill = guide_legend(override.aes = list(size=5)))
save_plts(GO_KC_vs_recentrcruti, "GO_KC_vs_recently_recruited", w=20,h=10)




#########
## KC versus macrophage
#########
# in the liver fetal atlas the say 
#"Tissue-resident Kupffer cells could be distinguished from monocyte-derived macrophages 
# by the expression of MARCO, CD163, FCGR3A and CD5L and the absence of LSP1 and CD48 expression


KC_vs_macrophage<-FeaturePlot(d10x.combined_myeloid, features = c("MARCO","CD163","FCGR3A","CD5L","LSP1","CD48"), min.cutoff = "q9", pt.size=1)
KC_vs_macrophage
save_plts(KC_vs_macrophage, "KC_vs_macrophage_markers", w=8,h=10)

