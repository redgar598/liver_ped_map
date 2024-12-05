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


source("scripts/00_pretty_plots.R")
source("scripts/00_fanciest_UMAP.R")
source("scripts/00_plot_gene_exp.R")
source("scripts/00_entropy_d10x.R")



load(here("data","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

d10x.combined_hepatocytes<-subset(d10x.combined, subset = CellType_refined %in% c("Hepatocytes"))
rm(d10x.combined)
gc()
d10x.combined_hepatocytes <- RunPCA(d10x.combined_hepatocytes, npcs = 30, verbose = FALSE)
d10x.combined_hepatocytes <- RunUMAP(d10x.combined_hepatocytes, reduction = "pca", dims = 1:30)
d10x.combined_hepatocytes <- FindNeighbors(d10x.combined_hepatocytes, reduction = "pca", dims = 1:30)
d10x.combined_hepatocytes <- FindClusters(d10x.combined_hepatocytes, resolution = 0.3)

hepatocyte_subtype<-DimPlot(d10x.combined_hepatocytes, label=T)
hepatocyte_subtype
save_plts(hepatocyte_subtype, "IFALD_hepatocyte_map_clusters", w=7,h=6)




zonation_markers_PC<-c("CYP1A1","CYP3A4","CYP2E1","GLUL","AXIN2","TBX3")
zonation_markers_PP<-c("ASS1","ASL","ALB","CYP2F1","ARG1","CPS1","C2","MFSD2A")

zonation_markers_PC_key<-c("CYP1A1","CYP3A4","CYP2E1","GLUL")
zonation_markers_PP_key<-c("ASS1","ASL","ALB","SERPINA1")

FeaturePlot(d10x.combined_hepatocytes, features = zonation_markers_PC_key, min.cutoff = "q9", pt.size=1)
FeaturePlot(d10x.combined_hepatocytes, features = zonation_markers_PP, min.cutoff = "q9", pt.size=1)


hepatocytes_cluster_umap_cellstatus<-DimPlot(d10x.combined_hepatocytes, reduction = "umap", pt.size=0.25, group.by="cell_status",label=T)+scale_color_manual(values=c("grey","red","cornflowerblue"),name="Nuclear\nFraction")
hepatocytes_cluster_umap_cellstatus
save_plts(hepatocytes_cluster_umap_cellstatus, "IFALD_hepatocytes_cluster_umap_cellstatus", w=6,h=4)
umap_MThepatocytes<-FeaturePlot(d10x.combined_hepatocytes, features = "percent.mt", min.cutoff = "q9", pt.size=1)
save_plts(umap_MThepatocytes, "IFALD_umap_MThepatocytes", w=5,h=4)

phase_hepatocytes<-DimPlot(d10x.combined_hepatocytes, reduction = "umap", group.by = "Phase" )+scale_color_manual(values=c("#e6ab02","#386cb0","#1b9e77"))
save_plts(phase_hepatocytes, "IFALD_umap_phase_hepatocytes", w=5,h=4)

phase_marker<-FeaturePlot(d10x.combined_hepatocytes, reduction = "umap", features = c("MKI67","TOP2A"), ncol = 2)
save_plts(phase_marker, "IFALD_umap_phaseMarker_hepatocytes", w=10,h=4)

## see if one cluster is anythign specific
cluster4.markers <-  FindMarkers(d10x.combined_hepatocytes, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 10)
FeaturePlot(d10x.combined_hepatocytes, reduction = "umap", features = c("FABP1","MT2A","PEBP1","FTL"), ncol = 2)
plot_grid(plot_gene_UMAP(d10x.combined_hepatocytes,"FABP1", 0),
          plot_gene_UMAP(d10x.combined_hepatocytes,"MT2A", 0),
          plot_gene_UMAP(d10x.combined_hepatocytes,"PEBP1", 0),
          plot_gene_UMAP(d10x.combined_hepatocytes,"FTL", 0), ncol=4)

Macrophage_genes<-c( "PTPRC", "CD68", "MARCO","CD5L","VSIG4", "MAF", "LYZ", "CSTA", "S100A8", "S10049",
                     "CD14", "CD74", "GPBAR1", "ID3")

B_genes<-c("POU2F2","FCER2","MS4A1","LTB","CD37","CD79B","IGLC2","IGHG1","IGKC", "CD19")
NK_T_genes<-c("CD3D","TRDC","SELL","IL7R","CCR7","S100A4","CD8A","GNLY","NKG7")


LEC_genes<-c("CALCRL", "VWF", "RAMP2", "STAB2", "LYVE1", "PECAM1", "ENG", "FCGR2B", "F8", "SPARCL1",
             "ID1", "SOX18", "CD32B", "ID3")
Hepatocyte_genes<-c("ALB", "HAMP", "ARG1", "PCK1", "AFP", "BCHE", "HAL", "SCD", "CPS1", "CYP3A4",
                    "ELF3", "CRP", "GSTA2", "AKR1C1", "MGST1", "CYP3A5", "ALDH1A1", "ADH1A", "CYP2E1",
                    "GLS2", "SDS", "GLUL", "AKR1D1", "HPR",
                    "HMGCS1", "IGSF23", "ACSS2", "G6PC", "ID3")
Cholangiocytes_genes<-c( "EPCAM", "SOX9", "KRT1", "KRT7", "ANXA4", "KRT18", "ID3")

HSCs_genes<-c( "RBP1", "LRAT", "PDE3B", "ACTA2", "AOX1", "PDE3D", "PDE4D", "SPARC", "TAGLN", "COL1A1", "COL1A2", "COL3A1",
               "TIMP1", "DCN", "MYL9", "TPM2", "MEG3", "BGN", "IGFBP7", "IGFBP3", "CYR61", "IGFBP6", "CCL2", "COLEC11",
               "CTGF", "HGF", "ID3")

T_genes<-c("CD3D","IL7R","CD8A","IL32")
NK_genes<-c("NKG7","CD7")
gd_genes<-c("GNLY")
## Dot plot of all markers
DefaultAssay(d10x.combined_hepatocytes) <- "RNA"
pdf(file = here("figures/IFALD_dot_plots_hepatocytes.pdf"), w=10, h=10)
DotPlot(object = d10x.combined_hepatocytes, features = B_genes)+xlab("B Cell Marker")
DotPlot(object = d10x.combined_hepatocytes, features = T_genes)+xlab("T Cell Marker")
DotPlot(object = d10x.combined_hepatocytes, features = NK_genes)+xlab("NK Cell Marker")
DotPlot(object = d10x.combined_hepatocytes, features = LEC_genes)+xlab("LSEC Marker")
DotPlot(object = d10x.combined_hepatocytes, features = zonation_markers_PC)+xlab("Hepatocyte PC Marker")
DotPlot(object = d10x.combined_hepatocytes, features = zonation_markers_PP)+xlab("Hepatocyte PP Marker")
DotPlot(object = d10x.combined_hepatocytes, features = Cholangiocytes_genes)+xlab("Cholangiocyte Marker")
DotPlot(object = d10x.combined_hepatocytes, features = HSCs_genes[1:14])+xlab("HSC Marker")
DotPlot(object = d10x.combined_hepatocytes, features = HSCs_genes[15:27])+xlab("HSC Marker")
DotPlot(object = d10x.combined_hepatocytes, features = Macrophage_genes)+xlab("Macrophage Marker")
dev.off()
DefaultAssay(d10x.combined_hepatocytes) <- "integrated"







##############
## fibrosis
##############
fibrosis_markers<-plot_grid(plot_gene_UMAP(d10x.combined_hepatocytes,"PDGFRA", 0),
                       plot_gene_UMAP(d10x.combined_hepatocytes,"CXCL12", 0),
                       plot_gene_UMAP(d10x.combined_hepatocytes,"COL1A1", 0),
                       plot_gene_UMAP(d10x.combined_hepatocytes,"IGFBP3", 0), ncol=4)
fibrosis_markers


proliferation_markers<-plot_grid(plot_gene_UMAP(d10x.combined_hepatocytes,"TOP2A", 0),
                            plot_gene_UMAP(d10x.combined_hepatocytes,"MKI67", 0), ncol=2)
proliferation_markers







