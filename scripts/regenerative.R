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


source("scripts/00_pretty_plots.R")
source("scripts/00_plot_gene_exp.R")
source("scripts/00_fanciest_UMAP.R")

load(here("data","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

#################
## HSC
#################
d10x.combined_hsc<-subset(d10x.combined, subset = CellType_rough %in% c("HSC"))
rm(d10x.combined)
gc()
d10x.combined_hsc<-subset(d10x.combined_hsc, subset = Treatment %in% c("Healthy"))

d10x.combined_hsc <- RunPCA(d10x.combined_hsc, npcs = 30, verbose = FALSE)
d10x.combined_hsc <- RunUMAP(d10x.combined_hsc, reduction = "pca", dims = 1:30)


HSC_differential_violin<-plot_grid(plot_gene_violin(d10x.combined_hsc,"MKI67"),
                                   plot_gene_violin(d10x.combined_hsc,"TOP2A"),
                                   plot_gene_violin(d10x.combined_hsc,"SOX4"))
HSC_differential_violin

HSC_differential_umap<-plot_grid(plot_gene_UMAP(d10x.combined_hsc,"MKI67",0),
                                 plot_gene_UMAP(d10x.combined_hsc,"TOP2A",0),
                                 plot_gene_UMAP(d10x.combined_hsc,"SOX4",0))
HSC_differential_umap


#################
## LSEC
#################
d10x.combined_lsec<-subset(d10x.combined, subset = CellType_rough %in% c("LSEC"))
rm(d10x.combined)
gc()
d10x.combined_lsec<-subset(d10x.combined_lsec, subset = Treatment %in% c("Healthy"))

d10x.combined_lsec <- RunPCA(d10x.combined_lsec, npcs = 30, verbose = FALSE)
d10x.combined_lsec <- RunUMAP(d10x.combined_lsec, reduction = "pca", dims = 1:30)


LSEC_differential_violin<-plot_grid(plot_gene_violin(d10x.combined_lsec,"MKI67"),
                                   plot_gene_violin(d10x.combined_lsec,"TOP2A"),
                                   plot_gene_violin(d10x.combined_lsec,"SOX4"),
                                   plot_gene_violin(d10x.combined_lsec,"VEGFA"))
LSEC_differential_violin

LSEC_differential_umap<-plot_grid(plot_gene_UMAP(d10x.combined_lsec,"MKI67",0),
                                 plot_gene_UMAP(d10x.combined_lsec,"TOP2A",0),
                                 plot_gene_UMAP(d10x.combined_lsec,"SOX4",0),
                                 plot_gene_UMAP(d10x.combined_lsec,"VEGFA",0))
LSEC_differential_umap




#################
## cholangiocytes
#################
d10x.combined_cholangiocytes<-subset(d10x.combined, subset = CellType_rough %in% c("Cholangiocytes"))
rm(d10x.combined)
gc()
d10x.combined_cholangiocytes<-subset(d10x.combined_cholangiocytes, subset = Treatment %in% c("Healthy"))

d10x.combined_cholangiocytes <- RunPCA(d10x.combined_cholangiocytes, npcs = 30, verbose = FALSE)
d10x.combined_cholangiocytes <- RunUMAP(d10x.combined_cholangiocytes, reduction = "pca", dims = 1:30)


cholangiocytes_differential_violin<-plot_grid(plot_gene_violin(d10x.combined_cholangiocytes,"MKI67"),
                                    plot_gene_violin(d10x.combined_cholangiocytes,"TOP2A"),
                                    plot_gene_violin(d10x.combined_cholangiocytes,"SOX4"),
                                    plot_gene_violin(d10x.combined_cholangiocytes,"VEGFA"))
cholangiocytes_differential_violin

cholangiocytes_differential_umap<-plot_grid(plot_gene_UMAP(d10x.combined_cholangiocytes,"MKI67",0),
                                  plot_gene_UMAP(d10x.combined_cholangiocytes,"TOP2A",0),
                                  plot_gene_UMAP(d10x.combined_cholangiocytes,"SOX4",0),
                                  plot_gene_UMAP(d10x.combined_cholangiocytes,"VEGFA",0))
cholangiocytes_differential_umap

#################
## hepatocytes
#################
d10x.combined_hep<-subset(d10x.combined, subset = CellType_rough %in% c("Hepatocytes"))
rm(d10x.combined)
gc()
d10x.combined_hep<-subset(d10x.combined_hep, subset = Treatment %in% c("Healthy"))

d10x.combined_hep <- RunPCA(d10x.combined_hep, npcs = 30, verbose = FALSE)
d10x.combined_hep <- RunUMAP(d10x.combined_hep, reduction = "pca", dims = 1:30)


hep_differential_violin<-plot_grid(plot_gene_violin(d10x.combined_hep,"MKI67"),
                                              plot_gene_violin(d10x.combined_hep,"TOP2A"),
                                              plot_gene_violin(d10x.combined_hep,"SOX4"),
                                              plot_gene_violin(d10x.combined_hep,"VEGFA"))
hep_differential_violin

hep_differential_umap<-plot_grid(plot_gene_UMAP(d10x.combined_hep,"MKI67",0),
                                            plot_gene_UMAP(d10x.combined_hep,"TOP2A",0),
                                            plot_gene_UMAP(d10x.combined_hep,"SOX4",0),
                                            plot_gene_UMAP(d10x.combined_hep,"VEGFA",0))
hep_differential_umap
