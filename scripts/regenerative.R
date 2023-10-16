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




####################################
## Fancy dot plot of gene s from Michael/Olivia
####################################
d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

## add cell type labels
load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

## markers
regenerative<-c( "RELA", "YAP1", "STAT3", "MYC", "JUN", "FOS", "IGF2", "TGFA", "VEGFA")
proliferative<-c("MKI67","TOP2A", "CCNA2", "CDK1", "CENPM")


gene_exp<-FetchData(d10x, vars=c(regenerative,proliferative))
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)

plt<-merge(gene_exp, d10x@meta.data, by.x="cell", by.y="index")


## summarize
scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}

plt_summary<-plt %>% group_by(CellType_refined,age_condition, variable) %>% 
  summarise(mn=mean(value), count=length(which(value>0)), percent_exp=(length(which(value>0))/length(value))*100)
plt_summary <- plt_summary %>% group_by(variable) %>%
  dplyr::mutate(scaled = scale_this(mn))
plt_summary<-as.data.frame(plt_summary)

# remove dots where 0 cell expressing marker
plt_summary<-plt_summary[(which(plt_summary$count>0)),]

plt_summary$cell_condition<-paste(plt_summary$CellType_refined, plt_summary$age_condition)

plt_summary$variable<-factor(plt_summary$variable, levels=rev(c(regenerative,proliferative)))

fancy_dotplot<-plot_grid(
  ggplot(plt_summary, aes(cell_condition, variable, color=scaled, size=percent_exp))+geom_point()+
    th+theme_classic()+
    scale_color_continuous_sequential(palette = "Oslo", rev=F, name="Scaled\nMean\nExpression")+
    scale_size(name="Percent\nCells\nExpressing")+
    theme(axis.text.x = element_blank(),axis.title = element_blank(),axis.ticks.x = element_blank(),
          plot.margin = margin(0.25,0.25,0,0.25,"cm"))+
    geom_vline(xintercept = c(seq(3.5,48.5, 3),seq(50.5,70, 3))),
  ggplot(plt_summary, aes(cell_condition, y=1, fill=CellType_refined))+geom_tile(color="black")+
    th+theme_classic()+fillscale_cellType+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",axis.line  = element_blank(),
          plot.margin = margin(0,0,1,1,"cm")),
  ncol=1, rel_heights = c(4,1.3), align = "v", axis="lr")
fancy_dotplot



#############
## Heat Map of regenerative genes
#############
d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))

load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

d10x.combined_cholangiocytesHSC<-subset(d10x, subset = CellType_refined %in% c("Cholangiocytes","HSC","LSEC"))

sig_de<-read.csv(file=here("data","differential_age_Cholangiocytes.csv"))
sig_de$age_condition<-"Adult\nHealthy"
sig_de$CellType_refined<-"Cholangiocytes"

sig_de_IFALD<-read.csv(file=here("data","differential_IFALD_Cholangiocytes.csv"))
sig_de_IFALD$age_condition<-"Ped\nIFALD"
sig_de_IFALD$CellType_refined<-"Cholangiocytes"

sig_de_hsc<-read.csv(file=here("data","differential_ageHSC.csv"))
sig_de_hsc$age_condition<-"Adult\nHealthy"
sig_de_hsc$CellType_refined<-"HSC"

sig_de_IFALD_hsc<-read.csv(file=here("data","differential_IFALDHSC.csv"))
sig_de_IFALD_hsc$age_condition<-"Ped\nIFALD"
sig_de_IFALD_hsc$CellType_refined<-"HSC"

sig_de_lsec<-read.csv(file=here("data","differential_age_lsec.csv"))
sig_de_lsec$age_condition<-"Adult\nHealthy"
sig_de_lsec$CellType_refined<-"LSEC"

sig_de_IFALD_lsec<-read.csv(file=here("data","differential_IFALD_lsec.csv"))
sig_de_IFALD_lsec$age_condition<-"Ped\nIFALD"
sig_de_IFALD_lsec$CellType_refined<-"LSEC"

de<-rbind(sig_de, sig_de_IFALD, sig_de_hsc, sig_de_IFALD_hsc,sig_de_lsec,sig_de_IFALD_lsec)
colnames(de)[1]<-"variable"
de$label<-"*"

cholangiocytes_age_heat<-plot_heat_map(d10x.combined_cholangiocytesHSC,regenerative, 
                                c("Cholangiocytes","HSC","LSEC"),T)
cholangiocytes_age_heat
save_plts(cholangiocytes_age_heat, "cholangiocyte_LSEC_HSC_age_heat", w=7,h=5)



