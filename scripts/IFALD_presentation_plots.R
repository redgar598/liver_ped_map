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




source("scripts/00_pretty_plots.R")
source("scripts/00_entropy_d10x.R")


####################################
## Fancy dot plot
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
MacrophageCLEC9A_genes<-c("IDO1","CLEC9A")
platelet_genes<-c("PPBP","NRGN")


gene_exp<-FetchData(d10x, vars=c(Macrophage_genes,LEC_genes,Hepatocyte_genes,Cholangiocytes_genes,HSCs_genes,T_genes,NK_genes,gd_genes,RBC,
                                 MAST, recent_recruit_myeloid, kuffer_signature, neutro_gene, MHCII, b_genes_noIG, immunoglobins, 
                                 MacrophageCLEC9A_genes,platelet_genes))
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)

plt<-merge(gene_exp, d10x@meta.data, by.x="cell", by.y="index")


## summarize
scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}

plt_summary<-plt %>% group_by(CellType_refined, variable) %>% 
  summarise(mn=mean(value), count=length(which(value>0)), percent_exp=(length(which(value>0))/length(value))*100)
plt_summary <- plt_summary %>% group_by(variable) %>%
  dplyr::mutate(scaled = scale_this(mn))
plt_summary<-as.data.frame(plt_summary)

# remove dots where 0 cell expressing marker
plt_summary<-plt_summary[(which(plt_summary$count>0)),]

plt_summary$variable<-factor(plt_summary$variable, levels=rev(c(b_genes_noIG, immunoglobins, T_genes,gd_genes,NK_genes, Cholangiocytes_genes,LEC_genes,
                                                            Macrophage_genes,recent_recruit_myeloid,MHCII, kuffer_signature, MacrophageCLEC9A_genes, 
                                                            neutro_gene, MAST,platelet_genes, 
                                                            RBC,HSCs_genes,Hepatocyte_genes)))

fancy_dotplot<-plot_grid(
  ggplot(plt_summary, aes(CellType_refined, variable, color=scaled, size=percent_exp))+geom_point()+
    th+theme_classic()+
    scale_color_continuous_sequential(palette = "Oslo", rev=F, name="Scaled\nMean\nExpression")+
    scale_size(name="Percent\nCells\nExpressing")+
    theme(axis.text.x = element_blank(),axis.title = element_blank(),axis.ticks.x = element_blank(),
          plot.margin = margin(0.25,0.25,0,0.25,"cm")),
  ggplot(plt_summary, aes(CellType_refined, y=1, fill=CellType_refined))+geom_tile(color="black")+
    th+theme_classic()+fillscale_cellType+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",axis.line  = element_blank(),
          plot.margin = margin(0,0,1,1,"cm")),
  ncol=1, rel_heights = c(6,1.3), align = "v", axis="lr")
save_plts(fancy_dotplot, "IFALD_dot_plot_celltype", w=6,h=10)
