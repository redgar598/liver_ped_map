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
source("scripts/00_fanciest_UMAP.R")
source("scripts/00_plot_gene_exp.R")

####################################
## Fancy dot plot
####################################
d10x<-readRDS(file = here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_d10x_adult_ped_raw.rds"))
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

## add cell type labels
load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

## markers
CBA_panels<-read.csv(here("data","CBA_Panels.csv"))


## only myeloids
d10x_raw_myeloid<-subset(d10x, subset = CellType_refined %in% c("Mono-Mac","Macrophage\n(MHCII high)","KC Like"))
rm(d10x)
gc()

## only healthy
d10x_raw_myeloid<-subset(d10x_raw_myeloid, subset = age_condition %in% c("Ped Healthy","Adult Healthy"))


gene_exp<-FetchData(d10x_raw_myeloid, vars=CBA_panels$Gene)
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)

plt<-merge(gene_exp, d10x_raw_myeloid@meta.data, by.x="cell", by.y="index")


## summarize
scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}

plt_summary<-plt %>% group_by(CellType_refined, variable, age_condition) %>% 
  summarise(mn=mean(value), count=length(which(value>0)), percent_exp=(length(which(value>0))/length(value))*100)
plt_summary <- plt_summary %>% group_by(variable) %>%
  dplyr::mutate(scaled = scale_this(mn))
plt_summary<-as.data.frame(plt_summary)

# remove dots where 0 cell expressing marker
plt_summary<-plt_summary[(which(plt_summary$count>0)),]
# 
# plt_summary$variable<-factor(plt_summary$variable, levels=rev(c(b_genes_noIG, immunoglobins, T_genes,gd_genes,NK_genes, Cholangiocytes_genes,LEC_genes,
#                                                                 Macrophage_genes,recent_recruit_myeloid,MHCII, kuffer_signature, neutro_gene, MAST, 
#                                                                 RBC,HSCs_genes,Hepatocyte_genes)))

## gene order
CBA_panels_summary<- CBA_panels %>% select(Gene, Panel)
CBA_panels_summary<-CBA_panels_summary[which(CBA_panels_summary$Gene%in%plt_summary$variable),]
order_genes<-as.data.frame(CBA_panels_summary %>% 
                             group_by(Gene) %>%
                             summarise(y = paste0(Panel, collapse="_")))
order_genes<-order_genes[rev(order(order_genes$y)),]
CBA_panels_summary$Gene<-factor(CBA_panels_summary$Gene, levels=order_genes$Gene)

plt_summary$variable<-factor(plt_summary$variable, levels=order_genes$Gene)

## differential 
sig_de_age_RR<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_RR.csv"))
sig_de_age_RR$CellType_refined<-"Mono-Mac"
sig_de_age_KC<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_KC.csv"))
sig_de_age_KC$CellType_refined<-"KC Like"
sig_de_age_MHCII<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_MHCII.csv"))
sig_de_age_MHCII$CellType_refined<-"Macrophage\n(MHCII high)"

differential<-rbind(sig_de_age_RR, sig_de_age_KC,sig_de_age_MHCII)
colnames(differential)[1]<-"variable"

differential<-differential[which(differential$variable%in%plt_summary$variable),]




fancy_dotplot<-plot_grid(
  ggplot(CBA_panels_summary, aes(Panel, Gene, fill=Panel))+geom_tile(color="black")+
    th+theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",axis.line  = element_blank(),
          plot.margin = margin(0,0,1,1,"cm")),
  ggplot()+geom_point(aes(age_condition, variable, color=scaled, size=percent_exp), plt_summary)+
    th+theme_bw()+facet_grid(~CellType_refined)+
    scale_color_continuous_sequential(palette = "Blues 3", rev=T, name="Scaled\nMean\nExpression")+
    scale_size(name="Percent\nCells\nExpressing")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.margin = margin(0.25,0.25,0,0.25,"cm"))+
    geom_hline(yintercept = c(seq(1.5,26.5,1)), color="grey90")+
    geom_text(aes(x=1.5,variable, label= round(avg_log2FC,2)),differential, size=3),
  ncol=2, rel_widths = c(1,4), align = "h", axis="tb")
fancy_dotplot

save_plts(fancy_dotplot, "dot_plot_celltype_CBA_panel", w=8,h=10)







  
            