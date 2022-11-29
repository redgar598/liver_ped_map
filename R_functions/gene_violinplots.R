
#'### Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(here)
library(ggplot2)
library(reshape2)
library(cowplot)


##########
## Objects used
##########
# d10x is a seurat object
# CellType_rough is my cell type column
cell_types<-unique(d10x$CellType_rough)

# diff_exp_all is the output of seurat FindMarkers and used to draw a box around significant genes

myeloid_immune_supressive<-c("CTSB","CD163","MS4A7","FOLR2","GPNMB","VSIG4","HMOX1","MSR1")
inflammatory_macs<-c("CD74","HLA-DRA","TYROBP","C1QC","HLA-DPA1","HLA-DPB1","LYZ","S100A6")
exhausted_tcells<-c("TOX","PDCD1","LAG3","TNFRSF9","CXCL13","ENTPD1","HAVCR2","CD38")


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


## plot each signature
plot_key_genes(myeloid_immune_supressive, "myeloid_immune_supressive")
plot_key_genes(inflammatory_macs, "inflammatory_macs")
plot_key_genes(exhausted_tcells, "exhausted_tcells")
