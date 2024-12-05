## Load Libraries
library(here)
library(Seurat)

library(SCINA)
library(reshape2)
library(dplyr)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

library(sp)
library(rgeos)
library(scales)
library(viridis)
library(colorspace)

source("xenium_scripts/00_pretty_plots.R")
source("xenium_scripts/00_long_functions.R")



## the panel comes with cell type labels
xenium<-read.csv(file=here("/home/redgar/Documents/xenium_liver/data/Xenium_CombinedPanel.csv"))
custom_markers<-read.csv(here("data/Xenium_MacParlandGeneListUpdated.csv"))
custom_markers$Gene[grep("HAL", custom_markers$Gene)]<-"HAL"

load(file=here("data/cell_type_labels_BIDCell.RData"))



count_files<-c("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_2/cell_gene_matrices/2024_03_21_16_22_10/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_3/cell_gene_matrices/2024_03_21_16_38_48/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_4/cell_gene_matrices/2024_03_21_17_03_05/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C105/cell_gene_matrices/2024_05_17_13_32_03/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C85/cell_gene_matrices/2024_05_17_13_29_40/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C95/cell_gene_matrices/2024_07_04_12_35_29/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C101/cell_gene_matrices/2024_06_28_09_23_40/expr_mat.csv")



d10x.list <- sapply(count_files, function(file_path){
  counts<-read.csv(file_path)
  
  counts$X<-NULL
  rownames(counts) <- counts$cell_id
  counts$cell_id <- NULL
  
  # Transpose the data and convert to sparse matrix.
  mat <- as(t(as.matrix(counts)), "sparseMatrix")
  
  seu <- CreateSeuratObject(counts=mat)
  seu$sample<-strsplit(file_path,"/")[[1]][6]
  seu
})

d10x.list

xenium.obj <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "xenium_liver")
rm(d10x.list)
gc()

plt_umap_xenium<-plt_umap_xenium[match(colnames(xenium.obj), plt_umap_xenium$cell),]
identical(colnames(xenium.obj), plt_umap_xenium$cell)
rownames(plt_umap_xenium)<-plt_umap_xenium$cell

xenium.obj<-AddMetaData(xenium.obj, plt_umap_xenium)

####################################
## Fancy dot plot
####################################
unique(xenium.obj$CellType)

unique(xenium$Annotation)
xenium$Gene[grep("Kupf|kupf", xenium$Annotation)]
custom_markers[grep("Kupf|kupf", custom_markers$Cell.Type),]

xenium$Gene[grep("HA", xenium$Gene)]



## markers
hep_genes<-c("SERPINA1", "CYP1A2","CYP2A6","CYP2A7","CYP2E1")
Endothelial_genes<-c( "CD36", "LYVE1", "EGFL7","RAMP2")
act_mes_genes<-c( "COL1A1",  "PDGFRA","SPARC")
qui_mes_genes<-c( "AOX1",  "ADIPOR1")

MHCII_genes<-c( "HLA.DQA1",  "HLA.DQB1")
mono_genes<-c( "S100A8","S100A9","LYZ")
KC_genes<-c( "MARCO","VSIG4")

immune_gene<-c("PTPRC")
T_genes<-c( "CD3D","CD8A")
NK_genes<-c( "GNLY","NKG7")

chol_genes<-c( "EPCAM")

neutro_genes<-c( "FCGR3B")

plasma_genes<-c( "IGHG1")
B_genes<-c( "MS4A1")

cycle_genes<-c( "MKI67",  "TOP2A")



gene_exp<-FetchData(xenium.obj, vars=c(hep_genes,Endothelial_genes,act_mes_genes,qui_mes_genes,chol_genes,
                                       immune_gene,
                                       T_genes,NK_genes,B_genes,plasma_genes,neutro_genes,
                                       mono_genes,MHCII_genes, KC_genes,cycle_genes))
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)

rm(xenium.obj)
gc()

plt<-merge(gene_exp, plt_umap_xenium, by.x="cell", by.y="cell")
gc()

## summarize
scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}

plt_summary<-plt %>% group_by(CellType, variable) %>% 
  summarise(mn=mean(value), count=length(which(value>0)), percent_exp=(length(which(value>0))/length(value))*100)
plt_summary <- plt_summary %>% group_by(variable) %>%
  dplyr::mutate(scaled = scale_this(mn))
plt_summary<-as.data.frame(plt_summary)

# remove dots where 0 cell expressing marker
plt_summary<-plt_summary[(which(plt_summary$count>0)),]

plt_summary$variable<-factor(plt_summary$variable, levels=rev(c(hep_genes,Endothelial_genes,act_mes_genes,qui_mes_genes,chol_genes,
                                                                immune_gene,
                                                                T_genes,NK_genes,B_genes,plasma_genes,neutro_genes,
                                                                mono_genes,MHCII_genes, KC_genes,cycle_genes)))

plt_summary$CellType<-factor(plt_summary$CellType, levels=c("Hepatocyte (Pericentral)", "Hepatocyte (Periportal)", "Hepatocyte (Cycling)",
                                                            "LSEC","VEC",
                                                            "HSC (Activated)","HSC (Periportal)","HSC (Quiescent)",
                                                            "Cholangiocytes" ,  "Cholangiocyte (Biliary)" ,
                                                            "CD3+ T-cells", "gd T-cells", "Mature B-cells","Plasma Cells",  "Neutrophil", 
                                                            "Mono-Mac","Macrophage MHCII High","KC Like", "Cycling Myeloid","Erythrocytes"    ))

                                  
            
                             
                                               
          

gene_list_len<-sapply(list(hep_genes,Endothelial_genes,act_mes_genes,qui_mes_genes,chol_genes,
                           immune_gene,
                           T_genes,NK_genes,B_genes,plasma_genes,neutro_genes,
                           mono_genes,MHCII_genes, KC_genes,cycle_genes), function(x) length(x))    

fancy_dotplot<-plot_grid(
  ggplot(plt_summary, aes(CellType, variable, color=scaled, size=percent_exp))+geom_point()+
    theme_classic()+
    scale_color_continuous_sequential(palette = "Oslo", rev=F, name="Scaled\nMean\nExpression")+
    scale_size(name="Percent\nCells\nExpressing")+
    theme(axis.text.x = element_blank(),axis.title = element_blank(),axis.ticks.x = element_blank())+
    geom_hline(yintercept = cumsum(rev(gene_list_len))+0.5, color="grey70")+
    geom_vline(xintercept = c(3,6,9,11,15)+0.5, color="grey70"),
  ggplot(plt_summary, aes(CellType, y=1, fill=CellType))+geom_tile(color="black")+
    theme_classic()+fillscale_cellType+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",axis.line  = element_blank(),
          plot.margin = margin(t = 0,  # Top margin
                               r = 50,  # Right margin
                               b = 40,  # Bottom margin
                               l = 10)),
  ncol=1, rel_heights = c(6,1.5), align = "v", axis="lr")
fancy_dotplot

save_plts(fancy_dotplot, "xenium_dot_plot_celltype", w=8,h=12)

