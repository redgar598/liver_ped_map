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



#############
## Heat Map of differential genes
#############
d10x<-readRDS(file = here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_d10x_adult_ped_raw.rds"))
load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

d10x.combined_myeloid<-subset(d10x, subset = CellType_refined %in% c("KC Like"))

d10x.combined_myeloid <- NormalizeData(d10x.combined_myeloid)
d10x.combined_myeloid <- FindVariableFeatures(d10x.combined_myeloid, selection.method = "vst", nfeatures = 2000)
d10x.combined_myeloid <- ScaleData(d10x.combined_myeloid) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))
d10x.combined_myeloid <- RunPCA(d10x.combined_myeloid, npcs = 30, verbose = FALSE)
d10x.combined_myeloid <- RunUMAP(d10x.combined_myeloid, reduction = "pca", dims = 1:30)



d10x<-d10x.combined_myeloid
gene<-c("CCL3","APOE","IL10")
cellsubset<-"KC Like"
                        
  DefaultAssay(d10x) <- "RNA"
  meta_myeloid<-d10x@meta.data
  meta_myeloid$cell<-rownames(meta_myeloid)
  
  if(is.character(cellsubset)){meta_myeloid<-meta_myeloid[which(meta_myeloid$CellType_refined%in%cellsubset),]}
  
  cell_num_all<-as.data.frame(table(meta_myeloid$age_condition, meta_myeloid$CellType_refined))
  colnames(cell_num_all)<-c("age_condition","CellType_refined","CellCount")
  
  gene_exp<-FetchData(d10x, vars=gene)
  gene_exp$cell<-rownames(gene_exp)
  
  meta_myeloid<-meta_myeloid[,c("age_condition","CellType_refined","cell")]
  plt_myeloid<-merge(meta_myeloid, gene_exp, by='cell')
  
  melt_exp<-melt(plt_myeloid)
  
  
  ### scale each genen across cells then take mean for gene in each cell type
  scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}
  # 
  # plt_exp_scaled <- melt_exp %>% group_by(variable) %>%
  #   dplyr::mutate(scaled = scale_this(value))
  
  plt_exp_summary <- melt_exp %>% 
    group_by(variable, CellType_refined, age_condition) %>%
    dplyr::summarize(Mean = mean(value, na.rm=TRUE))
  
  plt_exp_scaled <- plt_exp_summary %>% group_by(variable) %>%
    dplyr::mutate(scaled = scale_this(Mean))
  plt_exp_summary<-as.data.frame(plt_exp_scaled)
  
  plt_exp_summary$variable<-factor(plt_exp_scaled$variable, levels=rev(gene))
  
  plt_exp_summary$age_condition<-as.factor(plt_exp_summary$age_condition)
  levels(plt_exp_summary$age_condition)<-c("Ped\nHealthy","Ped\nIFALD","Adult\nHealthy")
  
  plt_exp_summary$CellType_refined<-factor(plt_exp_scaled$CellType_refined, levels=cellsubset)
  

    heat_KC<- ggplot(plt_exp_summary, aes( age_condition,variable, fill=scaled))+
      geom_tile()+facet_grid(.~CellType_refined, scales = "free_y", space = "free_y")+
      th+theme_classic()+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8, 'RdBu')), space='Lab')(100), name="Scaled\nMean\nExpression")+
      ylab("")+xlab("")

    save_plts(heat_KC, "Heat_KC_only", w=7,h=6)
    
    
