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

dataset_loc <- here("/media/redgar/Seagate Portable Drive/ped_liver_map_raw")

samples<-list.files(dataset_loc)
print(samples)

meta<-read.table(here("data","data_transfer_updated_jan16_2023.csv"), header=T, sep=",")

y=3
  caud<-meta$Sample_ID[which(meta$file == samples[y])]
  print(caud)
  print(file.path(dataset_loc,paste(samples[y],"/outs", sep=""),"filtered_feature_bc_matrix"))
  d10x <- Read10X(file.path(dataset_loc,paste(samples[y],"/outs", sep=""),"filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),caud,sep="-")
  # print(dim(d10x))
  #' Initialize the Seurat object with the raw (non-normalized data).
  d10x<-CreateSeuratObject(counts = d10x, project = "ped_adult_map", min.cells = 0, min.features = 0)
  
  ## SoupX needs clusters so quickly make clusters for each sample
  d10x    <- SCTransform(d10x, verbose = F)
  d10x    <- RunPCA(d10x, verbose = F)
  d10x    <- RunUMAP(d10x, dims = 1:30, verbose = F)
  d10x    <- FindNeighbors(d10x, dims = 1:30, verbose = F)
  d10x    <- FindClusters(d10x, verbose = T)
  meta_clusters    <- d10x@meta.data
  
  sc = load10X(file.path(dataset_loc,paste(samples[y],"/outs", sep="")))
  sc = setClusters(sc, setNames(meta_clusters$seurat_clusters, rownames(meta_clusters)))

  umap_mat<-as.data.frame(Embeddings(object = d10x, reduction = "umap"))#
  umap_mat$cell<-rownames(umap_mat)
  
  sc = setDR(sc, umap_mat[, c("UMAP_1", "UMAP_2")])
  
  
  ######
  ## Load data and estimate soup profile
  ######
  # Estimate rho
  sc = autoEstCont(sc, forceAccept=TRUE)
  #Genes with highest expression in background. These are often enriched for ribosomal proteins.
  print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20))
  # Clean the data
  out = adjustCounts(sc)
  
  load(here("data","adult_ped_cellRefined_withDropletQC.rds"))
  cell_label$cell<-rownames(cell_label)
  cell_label$cell<-gsub("1_3","C93_caud3pr",cell_label$cell)

  dd = d10x@meta.data
  dd<-cbind(dd, umap_mat)
  dd<-merge(dd, cell_label, by="cell", all.x=T)
  dd$ALB = sc$toc["ALB", ]
  dd$PTPRC = sc$toc["PTPRC", ]
  dd$CALCRL = sc$toc["CALCRL", ]
  
soup_one<-plot_grid(ggplot(dd, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = CellType_refined), size=0.25)+theme_bw()+colscale_cellType+ guides(colour = guide_legend(override.aes = list(size=2))),
          ggplot(dd, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = ALB > 0), size=0.2)+theme_bw()+scale_color_manual(values=c("grey","darkred"))+ guides(colour = guide_legend(override.aes = list(size=2))),
          ggplot(dd, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = PTPRC > 0), size=0.2)+theme_bw()+scale_color_manual(values=c("grey","darkred"))+ guides(colour = guide_legend(override.aes = list(size=2))),
          ggplot(dd, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = CALCRL > 0), size=0.2)+theme_bw()+scale_color_manual(values=c("grey","darkred"))+ guides(colour = guide_legend(override.aes = list(size=2))), 
          ncol=2, align="v")
soup_one
save_plts(soup_one, "one_sample_soup_before", w=11,h=7)

  

  plotMarkerMap(sc, "ALB")+theme_bw()
  plotChangeMap(sc, out, "ALB")+theme_bw()
  

DR = sc$metaData[,sc$DR]
df = DR
old = colSums(sc$toc["ALB",rownames(df),drop=FALSE])
new = colSums(out["ALB",rownames(df),drop=FALSE])
relChange = (old-new)/old
df$old = old
df$new = new
df$relALBChange=relChange
df$cell<-rownames(df)
#Stick the NAs at the bottom
df = df[order(!is.na(df$relALBChange)),]


soup_change<-plot_grid(ggplot(dd, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = CellType_refined), size=0.5)+theme_bw()+colscale_cellType+ guides(colour = guide_legend(override.aes = list(size=2))),
                    ggplot(df, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = relALBChange), size=0.5)+theme_bw()+scale_colour_gradientn(colours=c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704'),limits=c(NA,NA), name="Relative\nALB\nChange")+ guides(colour = guide_legend(override.aes = list(size=2))),
                    ncol=2, align="v")
soup_change
save_plts(soup_change, "one_sample_soup_change", w=12,h=4)



d10x = CreateSeuratObject(out)

#add meta data to each seurat object
meta_cell<-data.frame(cell=colnames(d10x), individual=caud)
meta_cell_add<-merge(meta_cell, meta, by.x="individual", by.y="Sample_ID")
meta_cell_add<-merge(meta_cell_add, df, by="cell")
meta_cell_add<-meta_cell_add[match(colnames(d10x), meta_cell_add$cell),]
print(identical(meta_cell_add$cell, colnames(d10x)))
rownames(meta_cell_add)<-meta_cell_add$cell
d10x<- AddMetaData(d10x, meta_cell_add)
