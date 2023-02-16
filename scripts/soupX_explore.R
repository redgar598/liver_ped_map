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

#dataset_loc <- here("/media/redgar/Seagate Portable Drive/ped_liver_map_raw")
dataset_loc <- here("../../../projects/macparland/RE/PediatricAdult/ped_liver_map_raw")

samples<-list.files(dataset_loc)
samples<-samples[-grep("data_transfer",samples)]
print(samples)

#meta<-read.table(here("data/data_transfer_updated_jan16_2023.csv"), header=T, sep=",")
meta<-read.table(here(dataset_loc,"data_transfer_updated_jan16_2023.csv"), header=T, sep=",")

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



################
## Checksaved rel change in all samples
################
load(here("data","adult_ped_integrated_refinedlabels_withDropletQC.rds"))
head(d10x.combined)

umap_mat<-as.data.frame(Embeddings(object = d10x.combined, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x.combined@meta.data
meta$cell<-rownames(meta)

plt<-merge(meta, umap_mat, by="cell")


plt<-plt[order(plt$cell_status),]
soup_change<-plot_grid(ggplot(plt, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = CellType_refined), size=0.5)+theme_bw()+colscale_cellType+ guides(colour = guide_legend(override.aes = list(size=2))),
                       ggplot(plt, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = relALBChange), size=0.5)+theme_bw()+scale_colour_gradientn(colours=c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704'),limits=c(NA,NA), name="Relative\nALB\nChange"),
                       ggplot(plt, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = nuclear_fraction), size=0.5)+theme_bw(),
                       ggplot(plt, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = cell_status), size=0.5)+theme_bw()+ guides(colour = guide_legend(override.aes = list(size=2)))+scale_color_manual(values=c("grey","red","cornflowerblue"),name="Cell Status"),
                       ggplot(plt, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = percent.mt), size=0.5)+theme_bw(),
                       ggplot(plt, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = Phase), size=0.5)+theme_bw()+scale_color_manual(values=c("#e6ab02","#386cb0","#1b9e77"))+ guides(colour = guide_legend(override.aes = list(size=2))),
                       ncol=2, align="v")
soup_change
save_plts(soup_change, "all_sample_soup_change", w=12,h=12)



## grab legened from plot
get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

align_plots_seperate<-function(plt_file_name,UMAP_plt){
  leg_umap = get_leg(UMAP_plt)
  plt_UMAP<- grid.arrange(UMAP_plt+theme(legend.position = "none"),leg_umap, ncol=2, widths=c(0.8,0.2))
  save_plts(plt_UMAP, plt_file_name, w=7,h=5)
}


align_plots_seperate("cell_type_UMAP", ggplot(plt, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = CellType_refined), size=0.5)+theme_bw()+colscale_cellType+ guides(colour = guide_legend(override.aes = list(size=2))))
align_plots_seperate("soupX_UMAP", ggplot(plt, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = relALBChange), size=0.5)+theme_bw()+scale_colour_gradientn(colours=c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704'),limits=c(NA,NA), name="Relative\nALB\nChange"))
align_plots_seperate("nuclearfraction_UMAP",  ggplot(plt, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = nuclear_fraction), size=0.5)+theme_bw())
align_plots_seperate("cell_status_UMAP", ggplot(plt, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = cell_status), size=0.5)+theme_bw()+ guides(colour = guide_legend(override.aes = list(size=2)))+scale_color_manual(values=c("grey","red","cornflowerblue"),name="Cell Status"))
align_plots_seperate("percentMT_UMAP",  ggplot(plt, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = percent.mt), size=0.5)+theme_bw())
align_plots_seperate("phase_UMAP",ggplot(plt, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = Phase), size=0.5)+theme_bw()+scale_color_manual(values=c("#e6ab02","#386cb0","#1b9e77"))+ guides(colour = guide_legend(override.aes = list(size=2))) )



soup_dropletQC<-ggplot(plt, aes(cell_status, relALBChange))+geom_violin(fill="lightgrey",color="lightgrey")+geom_boxplot(aes(fill=cell_status), width=0.05)+
  th+theme_bw()+scale_fill_manual(values=c("grey","red","cornflowerblue"),name="Cell Status", guide="none")+
  xlab("Cell Status")+ylab("SoupX\n(relative ALB change)")
save_plts(soup_dropletQC, "soupXdropletQCassociation", w=4,h=3)





################
## Checksaved rel change in all samples
################
d10x<-readRDS(file = here("data","d10x_adult_ped_raw.rds"))
alb_exp<-FetchData(object = d10x, vars = c("ALB"))
alb_exp$cell<-rownames(alb_exp)
rm(d10x)

load(here("data","adult_ped_integrated_refinedlabels_withDropletQC.rds"))
head(d10x.combined)

umap_mat<-as.data.frame(Embeddings(object = d10x.combined, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x.combined@meta.data
meta$cell<-rownames(meta)

plt<-merge(meta, umap_mat, by="cell")
plt<-merge(plt, alb_exp, by="cell")


plt<-plt[order(plt$ALB),]
my_breaks = c(0, 10, 100, 1000, 8000)
ALB_count_after_soup<-ggplot(plt, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = ALB), size=0.5)+
  facet_wrap(~AgeGroup)+
  theme_bw()+scale_color_gradient2(name = "ALB\ncount", trans = "log",
                                  breaks = my_breaks, labels = my_breaks,
                                  low = "white",
                                  high = "blue",
                                  midpoint = 0,
                                  na.value = "grey80")

save_plts(ALB_count_after_soup, "ALB_soup_change", w=10,h=4.5)


#############################################
### Soup for all genes all samples
#############################################
soup_est.list <- lapply(1:length(samples), function(y){
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
  
  ######
  ## Load data and estimate soup profile
  ######
  # Estimate rho
  sc = autoEstCont(sc, forceAccept=TRUE)
  #Genes with highest expression in background. These are often enriched for ribosomal proteins.
  print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20))
  
  soup_est<-sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ]
  soup_est$gene<-rownames(soup_est)
  soup_est$counts<-NULL
  colnames(soup_est)[1]<-paste(caud, colnames(soup_est)[1], sep="_")
  soup_est})


soup_est_allsamp<-Reduce(function(x, y) merge(x, y, by="gene"), soup_est.list)
head(soup_est_allsamp)

save(soup_est_allsamp, file=paste(here("data/"),"soup_est_allsamp.rds", sep=""))



load(here("data","soup_est_allsamp.rds"))


## Rank sum of each gene to get top soup
soup_rank<-lapply(2:12, function(x){
  soup_est_allsamp$gene[order(soup_est_allsamp[,x], decreasing = T)]
  data.frame(gene=soup_est_allsamp$gene[order(soup_est_allsamp[,x], decreasing = T)], rank=1:nrow(soup_est_allsamp))})
soup_rank<-Reduce(function(x, y) merge(x, y, by="gene"), soup_rank)
soup_rank_sum<-data.frame(gene=soup_rank$gene, rank_sum=rowSums(soup_rank[,2:12]))
soup_rank_sum<-soup_rank_sum[order(soup_rank_sum$rank_sum),]
head(soup_rank_sum, n=20)



#########
## pathway on the mean est
#########
soup_est_allsamp$mean_est<-rowMeans(soup_est_allsamp[,2:12])

source("scripts/00_GSEA_function.R")
GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")


  gene_list = soup_est_allsamp$mean_est
  names(gene_list) = soup_est_allsamp$gene
  gene_list = sort(gene_list, decreasing = TRUE)
  gene_list = gene_list[!duplicated(names(gene_list))]
  
  res = GSEA(gene_list, GO_file, pval = 0.05)
  ## Nothing significant
  
  
#   
#   plt_path<-res$Results
#   plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
#   plt_path$Enrichment_Cell<-"Up-regulated in \nAdult"
#   plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \nPediatric"
#   
#   plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))
#   
#   plt_path$direction_label<-as.factor(plt_path$Enrichment)
#   levels(plt_path$direction_label)<-c(0.1,-0.1)
#   plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))
#   
#   # top and bottom 15
#   plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])
#   
#   ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
#     theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
#     geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
#     geom_vline(xintercept = 0, color="grey40")+scale_fill_manual(values=c("#fd8d3c","#6baed6"))+ 
#     guides(fill = guide_legend(override.aes = list(size=5)))
# 
# RR_GSEA<-pathway_plt(diff_exp_all[which(diff_exp_all$cell.1=="RR Myeloid_Adult"),])
# save_plts(RR_GSEA, "GSEA_adult_ped_recently_recruited", w=20,h=10)
# 
