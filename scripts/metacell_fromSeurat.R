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

#library("metacell")



source("scripts/00_pretty_plots.R")
source("scripts/00_entropy_d10x.R")

#dataset_loc <- here("/media/redgar/Seagate Portable Drive/ped_liver_map_raw")
dataset_loc <- here("../../../projects/macparland/RE/PediatricAdult/ped_liver_map_raw")

samples<-list.files(dataset_loc)
samples<-samples[-grep("data_transfer",samples)]
print(samples)

#meta<-read.table(here("data/data_transfer_updated_jan16_2023.csv"), header=T, sep=",")
meta<-read.table(here(dataset_loc,"data_transfer_updated_jan16_2023.csv"), header=T, sep=",")

d10x.list <- sapply(1:length(samples), function(y){
  caud<-meta$Sample_ID[which(meta$file == samples[y])]
  print(caud)
  print(file.path(dataset_loc,paste(samples[y],"/outs", sep=""),"filtered_feature_bc_matrix"))
  d10x <- Read10X(file.path(dataset_loc,paste(samples[y],"/outs", sep=""),"filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),caud,sep="-")
  # print(dim(d10x))
  #' Initialize the Seurat object with the raw (non-normalized data).
  d10x<-CreateSeuratObject(counts = d10x, project = "ped_adult_map", min.cells = 0, min.features = 0)
})

d10x.list

##########
### Basic QC
##########
invisible(lapply(1:length(d10x.list), function(x){
  d10x.list[[x]][["percent.mt"]] <<- PercentageFeatureSet(d10x.list[[x]], pattern = "^MT-")}))
invisible(lapply(1:length(d10x.list), function(x){
  d10x.list[[x]] <<- subset(d10x.list[[x]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 25)
}))

d10x <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "adult_ped_map")#add.cell.ids = alldata_names2,

d10x
saveRDS(d10x, file = here("data","d10x_adult_ped_raw_noSoupX.rds"))
# 
# ##########
# ## Metacell outputs
# ##########
# scdb_init(here("data/metacell"), force_reinit=T)
# scfigs_init(here("figures/metacell"))
# 
# 
# sce = as.SingleCellExperiment(d10x)
# 
# mat = scm_import_sce_to_mat(sce)
# scdb_add_mat('liver', mat)
# 
# ## Exploring and filtering the UMI matrix 
# #To get a basic understanding of the new data, we will plot the distribution of UMI  count per cell (the plot is thresholded at 800 UMIs by default):
# mcell_plot_umis_per_cell("liver")
# 
# 
# 
# #We want to clean some known issues from the matrix before starting to work with it. We generate a list of mitochondrial genes that typically mark cells as being stressed or dying, as well as immunoglobulin genes that may represent strong clonal signatures in plasma cells, rather than cellular identity.
# mat = scdb_mat("liver")
# nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
# ig_genes = c(grep("^IGJ", nms, v=T), 
#              grep("^IGH",nms,v=T),
#              grep("^IGK", nms, v=T), 
#              grep("^IGL", nms, v=T))
# 
# bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes))
# bad_genes
# 
# #We will next ask the package to ignore the above genes:
# mcell_mat_ignore_genes(new_mat_id="liver", mat_id="liver", bad_genes, reverse=F) 
# 
# 
# 
# ## Selecting feature genes
# #We move on to computing statistics on the distributions of each gene in the data, 
# #which are going to be our main tool for selecting feature genes for MetaCell analysis:
# mcell_add_gene_stat(gstat_id="liver", mat_id="liver", force=T)
# 
# #This generates a new object of type gstat under the name "liver", by analyzing the count matrix with id "liver". 
# #We can explore interesting genes and their distributions, or move directly to select a gene set for downstream analysis. For now, let's to the latter. 
# #We create a new object of type gset (gene set), to which all genes whose scaled variance (variance divided by mean) exceeds a given threshold are added:
# mcell_gset_filter_varmean(gset_id="liver_feats", gstat_id="liver", T_vm=0.08, force_new=T)
# mcell_gset_filter_cov(gset_id = "liver_feats", gstat_id="liver", T_tot=100, T_top3=2)
# 
# mcell_plot_gstats(gstat_id="liver", gset_id="liver_feats")
# 
# ## Building the balanced cell graph
# #Assuming we are happy with the selected genes,
# #we will move forward to create a similarity graph (cgraph), using a construction called balanced K-nn graph:
# 
# mcell_add_cgraph_from_mat_bknn(mat_id="liver", 
#                                gset_id = "liver_feats", 
#                                graph_id="liver_graph",
#                                K=100,
#                                dsamp=T)
# 
# 
# 
# ## Resampling and generating the co-clustering graph
# #The next step will use the cgraph to sample five hundred metacell partitions, each covering 75% of the cells and organizing them in dense subgraphs:
# 
# mcell_coclust_from_graph_resamp(
#   coc_id="liver_coc500", 
#   graph_id="liver_graph",
#   min_mc_size=20, 
#   p_resamp=0.75, n_resamp=500)
# 
# #The resampling procedure creates a new coclust object in the database named _liver_coc500_, and stores the number of times each pair of cells ended up being part of the same metacell. The co-clustering statistics are used to generate a new similarity graph, based on which accurate calling of the final set of metacells is done:
# mcell_mc_from_coclust_balanced(
#   coc_id="liver_coc500", 
#   mat_id= "liver",
#   mc_id= "liver_mc", 
#   K=30, min_mc_size=30, alpha=2)
# 
# 
# ## Removing outlier cells
# #We now have a preliminary metacell object. It is a good practice to make sure all metacells within it are homogeneous. This is done by the outlier scan procedure, which splits metacells whose underlying similarity structure supports the existence of multiple sub-clusters, and removes outlier cells that strongly deviate from their metacell's expression profile.
#             
#             mcell_plot_outlier_heatmap(mc_id="liver_mc", mat_id = "liver", T_lfc=3)
#             mcell_mc_split_filt(new_mc_id="liver_mc_f", 
#                                 mc_id="liver_mc",
#                                 mat_id="liver",
#                                 T_lfc=3, plot_mats=F)
# 
# 
# ## Selecting markers and coloring metacells
# mcell_gset_from_mc_markers(gset_id="liver_markers", mc_id="liver_mc_f")
# # marks_colors = read.table(system.file("extdata", "pbmc_mc_colorize.txt", package="metacell"), sep="\t", h=T, stringsAsFactors=F)
# # mc_colorize("liver_mc_f", marker_colors=marks_colors)
# 
# 
# mc = scdb_mc("liver_mc_f")
# # table(mc@colors)
# 
# 
# 
# # 
# # ## Creating a heatmap of genes and metacells
# # #We can use the colors to produce a labeled heat map, showing selected genes and their distributions over metacells, with the colored annotation shown at the bottom:
# # mcell_mc_plot_marks(mc_id="liver_mc_f", gset_id="liver_markers", mat_id="liver")
# # 
# # lfp = log2(mc@mc_fp)
# # tail(sort(lfp["CD8A",]))
# # 
# # 
# # ## Projecting metacells and cells in 2D
# # #Heat maps are useful but sometimes hard to interprets, and so we may want to visualize the similarity structure among metacells (or among cells within metacells). To this end we construct a 2D projection of the metacells, and use it to plot the metacells and key similarities between them (shown as connecting edges), as well as the cells. This plot will use the same metacell coloring we established before (and in case we improve the coloring based on additional analysis, the plot can be regenerated): 
# # mcell_mc2d_force_knn(mc2d_id="liver_2dproj",mc_id="liver_mc_f", graph_id="liver_graph")
# # tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
# # tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
# # mcell_mc2d_plot(mc2d_id="liver_2dproj")
# # 
# # 
# # 
# 
# 
# # 
# # 
# # ## Visualizing the MC confusion matrix
# # 
# # #While 2D projections are popular and intuitive (albeit sometimes misleading) ways to visualize scRNA-seq results, we can also summarize the similarity structure among metacells using a "confusion matrix" which encodes the pairwise similarities between all metacells. This matrix may capture hierarchical structures or other complex organizations among metacells. 
# # #We first create a hierarchical clustering of metacells, based on the number of similarity relations between their cells:
# # mc_hc = mcell_mc_hclust_confu(mc_id="liver_mc_f", 
# #                               graph_id="liver_graph")
# # 
# # #Next, we generate clusters of metacells based on this hierarchy, and visualize the confusion matrix and these clusters. The confusion matrix is shown at the bottom, and the top panel encodes the cluster hierarchy (subtrees in blue, sibling subtrees in gray):
# # mc_sup = mcell_mc_hierarchy(mc_id="liver_mc_f",
# #                             mc_hc=mc_hc, T_gap=0.04)
# # mcell_mc_plot_hierarchy(mc_id="liver_mc_f", 
# #                         graph_id="liver_graph", 
# #                         mc_order=mc_hc$order, 
# #                         sup_mc = mc_sup, 
# #                         width=2800, heigh=2000, min_nmc=2)
# # 
# # #It is also useful to do scatter plots for specific genes. In the example below, we're interested if the data contain Treg cells, so we plot CTLA4 vs CD4, and indeed there is a group of metacells enriched in both genes.
# # lfp = log2(mc@mc_fp)
# # plt = function(gene1, gene2, lfp, colors) 
# # {
# #   plot(lfp[gene1, ], lfp[gene2, ], pch=21, cex=3, bg=colors, xlab=gene1, ylab=gene2)
# #   text(lfp[gene1, ], lfp[gene2, ], colnames(lfp))
# #   
# # }
# # plt(gene1 = 'ALB', gene2 = 'PTPRC', lfp = lfp, colors = mc@colors)
# # 
# # 
# # # It is also useful to browse the lfp table itself. *mcell_mc_export_tab* exports the lfp table to a tab-delimited file. 
# # #It filters genes with maximal lfp value above the *T_fold* parameter. The function reports the metacell size (n_cells) and mean UMIs per metacell (mean_umis), which is a proxy for the cell size. The group field contains the annotation of the metacell (this will be informative after the metacell coloring functions we'll present below). The metadata field names (columns in the mat object @cell_metadata table) can be supplied to the *metadata_fields* parameter. The function would then breakdown cells in each metacell on the values of the metadata field. Very useful to characterize the metacells if our data has informative metadata features.
# # mcell_mc_export_tab(mc_id = "liver_sample", gstat_id = mat_id, mat_id = mat_id, T_fold=2, metadata_fields=c('Patient', 'CD8_gate', 'Sex', 'Stage', 'Histology'))
# # lfp_tab = read.table(scfigs_fn(mc_id, "log2_mc_fp", ext = "txt"), header=F, sep="\t", stringsAsFactors = F)
# # 
# # 
# 
# 
# ######################
# ## Compare to louvain
# ######################
# metacells<-as.data.frame(mc@mc)
# colnames(metacells)<-"metacell"
# metacells$cell<-rownames(metacells)
# 
# 
# d10x    <- SCTransform(d10x, verbose = F)
# d10x    <- RunPCA(d10x, verbose = F)
# d10x    <- RunUMAP(d10x, dims = 1:30, verbose = F)
# d10x    <- FindNeighbors(d10x, dims = 1:30, verbose = F)
# 
# sapply(c(0.1,0.5,1,1.5,2,2.5), function(x){
#   d10x    <<- FindClusters(d10x, verbose = T, resolution = x)
# })
# 
# 
# umap_mat<-as.data.frame(Embeddings(object = d10x, reduction = "umap"))#
# umap_mat$cell<-rownames(umap_mat)
# 
# meta<-d10x@meta.data
# meta$cell<-rownames(meta)
# 
# plt<-merge(meta, umap_mat, by="cell")
# plt<-merge(plt, metacells, by="cell")
# 
# saveRDS(plt, file = here("data","metacell_louvain_compare.rds"))
# 
# # 
# # plot_grid(
# #   ggplot(plt, aes(x=UMAP_1, y=UMAP_2)) +
# #     geom_point(aes(color=seurat_clusters)) + 
# #     theme_bw(),
# #   ggplot(plt, aes(x=UMAP_1, y=UMAP_2)) +
# #     geom_point(aes(color=as.factor(metacell))) + 
# #     theme_bw(),
# #   ncol=1, align="v")
# # 
# # library(ggsankey)
# # 
# # plt$metacell<-as.factor(plt$metacell)
# # df <- plt %>%
# #   make_long(seurat_clusters, metacell)
# # df
# # 
# # ggplot(df, aes(x = x, 
# #                next_x = next_x, 
# #                node = node, 
# #                next_node = next_node,
# #                fill = factor(node))) +
# #   geom_sankey() +
# #   theme_sankey(base_size = 16)
