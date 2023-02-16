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

library("metacell")



source("scripts/00_pretty_plots.R")
source("scripts/00_entropy_d10x.R")


scdb_init("/home/redgar/Documents/liver_ped_map/data/testmetacell/", force_reinit=T)

## load one sample
mcell_import_scmat_10x("test", base_dir="/home/redgar/Documents/liver_ped_map/data/testmetacell/MacParland_Diana__C93_Frozen_Liver_220919_3pr_V3_1/filtered_feature_bc_matrix")
mat = scdb_mat("test")
print(dim(mat@mat))

scfigs_init("/home/redgar/Documents/liver_ped_map/figures/testmetacell/")

## Exploring and filtering the UMI matrix 
#To get a basic understanding of the new data, we will plot the distribution of UMI  count per cell (the plot is thresholded at 800 UMIs by default):
mcell_plot_umis_per_cell("test")



#We want to clean some known issues from the matrix before starting to work with it. We generate a list of mitochondrial genes that typically mark cells as being stressed or dying, as well as immunoglobulin genes that may represent strong clonal signatures in plasma cells, rather than cellular identity.
mat = scdb_mat("test")
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
ig_genes = c(grep("^IGJ", nms, v=T), 
             grep("^IGH",nms,v=T),
             grep("^IGK", nms, v=T), 
             grep("^IGL", nms, v=T))

bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes))
bad_genes

#We will next ask the package to ignore the above genes:
mcell_mat_ignore_genes(new_mat_id="test", mat_id="test", bad_genes, reverse=F) 


#In the current example we will also eliminate cells with less than 800 UMIs (threshold can be set based on examination of the UMI count distribution):
mcell_mat_ignore_small_cells("test", "test", 800)



## Selecting feature genes
#We move on to computing statistics on the distributions of each gene in the data, which are going to be our main tool for selecting feature genes for MetaCell analysis:
mcell_add_gene_stat(gstat_id="test", mat_id="test", force=T)

#This generates a new object of type gstat under the name "test", by analyzing the count matrix with id "test". We can explore interesting genes and their distributions, or move directly to select a gene set for downstream analysis. For now, let's to the latter. 
#We create a new object of type gset (gene set), to which all genes whose scaled variance (variance divided by mean) exceeds a given threshold are added:
mcell_gset_filter_varmean(gset_id="test_feats", gstat_id="test", T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id = "test_feats", gstat_id="test", T_tot=100, T_top3=2)

mcell_plot_gstats(gstat_id="test", gset_id="test_feats")

## Building the balanced cell graph
#Assuming we are happy with the selected genes,
#we will move forward to create a similarity graph (cgraph), using a construction called balanced K-nn graph:

mcell_add_cgraph_from_mat_bknn(mat_id="test", 
                               gset_id = "test_feats", 
                               graph_id="test_graph",
                               K=100,
                               dsamp=T)



## Resampling and generating the co-clustering graph
#The next step will use the cgraph to sample five hundred metacell partitions, each covering 75% of the cells and organizing them in dense subgraphs:

mcell_coclust_from_graph_resamp(
  coc_id="test_coc500", 
  graph_id="test_graph",
  min_mc_size=20, 
  p_resamp=0.75, n_resamp=500)

#The resampling procedure creates a new coclust object in the database named _test_coc500_, and stores the number of times each pair of cells ended up being part of the same metacell. The co-clustering statistics are used to generate a new similarity graph, based on which accurate calling of the final set of metacells is done:
mcell_mc_from_coclust_balanced(
  coc_id="test_coc500", 
  mat_id= "test",
  mc_id= "test_mc", 
  K=30, min_mc_size=30, alpha=2)


## Removing outlier cells
#We now have a preliminary metacell object. It is a good practice to make sure all metacells within it are homogeneous. This is done by the outlier scan procedure, which splits metacells whose underlying similarity structure supports the existence of multiple sub-clusters, and removes outlier cells that strongly deviate from their metacell's expression profile.

                    mcell_plot_outlier_heatmap(mc_id="test_mc", mat_id = "test", T_lfc=3)
                    mcell_mc_split_filt(new_mc_id="test_mc_f", 
                                mc_id="test_mc",
                                mat_id="test",
                                T_lfc=3, plot_mats=F)


## Selecting markers and coloring metacells
mcell_gset_from_mc_markers(gset_id="test_markers", mc_id="test_mc_f")
marks_colors = read.table(system.file("extdata", "pbmc_mc_colorize.txt", package="metacell"), sep="\t", h=T, stringsAsFactors=F)
mc_colorize("test_mc_f", marker_colors=marks_colors)
mc = scdb_mc("test_mc_f")
table(mc@colors)



## Creating a heatmap of genes and metacells
#We can use the colors to produce a labeled heat map, showing selected genes and their distributions over metacells, with the colored annotation shown at the bottom:
mcell_mc_plot_marks(mc_id="test_mc_f", gset_id="test_markers", mat_id="test")

lfp = log2(mc@mc_fp)
tail(sort(lfp["CD8A",]))


## Projecting metacells and cells in 2D
#Heat maps are useful but sometimes hard to interprets, and so we may want to visualize the similarity structure among metacells (or among cells within metacells). To this end we construct a 2D projection of the metacells, and use it to plot the metacells and key similarities between them (shown as connecting edges), as well as the cells. This plot will use the same metacell coloring we established before (and in case we improve the coloring based on additional analysis, the plot can be regenerated): 
mcell_mc2d_force_knn(mc2d_id="test_2dproj",mc_id="test_mc_f", graph_id="test_graph")
tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
mcell_mc2d_plot(mc2d_id="test_2dproj")







## Visualizing the MC confusion matrix

#While 2D projections are popular and intuitive (albeit sometimes misleading) ways to visualize scRNA-seq results, we can also summarize the similarity structure among metacells using a "confusion matrix" which encodes the pairwise similarities between all metacells. This matrix may capture hierarchical structures or other complex organizations among metacells. 
#We first create a hierarchical clustering of metacells, based on the number of similarity relations between their cells:
mc_hc = mcell_mc_hclust_confu(mc_id="test_mc_f", 
                              graph_id="test_graph")

#Next, we generate clusters of metacells based on this hierarchy, and visualize the confusion matrix and these clusters. The confusion matrix is shown at the bottom, and the top panel encodes the cluster hierarchy (subtrees in blue, sibling subtrees in gray):
mc_sup = mcell_mc_hierarchy(mc_id="test_mc_f",
                            mc_hc=mc_hc, T_gap=0.04)
mcell_mc_plot_hierarchy(mc_id="test_mc_f", 
                        graph_id="test_graph", 
                        mc_order=mc_hc$order, 
                        sup_mc = mc_sup, 
                        width=2800, heigh=2000, min_nmc=2)

#It is also useful to do scatter plots for specific genes. In the example below, we're interested if the data contain Treg cells, so we plot CTLA4 vs CD4, and indeed there is a group of metacells enriched in both genes.
lfp = log2(mc@mc_fp)
plt = function(gene1, gene2, lfp, colors) 
{
  plot(lfp[gene1, ], lfp[gene2, ], pch=21, cex=3, bg=colors, xlab=gene1, ylab=gene2)
  text(lfp[gene1, ], lfp[gene2, ], colnames(lfp))
  
}
plt(gene1 = 'ALB', gene2 = 'PTPRC', lfp = lfp, colors = mc@colors)


# It is also useful to browse the lfp table itself. *mcell_mc_export_tab* exports the lfp table to a tab-delimited file. It filters genes with maximal lfp value above the *T_fold* parameter. The function reports the metacell size (n_cells) and mean UMIs per metacell (mean_umis), which is a proxy for the cell size. The group field contains the annotation of the metacell (this will be informative after the metacell coloring functions we'll present below). The metadata field names (columns in the mat object @cell_metadata table) can be supplied to the *metadata_fields* parameter. The function would then breakdown cells in each metacell on the values of the metadata field. Very useful to characterize the metacells if our data has informative metadata features.
mcell_mc_export_tab(mc_id = "test_sample", gstat_id = mat_id, mat_id = mat_id, T_fold=2, metadata_fields=c('Patient', 'CD8_gate', 'Sex', 'Stage', 'Histology'))
lfp_tab = read.table(scfigs_fn(mc_id, "log2_mc_fp", ext = "txt"), header=F, sep="\t", stringsAsFactors = F)




######################
## Compare to louvain
######################
metacells<-as.data.frame(mc@mc)
colnames(metacells)<-"metacell"
metacells$cell<-rownames(metacells)


d10x <- Read10X(data.dir="/home/redgar/Documents/liver_ped_map/data/testmetacell/MacParland_Diana__C93_Frozen_Liver_220919_3pr_V3_1/filtered_feature_bc_matrix")
d10x<-CreateSeuratObject(counts = d10x, project = "ped_adult_map", min.cells = 0, min.features = 0)

d10x    <- SCTransform(d10x, verbose = F)
d10x    <- RunPCA(d10x, verbose = F)
d10x    <- RunUMAP(d10x, dims = 1:30, verbose = F)
d10x    <- FindNeighbors(d10x, dims = 1:30, verbose = F)
d10x    <- FindClusters(d10x, verbose = T, resolution = 0.8)

umap_mat<-as.data.frame(Embeddings(object = d10x, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x@meta.data
meta$cell<-rownames(meta)

plt<-merge(meta, umap_mat, by="cell")
plt<-merge(plt, metacells, by="cell")

plot_grid(
  ggplot(plt, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(aes(color=seurat_clusters)) + 
  theme_bw(),
  ggplot(plt, aes(x=UMAP_1, y=UMAP_2)) +
    geom_point(aes(color=as.factor(metacell))) + 
    theme_bw(),
  ncol=1, align="v")

remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)

plt$metacell<-as.factor(plt$metacell)
df <- plt %>%
  make_long(seurat_clusters, metacell)
df

ggplot(df, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node))) +
  geom_sankey() +
  theme_sankey(base_size = 16)
