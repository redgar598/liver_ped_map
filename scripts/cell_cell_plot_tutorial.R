source("scripts/generalized_cell_cell_interaction_UMAP.R")




load("data/example_myeloid_seurat.RData")

d10x <- NormalizeData(d10x)
d10x <- FindVariableFeatures(d10x, selection.method = "vst", nfeatures = 2000)
d10x <- ScaleData(d10x) 
d10x <- RunPCA(d10x, ndims.print = 1:10, nfeatures.print = 10)
d10x <- RunUMAP(d10x, dims = 1:30)
head(d10x)



## simple cell-cell interaction data structure
cell_cell_connections<-data.frame(Cell1=c("cDC2","Monocytes","Monocytes"), Cell2=c("cDC1","cDC1","Monocytes"))




#test examples 
plot_gene_UMAP_color_bycelltype(d10x, "annotation",cell_cell_connections)
plot_gene_UMAP_color_bycelltype(d10x, "annotation",cell_cell_connections, self_interactions = T)
plot_gene_UMAP_color_bycelltype(d10x, "annotation",cell_cell_connections,ligand_cell_type = "cDC2", label_cell_type=F)


plot_gene_UMAP_exp_colored(d10x,"annotation" ,cell_cell_connections, ligand = "CCL3", receptor="CCR1")
plot_gene_UMAP_exp_colored(d10x,"annotation" ,cell_cell_connections, ligand = "CCL3", receptor="CCR1", percentile = 0.9)
plot_gene_UMAP_exp_colored(d10x,"annotation" ,cell_cell_connections, ligand = "CCL3", receptor="CCR1",ligand_cell_type = "cDC2")



## generate cell-cell interactions from from cellphonedb output
cpdb_output<-read.table("data/example_statistical_analysis_significant_means_03_05_2024_15:06:27.txt", sep="\t", header=T)
cell_cell_connections<-cell_cell_format(cpdb_output, "CCR1","CCL3")

plot_gene_UMAP_exp_colored(d10x,"annotation" ,cell_cell_connections, ligand = "CCL3", receptor="CCR1")
plot_gene_UMAP_exp_colored(d10x,"annotation" ,cell_cell_connections,ligand_cell_type = "pDC",ligand = "CCL3", receptor="CCR1")




