d10x.combined_hsc <- RunPCA(d10x.combined_hsc, npcs = 20, verbose = FALSE)
d10x.combined_hsc <- RunUMAP(d10x.combined_hsc, reduction = "pca", dims = 1:20)
d10x.combined_hsc <- FindNeighbors(d10x.combined_hsc, reduction = "pca", dims = 1:20)
d10x.combined_hsc <- FindClusters(d10x.combined_hsc, resolution = 0.01)

DimPlot(d10x.combined_hsc, label=T)
