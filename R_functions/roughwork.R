
##########
## Hepatocyte refining
##########
d10x.combined_hep<-subset(d10x.combined, subset = CellType_rough %in% c("Hepatocytes"))
d10x.combined_hep <- RunPCA(d10x.combined_hep, npcs = 30, verbose = FALSE)
d10x.combined_hep <- RunUMAP(d10x.combined_hep, reduction = "pca", dims = 1:30)

hep_UMAP<-DimPlot(d10x.combined_hep, reduction = "umap", pt.size=0.25, label=T)
save_plts(hep_UMAP, "hep_UMAP", w=6,h=4)



DotPlot(object = d10x.combined_hep, features = B_genes)
DotPlot(object = d10x.combined_hep, features = T_genes)
DotPlot(object = d10x.combined_hep, features = NK_genes)
DotPlot(object = d10x.combined_hep, features = LEC_genes)
DotPlot(object = d10x.combined_hep, features = Hepatocyte_genes)
DotPlot(object = d10x.combined_hep, features = NK_T_B_genes)
DotPlot(object = d10x.combined_hep, features = Cholangiocytes_genes)
DotPlot(object = d10x.combined_hep, features = HSCs_genes)
DotPlot(object = d10x.combined_hep, features = Macrophage_genes)

FeaturePlot(d10x.combined_hep, reduction = "umap", features = NK_genes, ncol = 2)
FeaturePlot(d10x.combined_hep, reduction = "umap", features = gd_genes)
FeaturePlot(d10x.combined_hep, reduction = "umap", features = T_genes, ncol = 2)
FeaturePlot(d10x.combined_hep, reduction = "umap", features = B_genes, ncol = 2)
FeaturePlot(d10x.combined_hep, reduction = "umap", features = c("IGKC","IGHG1","IGLC2", "IGHA1"), ncol = 2)
FeaturePlot(d10x.combined_hep, reduction = "umap", features = Macrophage_genes[1:6], ncol = 2)
FeaturePlot(d10x.combined_hep, reduction = "umap", features = Macrophage_genes[7:13], ncol = 2)

## key LSEC
FeaturePlot(d10x.combined_hep, reduction = "umap", features = LEC_genes[1:6], ncol = 2)
FeaturePlot(d10x.combined_hep, reduction = "umap", features = LEC_genes[7:12], ncol = 2)

## key HSC
FeaturePlot(d10x.combined_hep, reduction = "umap", features = HSCs_genes[1:6], ncol = 2)
FeaturePlot(d10x.combined_hep, reduction = "umap", features = HSCs_genes[7:12], ncol = 2)
FeaturePlot(d10x.combined_hep, reduction = "umap", features = HSCs_genes[13:20], ncol = 2)
FeaturePlot(d10x.combined_hep, reduction = "umap", features = HSCs_genes[21:27], ncol = 2)


#28 Immunoglobin_Hep
#8 Tcell_Hep
#4 LSEC_Hep
