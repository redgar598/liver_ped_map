B_genes<-c("POU2F2","FCER2","MS4A1","LTB","CD37","CD79B","IGLC2","IGHG1","IGKC")
T_genes<-c("CD3D","IL7R","CD8A","IL32")
NK_genes<-c("NKG7","CD7")
gd_genes<-c("GNLY")

d10x.combined_NK_T_B<-subset(d10x.combined, subset = CellType_rough %in% c("NK_T_B"))
d10x.combined_NK_T_B <- RunPCA(d10x.combined_NK_T_B, npcs = 30, verbose = FALSE)
d10x.combined_NK_T_B <- RunUMAP(d10x.combined_NK_T_B, reduction = "pca", dims = 1:30)

DimPlot(d10x.combined_NK_T_B, reduction = "umap", pt.size=0.25)

FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = gd_genes, ncol = 2)
FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = T_genes, ncol = 2)
FeaturePlot(d10x.combined_NK_T_B, reduction = "umap", features = NK_genes, ncol = 2)


## same method on sub cluster
genes<-unique(c(T_genes, NK_genes,gd_genes))

d10x.exp<-as.data.frame(d10x.combined_NK_T_B[["RNA"]]@data)
d10x.exp.GOI<-d10x.exp[genes,]
d10x.exp.GOI$gene<-rownames(d10x.exp.GOI)
d10x.exp.GOI<-melt(d10x.exp.GOI)#

meta<-d10x.combined_NK_T_B@meta.data
meta$cell<-rownames(meta)

plt<-merge(d10x.exp.GOI, meta,by.x="variable", by.y="cell")

plt$variable<-as.character(plt$variable)
plt$seurat_clusters<-as.character(plt$seurat_clusters)


## possible cells types
cluster_marker_mean<-function(gene_list, type){
  plt_epi<-plt[which(plt$gene%in%gene_list),]
  mean_type<-as.data.frame(tapply(plt_epi$value, plt_epi$seurat_clusters, mean))
  colnames(mean_type)<-type
  mean_type
}

cell_rough<-cbind(cluster_marker_mean(T_genes, "CD3_Tcell"),
                  cluster_marker_mean(NK_genes, "nkTcell"),
                  cluster_marker_mean(gd_genes, "gdTcell"))


cell_rough$CellType_rough_sub<-sapply(1:nrow(cell_rough), function(x) {
  compart<-colnames(cell_rough)[which(cell_rough[x,] == max(cell_rough[x,]))]
  if(length(compart)==1){compart}else{"Unclear"}
})

cell_rough$seurat_clusters<-rownames(cell_rough)

meta<-d10x.combined_NK_T_B@meta.data
meta$cell<-rownames(meta)
plt_summary<-merge(meta, cell_rough[,c("seurat_clusters","CellType_rough_sub")], by="seurat_clusters")
plt_summary$CellType_rough<-plt_summary$CellType_rough_sub
plt_summary$CellType_rough_sub<-NULL
plt_summary<-plt_summary[match(rownames(d10x.combined_NK_T_B@meta.data),plt_summary$cell),]
identical(plt_summary$cell, rownames(d10x.combined_NK_T_B@meta.data))
rownames(plt_summary)<-plt_summary$cell
d10x.combined_NK_T_B<- AddMetaData(d10x.combined_NK_T_B, plt_summary)

DimPlot(d10x.combined_NK_T_B, group.by="CellType_rough_sub", reduction = "umap", pt.size=0.25)



meta<-d10x.combined@meta.data
meta$cell<-rownames(meta)
meta_not_ntb<-meta[which(meta$CellType_rough!="NK_T_B"),]

meta_ntb<-d10x.combined_NK_T_B@meta.data
meta_ntb$CellType_rough<-meta_ntb$CellType_rough_sub

cell_label<-rbind(meta_not_ntb, meta_ntb)


plt_summary<-merge(meta, cell_rough[,c("seurat_clusters","CellType_rough")], by="seurat_clusters")
plt_summary<-plt_summary[match(rownames(d10x.combined@meta.data),plt_summary$cell),]
identical(plt_summary$cell, rownames(d10x.combined@meta.data))

rownames(plt_summary)<-plt_summary$cell

d10x.combined<- AddMetaData(d10x.combined, plt_summary)
