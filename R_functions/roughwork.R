d10x.combined@meta.data$CellType_rough[which(d10x.combined@meta.data$CellType_rough=="macrophage")]<-"Myeloid"
d10x.combined@meta.data$CellType_rough[which(d10x.combined@meta.data$seurat_clusters==5)]<-"CD3_Tcell"
d10x.combined@meta.data$CellType_rough[which(d10x.combined@meta.data$seurat_clusters==11)]<-"nkTcell"
d10x.combined@meta.data$CellType_rough[which(d10x.combined@meta.data$seurat_clusters==16)]<-"gdTcell"
