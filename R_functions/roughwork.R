cell_label<-d10x@meta.data
cell_label$CellType_rough<-"Hepatocyte"

cell_label$CellType_rough[which(cell_label$seurat_clusters%in%c(2,14,26,25,22,15,19))]<-"macrophage"
cell_label$CellType_rough[which(cell_label$seurat_clusters%in%c(5,11,16))]<-"NK_T_B"
cell_label$CellType_rough[which(cell_label$seurat_clusters%in%c(3,28,23))]<-"LEC"
cell_label$CellType_rough[which(cell_label$seurat_clusters%in%c(20))]<-"Cholangiocytes"
cell_label$CellType_rough[which(cell_label$seurat_clusters%in%c(10,29))]<-"HSC"


cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

DimPlot(d10x, reduction = "umap",group.by="CellType_rough", pt.size=0.25, label=T)
