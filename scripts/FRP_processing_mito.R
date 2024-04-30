####################
## Random side analysis
####################

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
library(SCINA)



save_plts<-function(plt, name, w,h){
  ggsave(plt, file=paste(here("../FRP_mito/figures/"),name,".pdf", sep=""), w=w, h=h)
  ggsave(plt, file=paste(here("../FRP_mito/figures/jpeg/"),name,".jpeg", sep=""), w=w, h=h)
  ggsave(plt, file=paste(here("../FRP_mito/figures/png/"),name,".png", sep=""), w=w, h=h)}



myColors_celltype <- c("#660cc7","#8642d1","#5612a3","#4b911d","#7a4202",
                       "#2c6e02","#6bbce8","#0b3349","#b8d8dc",
                       "#3469ad","#d17906","#b01e87",
                       "#60ba5d","#207537","#a0c487","#d9a5a5","#87a4c9",
                       "#e8c392","#dea4ce","#79639a","#207537","#fa61ad",
                       "#b80783","#994676","#431039","#cb181d","maroon1",
                       "#b01629","grey","#ce1256","#a6d96a","#750c32",
                       "#d9667f","#1b4003","#e0a8ce","#8a68b0","#3d1b63",
                       "#c9a8ed","#c48db4","#a3588d","#6ca647","#3a7d31",
                       "#60ba5d","#3d1b63","#2dc918","#F4355B","#e4d5f2")
color_possibilities_celltype<-c("B-cells","Mature B-cells","Plasma cells","CD3+ T-cells","Cholangiocytes",
                                "gd T-cells","Hepatocytes","Hepatocytes (Central)","Hepatocytes (Portal)",
                                "HSC","LSEC","Myeloid cells",
                                "NK-like cells", "NK and T cells","NKT cells\n(Hepatocyte Like)","Cholangiocytes\n(Hepatocyte Like)","HSC\n(Hepatocyte Like)",
                                "LSEC\n(Hepatocyte Like)","Myeloid cells\n(Hepatocyte Like)","B-cells\n(Hepatocyte Like)","NKT cells","Macrophage",
                                "KC Like","Neutrophil","Neutrophil\n(DEFA+)","Erythrocytes","Mast cell",
                                "Myeloid Erythrocytes\n(phagocytosis)","Doublet","Macrophage\n(MHCII high)","Cycling T-cells","CDC1",
                                "Platelets","CLNK T-cells","Cycling Myeloid","Mature B-cells (104)","pDC",
                                "pre B-cell","KC Like\n(Hepatocyte Like)","KC Like (C97)","Naive CD4 T-cells","Memory CD4 T-cells",
                                "NK cells","DC","CD8 T-cells","CD14+ Mono","Cycling Plasma")
names(myColors_celltype) <- color_possibilities_celltype
fillscale_cellType <- scale_fill_manual(name="Cell Type",
                                        values = myColors_celltype, drop = T, limits=force)
colscale_cellType <- scale_color_manual(name="Cell Type",
                                        values = myColors_celltype, drop = T, limits=force)



cellranger_out_loc <- here("/media/redgar/Seagate Portable Drive/FRP_mito")
meta<-read.table(here("../FRP_mito/FRP_metadata.csv"), header=T, sep=",")

d10x.list <- sapply(1:nrow(meta), function(y){
  sampl<-meta$Sample[y]
  print(sampl)
  print(file.path(cellranger_out_loc, meta$Location[y]))
  d10x <- Read10X(file.path(cellranger_out_loc, meta$Location[y]))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),sampl,sep="-")
  # print(dim(d10x))
  d10x<-CreateSeuratObject(counts = d10x, project = "FRP", min.cells = 0, min.features = 0)

  ## soupX
  if(dir.exists(file.path(dataset_loc,paste(samples[y],"/outs/raw_feature_bc_matrix", sep="")))) {
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
    print(unique(sc$metaData$rho))
    # Clean the data
    out = adjustCounts(sc)
    
    ## Save a metric of soupness (ALB change after soupX)
    DR = sc$metaData[,sc$DR]
    df = DR
    old = colSums(sc$toc["ALB",rownames(df),drop=FALSE])
    new = colSums(out["ALB",rownames(df),drop=FALSE])
    relChange = (old-new)/old
    df$old = old
    df$new = new
    df$relALBChange=relChange
    df$cell<-rownames(df)

  #add meta data to each seurat object
  meta_cell<-data.frame(cell=colnames(d10x), Sample_ID=sampl)
  meta_cell_add<-merge(meta_cell, meta, by.x="Sample_ID", by.y="Sample")
  # meta_cell_add$relALBChange<-NA
  meta_cell_add<-meta_cell_add[match(colnames(d10x), meta_cell_add$cell),]
  print(identical(meta_cell_add$cell, colnames(d10x)))
  rownames(meta_cell_add)<-meta_cell_add$cell
  d10x<- AddMetaData(d10x, meta_cell_add)
  d10x
})

d10x.list



## cell counts
plt_count_raw<-lapply(1:length(d10x.list), function(x) {
df<-data.frame(raw_cell_count=nrow(d10x.list[[x]]@meta.data),Sample_ID=unique(d10x.list[[x]]@meta.data$Sample_ID))
df})
plt_count_raw<-do.call(rbind, plt_count_raw)
print(plt_count_raw)



#'## QC
#'The percentage of reads that map to the mitochondrial genome
#'Low-quality / dying cells often exhibit extensive mitochondrial contamination
#'We calculate mitochondrial QC metrics with the PercentageFeatureSet function, which calculates the percentage of counts originating from a set of features
#'We use the set of all genes starting with MT- as a set of mitochondrial genes

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
invisible(lapply(1:length(d10x.list), function(x){
d10x.list[[x]][["percent.mt"]] <<- PercentageFeatureSet(d10x.list[[x]], pattern = "^mt-")}))

# Show QC metrics for the first 5 cells
print(head(d10x.list[[2]]@meta.data, 5))

#'Low-quality cells or empty droplets will often have very few genes
#'Cell doublets or multiplets may exhibit an aberrantly high gene count


# Visualize QC metrics
#nFeature number of unique genes
#nCount number of total molecules
plt_QC_data<-do.call(rbind, lapply(1:length(d10x.list), function(x) d10x.list[[x]]@meta.data))




qc_plts<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +  
  geom_point() +
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th
save_plts(qc_plts, "FRP_intital_QC_plts", w=6,h=4)

qc_plts_chem<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
  geom_point() + facet_wrap(~Replicate)+
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th
save_plts(qc_plts_chem, "FRP_intital_QC_plts_Replicate", w=12,h=4)

qc_plts_Sample_ID<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
  geom_point() + facet_wrap(~Sample_ID)+
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th
save_plts(qc_plts_chem, "FRP_intital_QC_plts_Sample_ID", w=12,h=4)

MT_plt<-ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
  geom_vline(xintercept = 25)+ theme_bw()+xlab("Percent Mitochondrial")+th
save_plts(MT_plt, "FRP_percentMT_plt", w=6,h=4)

MT_plt_Sample_ID<-ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
  facet_wrap(~Sample_ID, scales="free_y")+
  geom_vline(xintercept = 25)+ theme_bw()+xlab("Percent Mitochondrial")+th
save_plts(MT_plt_Sample_ID, "FRP_percentMT_plt_Sample_ID", w=8,h=4)



#'We filter cells that have unique feature counts over 6,000 or less than 500
#'We filter cells that have >10% mitochondrial counts
#'we will also filter doublets as called by scrublet
invisible(lapply(1:length(d10x.list), function(x){
  d10x.list[[x]] <<- subset(d10x.list[[x]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 25)
  }))

d10x.list

## cell counts after QC
plt_count_QC<-lapply(1:length(d10x.list), function(x) {
  df<-data.frame(qc_cell_count=nrow(d10x.list[[x]]@meta.data),Sample_ID=unique(d10x.list[[x]]@meta.data$Sample_ID))
  df})
  plt_count_QC<-do.call(rbind, plt_count_QC)
print(plt_count_QC)

counts<-merge(plt_count_raw, plt_count_QC, by="Sample_ID")
meta<-merge(meta,counts,by.x="Sample", by.y="Sample_ID")

cell_count<-grid.arrange(ggplot(meta, aes(Type, raw_cell_count,fill=Type))+
                         geom_boxplot()+geom_point()+
                         theme_bw()+geom_text(aes(label=Sample), hjust=1.05, size=3)+xlab("Age Group")+
                         ylab("Total Cell Number")+th+fillscale_age+ylim(0,60000)+
                         theme(legend.position = "none")+ggtitle("Before Quality Control"),
                       ggplot(meta, aes(Type, qc_cell_count,fill=Type))+
                         geom_boxplot()+geom_point()+
                         theme_bw()+geom_text(aes(label=Sample), hjust=1.05, size=3)+xlab("Age Group")+
                         ylab("Total Cell Number")+th+fillscale_age+ylim(0,60000)+
                         theme(legend.position = "none")+ggtitle("After Quality Control"), ncol=2)

save_plts(cell_count, "FRP_QC_cellcount_age", w=8,h=4)


d10x <- merge(d10x.list[[1]], y= d10x.list[c(2,3,5)], merge.data=TRUE, project = "FRP")#add.cell.ids = alldata_names2,

d10x



save(d10x.list, file=here("/media/redgar/Seagate Portable Drive/FRP_mito/raw_FRP_for_integration.RData"))


load("/media/redgar/Seagate Portable Drive/FRP_mito/raw_FRP_for_integration.RData")
################
## Normalize scale and UMAP
################
d10x <- NormalizeData(d10x)
d10x <- FindVariableFeatures(d10x, selection.method = "vst", nfeatures = 2000)
d10x <- ScaleData(d10x) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))

# dimension reduction
d10x <- RunPCA(d10x, ndims.print = 1:10, nfeatures.print = 10)
d10x <- RunUMAP(d10x, dims = 1:30)

merged_cluster_umap<-DimPlot(d10x, reduction = "umap", pt.size=0.25, group.by = "Sample_ID")
save_plts(merged_cluster_umap, "FRP_merged_umap_Sample_ID", w=6,h=4)

merged_cluster_umap<-DimPlot(d10x, reduction = "umap", pt.size=0.25, group.by = "Condition")
save_plts(merged_cluster_umap, "FRP_merged_umap_Condition", w=6,h=4)

######################
## cell cycle gene expression
######################
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes



###############
## Integration
###############
load("/cluster/projects/macparland/RE/FRP_mito/raw_FRP_for_integration.RData")


#https://satijalab.org/seurat/articles/integration_rpca.html
print("RUNNING INTEGRATION")

## run integration across donor and hopefully that will also smooth out differences with chemistry?
## data is already split by donor

# normalize, identify variable features and score cell cycle for each dataset independently
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

d10x.list <- lapply(X = d10x.list[c(1,2,3,5)], FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = d10x.list)
d10x.list <- lapply(X = d10x.list, FUN = function(x) {
  #x <- ScaleData(x, features = features, verbose = FALSE)
  x <- ScaleData(x, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})



## Identify anchors
chem.anchors <- FindIntegrationAnchors(object.list = d10x.list, anchor.features = features, reduction = "rpca")
d10x.combined <- IntegrateData(anchorset = chem.anchors)

DefaultAssay(d10x.combined) <- "integrated"

print("INTEGRATED")


# Run the standard workflow for visualization and clustering
d10x.combined <- ScaleData(d10x.combined, verbose = FALSE)
d10x.combined <- RunPCA(d10x.combined, npcs = 30, verbose = FALSE)
d10x.combined <- RunUMAP(d10x.combined, reduction = "pca", dims = 1:30)
d10x.combined <- RunTSNE(d10x.combined, dims = 1:30)

d10x.combined <- FindNeighbors(d10x.combined, reduction = "pca", dims = 1:30)
d10x.combined <- FindClusters(d10x.combined, resolution = 0.5)

d10x.combined

save(d10x.combined, file="/cluster/projects/macparland/RE/FRP_mito/FRP_integrated.RData")








#load("/media/redgar/Seagate Portable Drive/FRP_mito/FRP_integrated.RData")
load("../FRP_mito/FRP_integrated.RData")


##########
# Visualize integration
##########
SCT_cluster_umap<-DimPlot(d10x.combined, reduction = "umap", pt.size=0.25, label=T)
save_plts(SCT_cluster_umap, "FRP_rPCA_cluster_umap", w=6,h=4)

chem_umap_sct<-DimPlot(d10x.combined, reduction = "umap", group.by = "Sample_ID", pt.size=0.25)
save_plts(chem_umap_sct, "FRP_Sample_ID_rPCA_umap", w=6,h=4)

Sample_ID_split<-DimPlot(d10x.combined, reduction = "umap", group.by = "Sample_ID", split.by="Sample_ID",pt.size=0.25)
save_plts(Sample_ID_split, "FRP_Sample_ID_facet_rPCA_UMAP", w=15,h=4)

age_umap_sct<-DimPlot(d10x.combined, reduction = "umap", group.by = "Condition", pt.size=0.25)+scale_color_manual(values=c("grey80","#20637b"))
save_plts(age_umap_sct, "FRP_Condition_rPCA_umap", w=6,h=4)

age_umap_sct<-DimPlot(d10x.combined, reduction = "umap", group.by = "Type", pt.size=0.25)+scale_color_manual(values=c("goldenrod","cornflowerblue"))
save_plts(age_umap_sct, "FRP_Type_rPCA_umap", w=6,h=4)

MT_umap_sct<-FeaturePlot(d10x.combined, reduction = "umap", features = "percent.mt", pt.size=0.25)
save_plts(MT_umap_sct, "FRP_MT_rPCA_umap", w=5,h=4)

alb_umap_sct<-FeaturePlot(d10x.combined, reduction = "umap", features = "Alb", pt.size=0.25)
save_plts(alb_umap_sct, "FRP_alb_rPCA_umap", w=5,h=4)





############################
##### Cell Annotation
############################
Macrophage_genes<-c( "Clec4f", "Vsig4", "Lyz1","CCr2")
neutro<-c("Retnig","Ccl4")
B_genes<-c("Cd79a","Ms4a1")
NK_T_genes<-c("Cd3g","Nkg7","Xcl1","Ncr1")
LEC_genes<-c("Clec4g","Fabp4")
Hepatocyte_genes<-c("Alb", "Serpina3k")
Cholangiocytes_genes<-c( "Epcam", "Tspan8")
HSCs_genes<-c( "Rein","Rgs5")
T_genes<-c("Cd3g","Nkg7")
NK_genes<-c("Xcl1","Ncr1")


########
## manual cell labels
########
macrophage<-c("Csf1r", "C1qb", "Clec4f", "Apoe", "Lyz2", "Ctss","Ly6c1")
kc_genes<-c( "Cd5l", "Timd4", "C1qa", "Vsig4", "Fcna", "C1qc","Adgre1")
neutro<-c("Retnig","Ccl4")
B_genes<-c("Cd79a","Ms4a1")
LEC_genes<-c("Clec4g","Fabp4")
Hepatocyte_genes<-c("Alb", "Serpina3k")
Cholangiocytes_genes<-c( "Epcam", "Tspan8")
HSCs_genes<-c( "Rein","Rgs5")
T_genes<-c("Cd3g","Nkg7")
NK_genes<-c("Xcl1","Ncr1")
zonation_markers_CV<-c("Cyp2e1","Glul","Axin2","Cyp1a2","Gstm3","Psmd4")
zonation_markers_PN<-c("Ass1","Asl","Cyp2f2","Arg1","Pck1","C2","Sdhd")



DefaultAssay(d10x.combined) <- "RNA"

pdf(file = here("../FRP_mito/figures/FRP_dot_plots.pdf"), w=10, h=10)
DotPlot(object = d10x.combined, features = B_genes)+xlab("B Cell Marker")
DotPlot(object = d10x.combined, features = T_genes)+xlab("T Cell Marker")
DotPlot(object = d10x.combined, features = NK_genes)+xlab("NK Cell Marker")
DotPlot(object = d10x.combined, features = LEC_genes)+xlab("LSEC Marker")
DotPlot(object = d10x.combined, features = Hepatocyte_genes)+xlab("Hepatocyte Marker")
DotPlot(object = d10x.combined, features = zonation_markers_CV)+xlab("hep central Marker")
DotPlot(object = d10x.combined, features = zonation_markers_PN)+xlab("hep portal Marker")
DotPlot(object = d10x.combined, features = Cholangiocytes_genes)+xlab("Cholangiocyte Marker")
DotPlot(object = d10x.combined, features = HSCs_genes)+xlab("HSC Marker")
DotPlot(object = d10x.combined, features = neutro)+xlab("Neutro Marker")
DotPlot(object = d10x.combined, features = macrophage)+xlab("mono-mac Marker")
DotPlot(object = d10x.combined, features = kc_genes)+xlab("kc Marker")
dev.off()


# ## donor integrated map, clusters needing checks: 6,9
mt_percent_box<-ggplot(d10x.combined@meta.data, aes(seurat_clusters, percent.mt))+geom_violin( fill='lightgrey')+theme_bw()
save_plts(mt_percent_box, "FRP_mt_percent_box", w=12,h=2)

# 12 Highly ALB contaminated myeloid
hep_QC<-FeaturePlot(d10x.combined, reduction = "umap",features="Alb",pt.size=0.15)
hep_QC
save_plts(hep_QC, "FRP_ALB_rPCA_UMAP", w=6,h=5)

# 21 seems to be red blood cells
RBC<-FeaturePlot(d10x.combined, features = c("Hbb-bt","Hba-a1","Hbb-bh2","Hba-a2"), min.cutoff = "q9", pt.size=0.15)
RBC
save_plts(RBC, "FRP_RBC_rPCA_UMAP", w=8,h=7)

DimPlot(d10x.combined, reduction = "umap", group.by = "Phase")
FeaturePlot(d10x.combined, features = c("Mki67", "Top2a"), min.cutoff = "q9", pt.size=0.15)


DimPlot(d10x.combined, reduction = "umap")+scale_color_manual(values=c(rep("grey",34),"red",rep("grey", 19)))



FeaturePlot(d10x.combined, features = c("Igfbp7", "Aqp1", "Bmp2", "Stab2", "Sparc", "Clec4g", "Kdr", "Plpp3", "Selenop"), min.cutoff = "q9", pt.size=0.15)
FeaturePlot(d10x.combined, features = c("Acta2", "Col3a1", "Dcn", "Tpm2", "Col14a1", "Ecm1", "Tmem56", "Bgn"), min.cutoff = "q9", pt.size=0.15)




rownames(d10x.combined)[grep("Serp", rownames(d10x.combined))]


cluster6.markers <-  FindMarkers(d10x.combined, ident.1 = 6, min.pct = 0.25)
head(cluster6.markers[which(cluster6.markers$avg_log2FC>0),], n = 20)
FeaturePlot(d10x.combined, reduction = "umap", features = c("Ccnl1","Aox3","Asl"), ncol = 2)

cluster20.markers <-  FindMarkers(d10x.combined, ident.1 = 20, min.pct = 0.25)
head(cluster20.markers[which(cluster20.markers$avg_log2FC>0),], n = 20)
FeaturePlot(d10x.combined, features = c("Pou2f2", "Sh2d1b1", "Cd244a"), min.cutoff = "q9", pt.size=0.15)
FeaturePlot(d10x.combined, features = c("Cd79b", "Cd19", "Ighm", "Ighd", "Ms4a1", "Ly6d", "H2-Ob"), min.cutoff = "q9", pt.size=0.15)





#############
# #relabel some clusters
#############
DefaultAssay(d10x.combined) <- "integrated"
d10x.combined@meta.data$CellType<-as.character(d10x.combined@meta.data$seurat_clusters)
d10x.combined@meta.data$CellType[which(d10x.combined@meta.data$seurat_clusters=="12")]<-"B-cells"
d10x.combined@meta.data$CellType[which(d10x.combined@meta.data$seurat_clusters=="10")]<-"CD3+ T-cells"
d10x.combined@meta.data$CellType[which(d10x.combined@meta.data$seurat_clusters=="17")]<-"NK-like cells"
d10x.combined@meta.data$CellType[which(d10x.combined@meta.data$seurat_clusters%in%c("5","20"))]<-"LSEC"
d10x.combined@meta.data$CellType[which(d10x.combined@meta.data$seurat_clusters%in%c("16","13"))]<-"Hepatocytes (Central)"
d10x.combined@meta.data$CellType[which(d10x.combined@meta.data$seurat_clusters%in%c("2","9"))]<-"Hepatocytes (Portal)"
d10x.combined@meta.data$CellType[which(d10x.combined@meta.data$seurat_clusters%in%c("1","3","0","4","18","6"))]<-"Hepatocytes"
d10x.combined@meta.data$CellType[which(d10x.combined@meta.data$seurat_clusters=="22")]<-"Cholangiocytes"
d10x.combined@meta.data$CellType[which(d10x.combined@meta.data$seurat_clusters%in%c("8"))]<-"HSC"
d10x.combined@meta.data$CellType[which(d10x.combined@meta.data$seurat_clusters%in%c("14","21"))]<-"Neutrophil"
d10x.combined@meta.data$CellType[which(d10x.combined@meta.data$seurat_clusters%in%c("15"))]<-"Macrophage"
d10x.combined@meta.data$CellType[which(d10x.combined@meta.data$seurat_clusters%in%c("19","11","7"))]<-"KC Like"
d10x.combined@meta.data$CellType[which(d10x.combined@meta.data$seurat_clusters%in%c("23"))]<-"Erythrocytes"


DimPlot(d10x.combined, reduction = "umap",group.by="CellType", pt.size=0.15, label=T)+colscale_cellType+ggtitle("")+
  annotate("text", x=-9, y=-14, label = paste0("n = ",comma(ncol(d10x.combined))))

########
## mean genes per cell
########
d10x.combined@meta.data %>%
  group_by(Sample_ID) %>%
  summarise(max = mean(nFeature_RNA, na.rm=TRUE))

#####
#fancy dot
#####
DefaultAssay(d10x.combined) <- "RNA"

gene_exp<-FetchData(d10x.combined, vars=c(macrophage,kc_genes,neutro, 
                                          B_genes,T_genes, NK_genes,
                                          LEC_genes,Hepatocyte_genes,zonation_markers_CV,zonation_markers_PN,
                                          Cholangiocytes_genes,HSCs_genes,c("Hbb-bt","Hba-a2")))
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)

plt<-merge(gene_exp, d10x.combined@meta.data, by.x="cell", by.y="cell")

## summarize
scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}

plt_summary<-plt %>% group_by(CellType, variable) %>% 
  summarise(mn=mean(value), count=length(which(value>0)), percent_exp=(length(which(value>0))/length(value))*100)
plt_summary <- plt_summary %>% group_by(variable) %>%
  dplyr::mutate(scaled = scale_this(mn))
plt_summary<-as.data.frame(plt_summary)

# remove dots where 0 cell expressing marker
plt_summary<-plt_summary[(which(plt_summary$count>0)),]

plt_summary$variable<-factor(plt_summary$variable, levels=rev(c(macrophage,kc_genes,neutro, 
                                                                B_genes,T_genes, NK_genes,
                                                                LEC_genes,Hepatocyte_genes,zonation_markers_CV,zonation_markers_PN,
                                                                Cholangiocytes_genes,HSCs_genes,c("Hbb-bt","Hba-a2"))))

plt_summary$CellType<-factor(plt_summary$CellType, levels=c("Macrophage","KC Like","Neutrophil","B-cells","CD3+ T-cells","NK-like cells",
                                                            "LSEC", "Hepatocytes", "Hepatocytes (Central)","Hepatocytes (Portal)",
                                                            "Cholangiocytes","HSC","Erythrocytes"))

fancy_dotplot<-plot_grid(
  ggplot(plt_summary[which(!(is.na(plt_summary$variable))),], aes(CellType, variable, color=scaled, size=percent_exp))+geom_point()+
    th+theme_classic()+
    scale_color_continuous_sequential(palette = "Oslo", rev=F, name="Scaled\nMean\nExpression")+
    scale_size(name="Percent\nCells\nExpressing")+
    theme(axis.text.x = element_blank(),axis.title = element_blank(),axis.ticks.x = element_blank(),
          plot.margin = margin(0.25,0.25,0,0.25,"cm"))+
    geom_hline(yintercept = 36.5)+ geom_hline(yintercept = 26.5)+ 
    geom_hline(yintercept = 29.5)+ geom_hline(yintercept = 28.5)+ 
    geom_hline(yintercept = 22.5)+geom_hline(yintercept = 20.5)+
    geom_hline(yintercept = 18.5)+  geom_hline(yintercept = 12.5)+
    geom_hline(yintercept = 5.5)+ geom_hline(yintercept = 2.5) +geom_hline(yintercept = 3.5)  ,
  ggplot(plt_summary, aes(CellType, y=1, fill=CellType))+geom_tile(color="black")+
    th+theme_classic()+fillscale_cellType+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",axis.line  = element_blank(),
          plot.margin = margin(0,0,1,1,"cm")),
  ncol=1, rel_heights = c(6,1.3), align = "v", axis="lr")
fancy_dotplot
save_plts(fancy_dotplot, "marker_dotplot", h=12,w=6)


##############
## KC marcophage differential genes
##############
DefaultAssay(d10x.combined) <- "RNA"
Idents(d10x.combined)<-"CellType"

de<-FindMarkers(d10x.combined, ident.1 = "Macrophage", ident.2 = "KC Like", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F, logfc.threshold=0)
KC_de<-de[which(de$p_val_adj < 0.005 & de$avg_log2FC < 0),]
mono_mac_de<-de[which(de$p_val_adj < 0.005 & de$avg_log2FC > 0),]
  
FeaturePlot(d10x.combined, features = "Marco")

all_refined_cluster_umap<-DimPlot(d10x.combined, reduction = "umap", pt.size=0.25, label=T,label.size = 3, group.by = "CellType")+
  colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-13, y=-13, label=paste0("n = ",comma(ncol(d10x.combined))))
all_refined_cluster_umap
save_plts(all_refined_cluster_umap, "FRP_refined_cellType_map", w=12,h=8)


all_refined_cluster_umap_nolab<-DimPlot(d10x.combined, reduction = "umap", pt.size=0.25, group.by = "CellType")+colscale_cellType+ggtitle("")+xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-13, y=-13, label=paste0("n = ",comma(ncol(d10x.combined))))
all_refined_cluster_umap_nolab
save_plts(all_refined_cluster_umap_nolab, "FRP_refined_cellType_map_nolabel", w=12,h=8)

Sample_ID_split<-DimPlot(d10x.combined, reduction = "umap", group.by = "CellType", split.by="Sample_ID",pt.size=0.25, ncol=2)+colscale_cellType+ggtitle("")
save_plts(Sample_ID_split, "FRP_Sample_ID_roughCell_facet_rPCA_UMAP_refined", w=14,h=14)



#############
## Cell type count table
#############
count_plt<-as.data.frame(d10x.combined@meta.data %>% dplyr::select(CellType,Sample_ID) %>% group_by(Sample_ID) %>% count(CellType))

count_plt$CellType_rough<-as.factor(count_plt$CellType)
levels(count_plt$CellType_rough)<-c("B-cells","NK and T cells","Cholangiocytes","Erythrocytes","Hepatocytes","Hepatocytes", "Hepatocytes" ,
                                    "HSC",    "Myeloid cells","LSEC",   "Myeloid cells","Myeloid cells","NK and T cells" )

count_plt$CellType_rough<-factor(count_plt$CellType_rough, levels=c(
  "Hepatocytes","LSEC","HSC","Cholangiocytes",
  "Myeloid cells","NK and T cells",  "B-cells",
  "Erythrocytes"))

cell_count<-ggplot(count_plt, aes(CellType_rough,n,  fill=CellType))+
  geom_bar(stat="identity", color="black")+
  facet_wrap(~Sample_ID)+fillscale_cellType+th_present+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position="bottom")+ylab("Cell Count")+xlab("Broad Cell Category")
save_plts(cell_count, "cell_count_age_condiditon_bar", w=9,h=7)






##############
## Differential Expression
##############
DefaultAssay(d10x.combined) <- "RNA"


celltypes<-c("Hepatocytes","HSC","LSEC","KC Like","CD3+ T-cells","Macrophage","B-cells", "Hepatocytes (Portal)", 
                 "Neutrophil","Hepatocytes (Central)","Erythrocytes","NK-like cells","Cholangiocytes" )

source("scripts/00_GSEA_function.R")
#https://download.baderlab.org/EM_Genesets/current_release/Mouse/symbol/
GO_file = here("../FRP_mito/Mouse_GOBP_AllPathways_with_GO_iea_January_01_2024_symbol.gmt")

#HEP and mon mac

lapply(1:length(celltypes), function(x){
  print(celltypes[x])
  d10x_subset<-subset(d10x.combined, subset = CellType %in% celltypes[x])
  
  Idents(d10x_subset)<-d10x_subset$Condition
  table(d10x_subset$Condition)
  
  ## differential
  print("DGE")
  de<-FindMarkers(d10x_subset, ident.1 = "Control", ident.2 = "Treated", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F, logfc.threshold=0)
  sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]
  print(nrow(sig_de))
  
  write.csv(sig_de, file=paste("../FRP_mito/DGE/DEG_",celltypes[x],".csv", sep=""))
  
  ### GSEA
  print("GSEA")
  de$gene<-rownames(de)
  gene_list = de$avg_log2FC
  names(gene_list) = de$gene
  gene_list = sort(gene_list, decreasing = TRUE)
  gene_list = gene_list[!duplicated(names(gene_list))]
  
  res = GSEA(gene_list, GO_file, pval = 0.05)
  
  plt_path<-res$Results
  plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
  plt_path$Enrichment_Cell<-"Down-regulated in \nTreated"
  plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \nTreated"
  
  plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))
  
  plt_path$direction_label<-as.factor(plt_path$Enrichment)
  levels(plt_path$direction_label)<-c(0.1,-0.1)
  plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))
  
  # top and bottom 15
  plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])
  
  celltype_GSEA<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
    theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
    geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
    scale_fill_manual(values=c("#D64A56","cornflowerblue"))
  celltype_GSEA
  save_plts(celltype_GSEA, paste("GSEA_",celltypes[x],sep=""), w=15,h=7)
  gc()
})


## plot genes

plt_gene<-function(gene){
  gene_exp<-FetchData(d10x.combined, vars=gene)
  gene_exp$cell<-rownames(gene_exp)
  plt_box<-merge(d10x.combined@meta.data, gene_exp, by='cell')
  
  plt_box$Condition<-as.factor(plt_box$Condition)
  levels(plt_box$Condition)<-c("Ischemia\nreperfusion","Ischemia\nreperfusion\nplus\nmitochondria\ntransplant")
  
 box<- ggplot(plt_box, aes(Condition, eval(parse(text = gene))))+geom_violin(fill="lightgrey", color="lightgrey")+
    geom_boxplot(aes(fill=Condition),width=0.1, outlier.shape = NA)+facet_wrap(~CellType)+theme_bw()+scale_fill_manual(values=c("#3094b8","#b83050"))+
    theme(legend.position = "none")+  ylim(0,max(plt_box[,18])+1)+ylab(paste(gene, "Expression"))
  save_plts(box, paste(gene,"_all_cell_types", sep=""), w=10, h=10)
  
  Cyp2e1_umap<-FeaturePlot(d10x.combined, features = gene, split.by = "Condition")
  save_plts(Cyp2e1_umap, paste(gene,"_UMAP", sep=""), w=9, h=4)
}

plt_gene("Cd38")
plt_gene("Cd40")
plt_gene("Cd163")
#plt_gene("Cd206")
plt_gene("Arg1")
plt_gene("Nos2")
plt_gene("Il6")
plt_gene("Tnf")
plt_gene("Il10")
plt_gene("Hgf")
plt_gene("Vsig4")
plt_gene("Mki67")
plt_gene("Top2a")
plt_gene("Ccna2")
plt_gene("Cdkn1a")
plt_gene("Chek1")
plt_gene("Mcm5")

plt_gene("Il6")
plt_gene("Il2")
plt_gene("Il3")
plt_gene("Il4")
plt_gene("Il5")
plt_gene("Ifng")
plt_gene("Egr1")

plt_gene("Egr2")


plt_gene("Apoa4")
plt_gene("Nfkb1")






Cyp2e1<-VlnPlot(d10x.combined, features = "Cyp2e1", split.by = "Condition", group.by = "CellType")+scale_fill_manual(values=c("grey80","#20637b"))
save_plts(Cyp2e1, "Cyp2e1_all_cell_types", w=10, h=4)

Cyp2e1_umap<-FeaturePlot(d10x.combined, features = "Cyp2e1", split.by = "Condition")
save_plts(Cyp2e1_umap, "Cyp2e1_UMAP", w=9, h=4)

Vsig4_umap<-FeaturePlot(d10x.combined, features = "Vsig4", split.by = "Condition")
save_plts(Vsig4_umap, "Vsig4_UMAP", w=9, h=4)



##############
## Cell Cycle
##############

head(d10x.combined@meta.data)

phase_propotion<- d10x.combined@meta.data %>% group_by(CellType, Condition, Phase) %>% summarise(Percentage = n() / nrow(iris) * 100)

phase_propotion<- d10x.combined@meta.data %>%  group_by(CellType, Condition,Phase) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))



cell_cycle<-ggplot(phase_propotion, aes(Condition, freq, fill=Phase))+
  geom_bar(stat="identity", color="black")+
  facet_wrap(~CellType, scale="free_x")+th_present+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position="bottom")+ylab("Cell Count")+xlab("Cell Cycle Phase")+
  scale_fill_manual(values=c("#1c9099","#016c59","grey"))
save_plts(cell_cycle, "cell_cycle_mito", h=8,w=5)






###########
## Volcano
###########
d10x_subset<-subset(d10x.combined, subset = CellType %in% celltypes[4])

Idents(d10x_subset)<-d10x_subset$Condition
table(d10x_subset$Condition)

## differential
print("DGE")
de<-FindMarkers(d10x_subset, ident.1 = "Control", ident.2 = "Treated", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F, logfc.threshold=0)



volcano<-data.frame(gene=rownames(de),Pvalue=de$p_val, Delta_Beta=de$avg_log2FC)
  
#Thresholds
dB<-0.5 #delta beta cutoff
Pv<- 0.005 #Pvalue cutoff
 
volcano<-volcano[complete.cases(volcano),]

## positive delta beta is hypomethylated (code for volcano should be right now, should colors change?)
color3<-sapply(1:nrow(volcano), function(x) if(volcano$Pvalue[x]<=Pv){
  if(abs(volcano$Delta_Beta[x])>dB){
    if(volcano$Delta_Beta[x]>dB){"Higher Expression in Controls\n(above FC)"}else{"Higher Expression in Treated\n(above FC)"}
  }else{if(volcano$Delta_Beta[x]>0){"Higher Expression in Controls"}else{"Higher Expression in Treated"}}}else{"Not Significantly Different"})

volcano$Interesting_CpG3<-color3


# COLORS! define here so they are consistent between plots
# so even if you don't have CpGs in a color catagory the pattern will be maintained
myColors <- c(muted("red", l=80, c=30),"red",muted("blue", l=70, c=40),"blue", "grey")

color_possibilities<-c("Higher Expression in Controls",
                       "Higher Expression in Controls\n(above FC)",
                       "Higher Expression in Treated",
                       "Higher Expression in Treated\n(above FC)",
                       "Not Significantly Different")

names(myColors) <- color_possibilities
colscale <- scale_color_manual(name = "Direction of Change",
                               values = myColors, drop = FALSE)
fillscale <- scale_fill_manual(name = "Direction of Change",
                               values = myColors, drop = FALSE)


sig_MC_celltype<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 0.8),]

volcano_label<-volcano[which(volcano$gene%in%rownames(sig_MC_celltype)),]
volcano$sig<-"not_sig"
volcano$sig[which(volcano$gene%in%rownames(sig_MC_celltype))]<-"MC_sig"

vol<-ggplot(volcano, aes(Delta_Beta, -log10(Pvalue), fill=Interesting_CpG3, color=sig))+
  geom_point(shape=21, size=1)+theme_bw()+
  fillscale+scale_color_manual(values=c("black","white"))+guides(color = "none")+
  geom_vline(xintercept=c(-dB,dB), color="grey60")+
  geom_hline(yintercept=-log10(Pv), color="grey60")+
  ylab("P Value (-log10)")+xlab("Differential Expression\n(Fold change)")+
  theme(plot.margin=unit(c(1,1,1,2),"cm"))+
  guides(fill = guide_legend(override.aes = list(size = 4)))+
  geom_text(aes(label=gene),volcano_label,color="black",vjust=1.1, hjust=1,size=3)

save_plts(vol, "KC_volcano_relaxed_Fold_change", w=10,h=7)


sig_MC_celltype<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 0.5),]
sig_MC_celltype$direction<-"Higher in Treated"
sig_MC_celltype$direction[which(sig_MC_celltype$avg_log2FC >0 )]<-"Higher in Controls"

write.csv(sig_MC_celltype, file="../FRP_mito/DGE/KC_cell_differential_expression_relaxed_fold_change.csv")



############
## Scores
############
d10x_subset<-subset(d10x.combined, subset = CellType %in% celltypes[c(4,6)])
DefaultAssay(d10x_subset)<-"RNA"
# 
# a. M1-like Score = TNF-alpha, IL-6, CD38, CD40, Nos2, IL1-beta, Fpr2
# b. M2-like Score = Arg1, Egr2, CD163, CD206, c-Myc, IL-10, Fn1
# c. Pro-efferocytosis/pro-resolution Score = MerTK, TIM4 (also known as TIMD4), Axl, TYRO, PPAR-gamma, 

M1 <- list(c(
  "Tnf",  "Il6",  "Cd38",  "Cd40",  "Nos2",  "Il1b",  "Fpr2"))
M2 <- list(c(
   "Cd163",  "Mrc1",  "Arg1",  "Egr2",  "Myc",  "Fn1",  "Il10"))
Pro_efferocytosis <- list(c(
  "Mertk",  "Timd4",  "Axl",  "Tyrobp",  "Pparg"))


d10x_subset <- AddModuleScore(
    object = d10x_subset,
    features = M1,
    ctrl = 5,
    name = 'M1_score'
  )


d10x_subset <- AddModuleScore(
  object = d10x_subset,
  features = M2,
  ctrl = 5,
  name = 'M2_score'
)

d10x_subset <- AddModuleScore(
  object = d10x_subset,
  features = Pro_efferocytosis,
  ctrl = 5,
  name = 'pro_efferocytosis_score'
)


d10x_subset@meta.data$Condition<-as.factor(d10x_subset@meta.data$Condition)
levels(d10x_subset@meta.data$Condition)<-c("Ischemia\nreperfusion","Ischemia\nreperfusion\nplus\nmitochondria\ntransplant")

M1_score<-ggplot(d10x_subset@meta.data, aes(Condition, M1_score1))+geom_violin(fill="lightgrey", color="lightgrey")+
  geom_boxplot(aes(fill=Condition),width=0.25, outlier.shape = NA)+facet_wrap(~CellType)+theme_bw()+scale_fill_manual(values=c("#3094b8","#b83050"))+ ylim(-1,2.25)+
  ylab("M1-like Score")+ theme(legend.position = "none")+  
  geom_signif(comparisons = list(c("Ischemia\nreperfusion","Ischemia\nreperfusion\nplus\nmitochondria\ntransplant")),map_signif_level = TRUE, color="grey40")
M1_score
save_plts(M1_score, "M1_score_KC_macrophage", w=4,h=4)

M2_score<-ggplot(d10x_subset@meta.data, aes(Condition, M2_score1))+geom_violin(fill="lightgrey", color="lightgrey")+
  geom_boxplot(aes(fill=Condition),width=0.25, outlier.shape = NA)+facet_wrap(~CellType)+theme_bw()+scale_fill_manual(values=c("#3094b8","#b83050"))+ ylim(-1,2.25)+
  ylab("M2-like Score")+ theme(legend.position = "none")+  
  geom_signif(comparisons = list(c("Ischemia\nreperfusion","Ischemia\nreperfusion\nplus\nmitochondria\ntransplant")),map_signif_level = TRUE, color="grey40")
M2_score
save_plts(M2_score, "M2_score_KC_macrophage", w=4,h=4)

pro_efferocytosis_score<-ggplot(d10x_subset@meta.data, aes(Condition, pro_efferocytosis_score1))+geom_violin(fill="lightgrey", color="lightgrey")+
  geom_boxplot(aes(fill=Condition),width=0.25, outlier.shape = NA)+facet_wrap(~CellType)+theme_bw()+scale_fill_manual(values=c("#3094b8","#b83050"))+
  ylab("Pro-efferocytosis or pro-resolution Score")+ theme(legend.position = "none")+  ylim(-1,2.75)+
  geom_signif(comparisons = list(c("Ischemia\nreperfusion","Ischemia\nreperfusion\nplus\nmitochondria\ntransplant")),map_signif_level = TRUE, color="grey40")
pro_efferocytosis_score
save_plts(pro_efferocytosis_score, "pro_efferocytosis_score_KC_macrophage", w=4,h=4)

mt_percent<-ggplot(d10x.combined@meta.data, aes(Condition, percent.mt))+geom_violin(fill="lightgrey", color="lightgrey")+
  geom_boxplot(aes(fill=Condition),width=0.25)+facet_wrap(~CellType)+theme_bw()+scale_fill_manual(values=c("#3094b8","#b83050"))
mt_percent
save_plts(mt_percent, "mt_percent_celltype", w=7,h=12)


t.test(d10x_subset@meta.data$M1_M21[which(d10x_subset@meta.data$CellType=="KC Like")] ~ d10x_subset@meta.data$Condition[which(d10x_subset@meta.data$CellType=="KC Like")])
t.test(d10x_subset@meta.data$M1_M21[which(d10x_subset@meta.data$CellType=="Macrophage")] ~ d10x_subset@meta.data$Condition[which(d10x_subset@meta.data$CellType=="Macrophage")])


################
## M1/M2 ratio
################

d10x_subset$M1_score1_shift<-d10x_subset$M1_score1 + 1
d10x_subset$M2_score1_shift<-d10x_subset$M2_score1 + 1
d10x_subset$M2_M1_ratio_shift<-d10x_subset$M2_score1_shift/d10x_subset$M1_score1_shift


ratio<-ggplot(d10x_subset@meta.data, aes(Condition, M2_M1_ratio_shift))+geom_violin(fill="lightgrey", color="lightgrey")+
  geom_boxplot(aes(fill=Condition),width=0.25, outlier.shape = NA)+facet_wrap(~CellType)+theme_bw()+scale_fill_manual(values=c("#3094b8","#b83050"))+
  ylab("M2/M1")+ theme(legend.position = "none")+  ylim(0,4)+
geom_signif(comparisons = list(c("Ischemia\nreperfusion","Ischemia\nreperfusion\nplus\nmitochondria\ntransplant")),map_signif_level = TRUE, color="grey40")
ratio
save_plts(ratio, "polarization_M2_M1_ration", w=4,h=4)

ggplot(d10x_subset@meta.data, aes(M1_score1, M2_score1, color=M2_M1_ratio_shift))+geom_point()+facet_grid(CellType~Condition)+theme_bw()+
  scale_color_continuous_diverging(palette="Blue-Red 3", mid=1.5)





d10x_subset <- NormalizeData(d10x_subset)
d10x_subset <- FindVariableFeatures(d10x_subset, selection.method = "vst", nfeatures = 2000)
d10x_subset <- ScaleData(d10x_subset, verbose = FALSE)
d10x_subset <- RunPCA(d10x_subset, npcs = 30, verbose = FALSE)
d10x_subset <- RunUMAP(d10x_subset, reduction = "pca", dims = 1:30)

fancyUMAP_all<-fanciest_UMAP(d10x_subset,NA,F)
save_plts(fancyUMAP_all, "FRP_umap_macrophage_fancy", w=4,h=3)








## Score UMAPS

d10x_subset <- NormalizeData(d10x_subset)
d10x_subset <- FindVariableFeatures(d10x_subset, selection.method = "vst", nfeatures = 2000)
d10x_subset <- ScaleData(d10x_subset, verbose = FALSE)
d10x_subset <- RunPCA(d10x_subset, npcs = 30, verbose = FALSE)
d10x_subset <- RunUMAP(d10x_subset, reduction = "pca", dims = 1:30)

DefaultAssay(d10x_subset) <- "RNA"
umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x_subset, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x_subset@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

cell_num_all<-as.data.frame(table(plt_myeloid$Condition, plt_myeloid$CellType))
colnames(cell_num_all)<-c("Condition","CellType","CellCount")


ratio_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(color="black",size=0.5)+
  geom_point(aes(color=M2_M1_ratio_shift),size=0.2)+xlab("UMAP 1")+ylab("UMAP 2")+
  facet_grid(CellType~Condition)+
  geom_text(aes(x = 3, y = -10, label=paste0("n = ",comma(CellCount))), cell_num_all)+
  scale_color_continuous_diverging(palette="Blue-Red 3", mid=1.5, name="M2/M1 Ratio")+
  theme_bw()+theme(legend.text=element_text(size=6),
                   legend.title=element_text(size=8), 
                   plot.margin = margin(0.25,0.25,0.25,0.25, "cm"))
ratio_UMAP
save_plts(ratio_UMAP, "M2_M1_ratio_UMAP_KC_macrophage_UMAP_split", w=5.5,h=5)






## M1
exp_limit<-quantile(plt_myeloid[, "M1_score1"], 0.5)
plt_myeloid$gene_exp_limited<-NA
over_limit<-which(plt_myeloid[, "M1_score1"]>exp_limit)
plt_myeloid$gene_exp_limited[over_limit]<-plt_myeloid[over_limit, "M1_score1"]
plt_myeloid<-plt_myeloid[rev(order(plt_myeloid$gene_exp_limited)),]

plt_myeloid<-rbind(plt_myeloid[which(is.na(plt_myeloid$gene_exp_limited)),],
                   plt_myeloid[which(!(is.na(plt_myeloid$gene_exp_limited))),][(order(plt_myeloid[which(!(is.na(plt_myeloid$gene_exp_limited))),]$gene_exp_limited)),])

M1_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=gene_exp_limited),size=0.75)+xlab("UMAP 1")+ylab("UMAP 2")+
  facet_wrap(CellType~Condition,ncol=2)+
  geom_text(aes(x = 3, y = -10, label=paste0("n = ",comma(CellCount))), cell_num_all)+
  scale_color_continuous_sequential(palette = "Mako", rev=T, 
                                    name="M1-like Score",na.value = "grey80")+
  theme_bw()+theme(legend.text=element_text(size=6),
                     legend.title=element_text(size=8), 
                     plot.margin = margin(0.25,0.25,0.25,0.25, "cm"))
M1_UMAP
save_plts(M1_UMAP, "M1_UMAP_KC_macrophage_UMAP_split", w=5,h=5)



## M2
exp_limit<-quantile(plt_myeloid[, "M2_score1"], 0.5)
plt_myeloid$gene_exp_limited<-NA
over_limit<-which(plt_myeloid[, "M2_score1"]>exp_limit)
plt_myeloid$gene_exp_limited[over_limit]<-plt_myeloid[over_limit, "M1_score1"]
plt_myeloid<-plt_myeloid[rev(order(plt_myeloid$gene_exp_limited)),]

plt_myeloid<-rbind(plt_myeloid[which(is.na(plt_myeloid$gene_exp_limited)),],
                   plt_myeloid[which(!(is.na(plt_myeloid$gene_exp_limited))),][(order(plt_myeloid[which(!(is.na(plt_myeloid$gene_exp_limited))),]$gene_exp_limited)),])

M2_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=gene_exp_limited),size=0.75)+xlab("UMAP 1")+ylab("UMAP 2")+
  facet_wrap(CellType~Condition,ncol=2)+
  geom_text(aes(x = 3, y = -10, label=paste0("n = ",comma(CellCount))), cell_num_all)+
  scale_color_continuous_sequential(palette = "Mako", rev=T, 
                                    name="M2-like Score",na.value = "grey90")+
  theme_bw()+theme(legend.text=element_text(size=6),
                   legend.title=element_text(size=8), 
                   plot.margin = margin(0.25,0.25,0.25,0.25, "cm"))
save_plts(M2_UMAP, "M2_UMAP_KC_macrophage_UMAP_split", w=5,h=5)



## pro_efferocytosis_score
exp_limit<-quantile(plt_myeloid[, "pro_efferocytosis_score1"], 0.5)
plt_myeloid$gene_exp_limited<-NA
over_limit<-which(plt_myeloid[, "pro_efferocytosis_score1"]>exp_limit)
plt_myeloid$gene_exp_limited[over_limit]<-plt_myeloid[over_limit, "M1_score1"]
plt_myeloid<-plt_myeloid[rev(order(plt_myeloid$gene_exp_limited)),]

plt_myeloid<-rbind(plt_myeloid[which(is.na(plt_myeloid$gene_exp_limited)),],
                   plt_myeloid[which(!(is.na(plt_myeloid$gene_exp_limited))),][(order(plt_myeloid[which(!(is.na(plt_myeloid$gene_exp_limited))),]$gene_exp_limited)),])

effector_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(color=gene_exp_limited),size=0.75)+xlab("UMAP 1")+ylab("UMAP 2")+
  facet_wrap(CellType~Condition,ncol=2)+
  geom_text(aes(x = 3, y = -10, label=paste0("n = ",comma(CellCount))), cell_num_all)+
  scale_color_continuous_sequential(palette = "Mako", rev=T, 
                                    name="Pro-efferocytosis\npro-resolution Score",na.value = "grey90")+
  theme_bw()+theme(legend.text=element_text(size=6),
                   legend.title=element_text(size=8), 
                   plot.margin = margin(0.25,0.25,0.25,0.25, "cm"))

save_plts(effector_UMAP, "effector_UMAP_KC_macrophage_UMAP_split", w=6,h=5)


print(sessionInfo())




