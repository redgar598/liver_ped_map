                            #' #'### Load libraries
                            #' library(dplyr)
                            #' library(Seurat)
                            #' library(patchwork)
                            #' library(here)
                            #' library(ggplot2)
                            #' library(reshape2)
                            #' library(gridExtra)
                            #' library(limma)
                            #' library(cowplot)
                            #' library(gtools)
                            #' library(ggsignif)
                            #' 
                            #' 
                            #' 
                            #' #'# Functions
                            #' source(here("pretty_plots.R"))
                            #' 
                            #' #'color pallettes
                            #' myColors_copd_diagnosis <- c("grey","#659FA4","#AD7185")
                            #' color_possibilities_copd_diagnosis<-c("Control","COPD","IPF")
                            #' names(myColors_copd_diagnosis) <- color_possibilities_copd_diagnosis
                            #' fillscale_copd_diagnosis <- scale_fill_manual(name="Diagnosis", values = myColors_copd_diagnosis, drop = T)
                            #' colscale_copd_diagnosis <- scale_color_manual(name="Diagnosis",values = myColors_copd_diagnosis, drop = T)
                            #' 
                            #' 
                            #' # Downloaded
                            #' #GSE136831_AllCells.cellBarcodes.txt
                            #' #GSE136831_AllCells.GeneIDs.txt
                            #' #GSE136831_RawCounts_Sparse.mtx.gz
                            #' #GSE136831_AllCells.Samples.CellType.MetadataTable.txt
                            #' #from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136831
                            #' 
                            #' # Will rename to features.tsv and barcodes.tsv for seurat and 
                            #' # GSE136831_RawCounts_Sparse.mtx.gz to matrix.mtx.gz
                            #' 
                            #' # cd /Documents/ibd/data/public_other/scRNA/lung/GSE136831
                            #' # gzip *.tsv
                            #' 
                            #' # also removed the header from the features matrix before loading here
                            #' 
                            #' #' ## processed data?
                            #' dataset_loc <- here("../../data/public_other/scRNA/lung/GSE136831/")
                            #' 
                            #' d10x.copd <- Read10X(dataset_loc)
                            #' #' Initialize the Seurat object with the raw (non-normalized data).
                            #' d10x.copd <- CreateSeuratObject(counts = d10x.copd, min.cells = 3, min.features = 0)
                            #' d10x.copd
                            #' 
                            #' meta<-read.csv(here("../../data/public_other/scRNA/lung/GSE136831/GSE136831_AllCells.Samples.CellType.MetadataTable.txt"), sep="\t", header=T)
                            #' 
                            #' #meta_cell_add<-meta_cell_add[match(meta_cell_add$cell, colnames(d10x)),]
                            #' meta_cell_add<-meta[match(colnames(d10x.copd), meta$CellBarcode_Identity),]
                            #' rm(meta)
                            #' 
                            #' head(d10x.copd@meta.data)
                            #' 
                            #' identical(as.character(meta_cell_add$CellBarcode_Identity), rownames(d10x.copd@meta.data))
                            #' head(meta_cell_add$CellBarcode_Identity)
                            #' head(rownames(d10x.copd@meta.data))
                            #' 
                            #' rownames(meta_cell_add)<-as.character(meta_cell_add$CellBarcode_Identity)
                            #' 
                            #' d10x.copd<- AddMetaData(d10x.copd, meta_cell_add)
                            #' save(d10x.copd, file=here("../../data/public_other/scRNA/lung/GSE136831/GSE136831_raw.RData"))
                            #' 
                            #' 
                            #' ## sample cells for testing locally
                            #' set.seed(111)
                            #' d10x.copd.mini <- d10x.copd[, sample(colnames(d10x.copd), size = 1500, replace=F)]
                            #' 
                            #' save(d10x.copd.mini, file=here("../../data/public_other/scRNA/lung/GSE136831/GSE136831_mini.RData"))
                            #' 


################## COPD GSE136831 MINI DATA
load(here("../../data/public_other/scRNA/lung/GSE136831/GSE136831_mini.RData"))
d10x.copd.mini

d10x.copd<-d10x.copd.mini

#'# QC cell removal

load(here("../../data/public_other/scRNA/lung/GSE136831/GSE136831_raw.RData"))


## cell counts
plt_count_raw<-as.data.frame(tapply(rownames(d10x.copd@meta.data), d10x.copd@meta.data$Subject_Identity, length))
plt_count_raw$individual<-rownames(plt_count_raw)
colnames(plt_count_raw)<-c("cell_count","individual")

meta_disease<-d10x.copd@meta.data[,c("Subject_Identity","Disease_Identity")]
meta_disease$Subject_Identity<-as.character(meta_disease$Subject_Identity)
meta_disease$Disease_Identity<-as.character(meta_disease$Disease_Identity)
rownames(meta_disease)<-NULL
meta_disease<-meta_disease[!duplicated(meta_disease),]

plt_count_raw<-merge(plt_count_raw,meta_disease , by.x="individual",by.y="Subject_Identity")

plt_count_raw

ggplot(plt_count_raw, aes(Disease_Identity, cell_count))+
  geom_boxplot(fill="lightgrey")+geom_point()+
  theme_bw()+ylab("Total Cell Number")+th

ggsave(file=here("../../figs/scRNAseq/lung_COPD_cell_counts.pdf"), w=6,h=4)
ggsave(file=here("../../figs/scRNAseq/jpeg/lung_COPD_cell_counts.jpeg"), w=6,h=4)

#'## QC
#'The percentage of reads that map to the mitochondrial genome
#'Low-quality / dying cells often exhibit extensive mitochondrial contamination
#'We calculate mitochondrial QC metrics with the PercentageFeatureSet function, which calculates the percentage of counts originating from a set of features
#'We use the set of all genes starting with MT- as a set of mitochondrial genes

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
d10x.copd[["percent.mt"]] <- PercentageFeatureSet(d10x.copd, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(d10x.copd@meta.data, 5)

#'Low-quality cells or empty droplets will often have very few genes
#'Cell doublets or multiplets may exhibit an aberrantly high gene count


# Visualize QC metrics
#nFeature number of unique genes
#nCount number of total molecules: already filtered at min 1000 unique genes
#from geo page
#Cells filtered with > 12% of transcriptome from intron-spanning reads; < 20% mitochondrial origin; at least 1000 unique genes
ggplot(d10x.copd@meta.data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th

ggsave(file=here("../../figs/scRNAseq/lung_COPD_ncount_nfeatures.pdf"), w=6,h=4)
ggsave(file=here("../../figs/scRNAseq/jpeg/lung_COPD_ncount_nfeatures.jpeg"), w=6,h=4)

## They have also already filtered greater than 20% MT
ggplot(d10x.copd@meta.data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
  geom_vline(xintercept = 10)+ theme_bw()+xlab("Percent Mitochondrial")+th

ggsave(file=here("../../figs/scRNAseq/lung_COPD_MT_features.pdf"), w=4,h=2)
ggsave(file=here("../../figs/scRNAseq/jpeg/lung_COPD_MT_features.jpeg"), w=4,h=2)





#' #### they gave a multiple label
d10x.copd@meta.data$multiplet<-"Singlet"
d10x.copd@meta.data$multiplet[which(d10x.copd@meta.data$CellType_Category=="Multiplet")]<-"Multiplet"

ggplot(d10x.copd@meta.data, aes(nCount_RNA,nFeature_RNA,colour=multiplet)) + 
  geom_point(size=0.75) + 
  scale_color_manual(values=c("green2","black"),name="Multiplets") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th

ggsave(file=here("../../figs/scRNAseq/lung_COPD_doublet_features.pdf"), w=6,h=4)
ggsave(file=here("../../figs/scRNAseq/jpeg/lung_COPD__doublet_features.jpeg"), w=6,h=4)



#'will use the data as supplied and filtered by original authors (saved above)




##########################
#' ## MHC scoring
##########################
MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","MR1","CD1D","IRF1","NLRC5")
crypt_villis = c("SEPP1", "CEACAM7", "PLAC8", "CEACAM1", "TSPAN1", "CEACAM5", "CEACAM6", "IFI27", "DHRS9", "KRT20", "RHOC", "CD177", "PKIB", "HPGD", "LYPD8", "APOBEC1", "APOB", "APOA4", "APOA1", "NPC1L1", "EGFR", "KLF4", "ENPP3", "NT5E", "SLC28A2", "ADA")


##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x.copd <- NormalizeData(d10x.copd,scale.factor = 10000, normalization.method = "LogNormalize")

d10x.copd <- AddModuleScore(
  object = d10x.copd,
  features = list(MHCI),
  ctrl = 5,
  name = 'MHCI_score'
)

d10x.copd <- AddModuleScore(
  object = d10x.copd,
  features = list(crypt_villis),
  ctrl = 5,
  name = 'crypt_villis_score'
)



d10x.copd <- NormalizeData(d10x.copd)
d10x.copd <- FindVariableFeatures(d10x.copd, selection.method = "vst", nfeatures = 2000)
d10x.copd <- ScaleData(d10x.copd) 

# dimension reduction
d10x.copd <- RunPCA(d10x.copd, ndims.print = 1:10, nfeatures.print = 10)
d10x.copd <- RunUMAP(d10x.copd, dims = 1:30)

DimPlot(d10x.copd, reduction = "umap", group.by = "CellType_Category", pt.size=0.5)
ggsave(file=here("../../figs/scRNAseq/primary", "lung_COPD_celltypebroad_UMAP.pdf"), w=15,h=12)
ggsave(file=here("../../figs/scRNAseq/primary/jpeg", "lung_COPD_celltypebroad_UMAP.jpeg"), w=15,h=12)

DimPlot(d10x.copd, reduction = "umap", group.by = "Manuscript_Identity", pt.size=0.5, label=T)
ggsave(file=here("../../figs/scRNAseq/primary", "lung_COPD_celltype_UMAP.pdf"), w=20,h=12)
ggsave(file=here("../../figs/scRNAseq/primary/jpeg", "lung_COPD_celltype_UMAP.jpeg"), w=20,h=12)

DimPlot(d10x.copd, reduction = "umap", group.by = "Disease_Identity", pt.size=0.5)
ggsave(file=here("../../figs/scRNAseq/primary", "lung_COPD_disease_UMAP.pdf"), w=15,h=12)
ggsave(file=here("../../figs/scRNAseq/primary/jpeg", "lung_COPD_disease_UMAP.jpeg"), w=15,h=12)

DimPlot(d10x.copd, reduction = "umap", group.by = "Subject_Identity", pt.size=0.5)

FeaturePlot(d10x.copd, features = "MHCI_score1",reduction = "umap", min.cutoff = "q9", pt.size=1)
ggsave(file=here("../../figs/scRNAseq/primary", "lung_COPD_MHC1_UMAP.pdf"), w=15,h=12)
ggsave(file=here("../../figs/scRNAseq/primary/jpeg", "lung_COPD_MHC1_UMAP.jpeg"), w=15,h=12)

FeaturePlot(d10x.copd, features = "crypt_villis_score1",reduction = "umap", min.cutoff = "q9", pt.size=1)
ggsave(file=here("../../figs/scRNAseq/primary", "lung_COPD_crypt_villis_UMAP.pdf"), w=15,h=12)
ggsave(file=here("../../figs/scRNAseq/primary/jpeg", "lung_COPD_crypt_villis_UMAP.jpeg"), w=15,h=12)


### box plotting
plt<-d10x.copd@meta.data[,c("MHCI_score1","crypt_villis_score1","CellType_Category","Manuscript_Identity","Disease_Identity","Subject_Identity")]


ggplot(plt, aes(Disease_Identity,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1)+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file=here("../../figs/scRNAseq/primary", "lung_COPD_MHC1_disease.pdf"), w=15,h=6)
ggsave(file=here("../../figs/scRNAseq/primary/jpeg", "lung_COPD_MHC1_disease.jpeg"), w=15,h=6)

ggplot(plt, aes(Disease_Identity,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1)+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+facet_wrap(~CellType_Category)
ggsave(file=here("../../figs/scRNAseq/primary", "lung_COPD_MHC1_disease_broadcell.pdf"), w=10,h=6)
ggsave(file=here("../../figs/scRNAseq/primary/jpeg", "lung_COPD_MHC1_disease_broadcell.jpeg"), w=10,h=6)

ggplot(plt, aes(Disease_Identity,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(aes(fill=Disease_Identity),width=0.1)+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+facet_wrap(~Manuscript_Identity)+
  fillscale_copd_diagnosis
ggsave(file=here("../../figs/scRNAseq/primary", "lung_COPD_MHC1_disease_celltype.pdf"), w=15,h=20)
ggsave(file=here("../../figs/scRNAseq/primary/jpeg", "lung_COPD_MHC1_disease_celltype.jpeg"), w=15,h=20)




#######################
#'# Epithelial only
#######################
load(here("../../data/public_other/scRNA/lung/GSE136831/GSE136831_raw.RData"))

d10x.copd.epi<-subset(d10x.copd, subset = CellType_Category == "Epithelial")
d10x.copd.epi

MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","MR1","CD1D","IRF1","NLRC5")
crypt_villis = c("SEPP1", "CEACAM7", "PLAC8", "CEACAM1", "TSPAN1", "CEACAM5", "CEACAM6", "IFI27", "DHRS9", "KRT20", "RHOC", "CD177", "PKIB", "HPGD", "LYPD8", "APOBEC1", "APOB", "APOA4", "APOA1", "NPC1L1", "EGFR", "KLF4", "ENPP3", "NT5E", "SLC28A2", "ADA")

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x.copd.epi <- NormalizeData(d10x.copd.epi,scale.factor = 10000, normalization.method = "LogNormalize")

d10x.copd.epi <- AddModuleScore(
  object = d10x.copd.epi,
  features = list(MHCI),
  ctrl = 5,
  name = 'MHCI_score'
)

d10x.copd.epi <- AddModuleScore(
  object = d10x.copd.epi,
  features = list(crypt_villis),
  ctrl = 5,
  name = 'crypt_villis_score'
)



d10x.copd.epi <- NormalizeData(d10x.copd.epi)
d10x.copd.epi <- FindVariableFeatures(d10x.copd.epi, selection.method = "vst", nfeatures = 2000)
d10x.copd.epi <- ScaleData(d10x.copd.epi) 

# dimension reduction
d10x.copd.epi <- RunPCA(d10x.copd.epi, ndims.print = 1:10, nfeatures.print = 10)
d10x.copd.epi <- RunUMAP(d10x.copd.epi, dims = 1:30)



DimPlot(d10x.copd.epi, reduction = "umap", group.by = "Manuscript_Identity", pt.size=0.5, label=T)
ggsave(file=here("../../figs/scRNAseq/primary", "lung_COPD_epicelltype_UMAP.pdf"), w=20,h=12)
ggsave(file=here("../../figs/scRNAseq/primary/jpeg", "lung_COPD_epicelltype_UMAP.jpeg"), w=20,h=12)

DimPlot(d10x.copd.epi, reduction = "umap", group.by = "Disease_Identity", pt.size=0.5)+
  colscale_copd_diagnosis
ggsave(file=here("../../figs/scRNAseq/primary", "lung_COPD_epidisease_UMAP.pdf"), w=15,h=12)
ggsave(file=here("../../figs/scRNAseq/primary/jpeg", "lung_COPD_epidisease_UMAP.jpeg"), w=15,h=12)

DimPlot(d10x.copd.epi, reduction = "umap", group.by = "Subject_Identity", pt.size=0.5)

FeaturePlot(d10x.copd.epi, features = "MHCI_score1",reduction = "umap", min.cutoff = "q9", pt.size=1)
ggsave(file=here("../../figs/scRNAseq/primary", "lung_COPD_epiMHC1_UMAP.pdf"), w=15,h=12)
ggsave(file=here("../../figs/scRNAseq/primary/jpeg", "lung_COPD_epiMHC1_UMAP.jpeg"), w=15,h=12)

FeaturePlot(d10x.copd.epi, features = "MHCI_score1",reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("../../figs/scRNAseq/primary", "lung_COPD_epiMHC1_pca.pdf"), w=15,h=12)
ggsave(file=here("../../figs/scRNAseq/primary/jpeg", "lung_COPD_epiMHC1_pca.jpeg"), w=15,h=12)

FeaturePlot(d10x.copd.epi, features = "crypt_villis_score1",reduction = "umap", min.cutoff = "q9", pt.size=1)
ggsave(file=here("../../figs/scRNAseq/primary", "lung_COPD_epicrypt_villis_UMAP.pdf"), w=15,h=12)
ggsave(file=here("../../figs/scRNAseq/primary/jpeg", "lung_COPD_epicrypt_villis_UMAP.jpeg"), w=15,h=12)


### box plotting
plt<-d10x.copd.epi@meta.data[,c("MHCI_score1","crypt_villis_score1","CellType_Category","Manuscript_Identity","Disease_Identity","Subject_Identity")]
save(plt, file="../../data/public_other/scRNA/lung/GSE136831/scores_each_cell.RData")

load(here("../../data/public_other/scRNA/lung/GSE136831/scores_each_cell.RData"))
plt$Manuscript_Identity<-as.character(plt$Manuscript_Identity)

tapply(plt$MHCI_score1, plt$Disease_Identity, median)

# sig in all epi cells?
print(summary(aov(MHCI_score1 ~ Disease_Identity, data = plt)))
pairwise<-pairwise.t.test(plt$MHCI_score1, plt$Disease_Identity,p.adjust.method = "BH", pool.sd = FALSE)
print(pairwise)
pairwise<-melt(pairwise$p.value)
pairwise<-pairwise[!is.na(pairwise$value),]
pairwise

comp_simple<-list(c("COPD","Control"),c("Control","IPF"),c("COPD","IPF"))

ggplot(plt, aes(Disease_Identity,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=Disease_Identity))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  fillscale_copd_diagnosis+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")#
ggsave(file=here("../../figs/scRNAseq/primary", "lung_COPD_epiMHC1_disease.pdf"), w=7,h=6)
ggsave(file=here("../../figs/scRNAseq/primary/jpeg", "lung_COPD_epiMHC1_disease.jpeg"), w=7,h=6)

##does crypt-villus relate?
cor(plt$MHCI_score1, plt$crypt_villis_score1)


#######
### statistics each cell type
#######
tapply(rownames(plt), list(plt$Disease_Identity,plt$Manuscript_Identity), length)
tapply(plt$MHCI_score1, list(plt$Disease_Identity,plt$Manuscript_Identity), median)

# Group the data by cell type do another anova
stats_plt_disease<-lapply(unique(plt$Manuscript_Identity), function(x){
  print("##################")
  cell_meta<-plt[which(plt$Manuscript_Identity==x),]
  print(unique(cell_meta$Manuscript_Identity))
  print(summary(aov(MHCI_score1 ~ Disease_Identity, data = cell_meta)))
  pairwise<-pairwise.t.test(cell_meta$MHCI_score1, cell_meta$Disease_Identity,p.adjust.method = "BH", pool.sd = FALSE)
  print(pairwise)
  pairwise<-melt(pairwise$p.value)
  pairwise<-pairwise[!is.na(pairwise$value),]
  pairwise$cell_type<-x
  pairwise
})

stats_plt_disease<-do.call(rbind, stats_plt_disease)
stats_plt_disease$Var1<-as.character(stats_plt_disease$Var1)
stats_plt_disease$Var2<-as.character(stats_plt_disease$Var2)


comp_simple<-list(c("COPD","Control"),c("Control","IPF"),c("COPD","IPF"))

### CONFIRM P vale match the pairwise above
ggplot(plt, aes(Disease_Identity,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=Disease_Identity))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~Manuscript_Identity, ncol=5)+fillscale_copd_diagnosis+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")#

ggsave(file=here("../../figs/scRNAseq/primary", "lung_COPD_epiMHC1_boxplot_stimulation.pdf"), w=12,h=8)
ggsave(file=here("../../figs/scRNAseq/primary/jpeg", "lung_COPD_epiMHC1_boxplot_stimulation.jpeg"), w=12,h=8)

stats_plt_disease[which(stats_plt_disease$value<0.05),]


# 
# # Group the data by stimulation do another anova
# stats_plt<-lapply(c("TNFa","NT","IFNg"), function(x){
#   print("##################")
#   cell_meta<-meta[which(meta$orig.ident==x),]
#   print(unique(cell_meta$orig.ident))
#   print(summary(aov(MHCI_score1 ~ cluster_ID, data = cell_meta)))
#   pairwise<-pairwise.t.test(cell_meta$MHCI_score1, cell_meta$cluster_ID,p.adjust.method = "BH", pool.sd = FALSE)
#   print(pairwise)
#   pairwise<-melt(pairwise$p.value)
#   pairwise[!is.na(pairwise$value),]
#   pairwise$orig.ident<-x
#   pairwise
# })
# 
# stats_plt<-do.call(rbind, stats_plt)
# stats_plt$Var1<-as.character(stats_plt$Var1)
# stats_plt$Var2<-as.character(stats_plt$Var2)
# 
# stats_plt<-stats_plt[-which(is.na(stats_plt$value)),]
# 
# comp_simple<-list(c("crypt","enterocyte"),c("crypt","TA"), c("enterocyte","TA"))
# 
# 
# ### CONFIRM P vale match the pairwise above
# ggplot(plt, aes(cluster_ID,MHCI_score1))+
#   geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=cluster_ID))+
#   theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#   facet_wrap(~old.ident)+scale_fill_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))+
#   geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
#               size = 0.3,vjust = 0.5,
#               textsize = 2,  map_signif_level = T, color="grey60")#
# 
# ggsave(file="../../figs/scRNAseq/jpeg/MHCI_score_boxplot_celltype.jpeg", w=12, h=5)
# ggsave(file="../../figs/scRNAseq/MHCI_score_boxplot_celltype.pdf", w=12, h=6)
# 
# 
# 

#############################
#'# individual genes
#############################
IPF_vs_COPD<-read.table(paste(dataset_loc,"manuscript_supplement/aba1983_Data_S10.txt",sep=""))
IPF_vs_ctrl<-read.table(paste(dataset_loc,"manuscript_supplement/aba1983_Data_S8.txt",sep=""))
ctrl_vs_COPD<-read.table(paste(dataset_loc,"manuscript_supplement/aba1983_Data_S9.txt",sep=""))

intersect(IPF_vs_COPD$gene, MHCI)
intersect(IPF_vs_ctrl$gene, MHCI)
intersect(ctrl_vs_COPD$gene, MHCI)

ctrl_vs_COPD[which(ctrl_vs_COPD$gene%in%MHCI),]
IPF_vs_ctrl[which(IPF_vs_ctrl$gene%in%MHCI),]
IPF_vs_COPD[which(IPF_vs_COPD$gene%in%MHCI),]

ctrl_vs_COPD[which(ctrl_vs_COPD$gene%in%MHCI & ctrl_vs_COPD$cellType%in%plt$Manuscript_Identity),]
IPF_vs_ctrl[which(IPF_vs_ctrl$gene%in%MHCI & IPF_vs_ctrl$cellType%in%plt$Manuscript_Identity),]
IPF_vs_COPD[which(IPF_vs_COPD$gene%in%MHCI & IPF_vs_COPD$cellType%in%plt$Manuscript_Identity),]
