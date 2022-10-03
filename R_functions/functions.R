# load seurat object
load_d10x <- function(file){
  load(here(file))
  d10x.copd.mini
}


## QC
#'The percentage of reads that map to the mitochondrial genome
#'Low-quality / dying cells often exhibit extensive mitochondrial contamination
#'We calculate mitochondrial QC metrics with the PercentageFeatureSet function, which calculates the percentage of counts originating from a set of features
#'We use the set of all genes starting with MT- as a set of mitochondrial genes
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
MT <- function(d10x){
  d10x[["percent.mt"]] <- PercentageFeatureSet(d10x, pattern = "^MT-")
  d10x
}


QC <- function(d10x){
  d10x_QC <- subset(d10x, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)# & predicted_doublet=="False")
  d10x_QC
}
