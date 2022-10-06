# load seurat object
load_d10x <- function(file){
  load(here(file))
  d10x.copd.mini
}

#load raw data
load_d10x_raw <- function(dataset_loc){
  #dataset_loc <- here("../../../projects/macparland/RE/PediatricAdult")
  
  samples<-list.files(dataset_loc)
  
  d10x.data <- sapply(samples, function(y){
    d10x <- Read10X(file.path(dataset_loc,paste(samples[y], sep=""),"filtered_feature_bc_matrix"))
    colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
    d10x
  })
  d10x.data<-do.call("cbind",d10x.data)
    
  #' Initialize the Seurat object with the raw (non-normalized data).
  d10x <- CreateSeuratObject(counts = d10x.data, project = "ped_adult_map", min.cells = 3, min.features = 0)
  d10x
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
