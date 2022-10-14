# load meta
load_meta <- function(file){
  meta<-read.table(file, header=T)
  meta$Sample_ID[which(meta$Sample_ID=="C85_caud3pr")]<-"C85_caud5pr"
  meta
}

# load seurat object
load_d10x <- function(file){
  load(here(file))
  d10x.copd.mini
}

#load raw data
load_d10x_raw <- function(dataset_loc){
  #dataset_loc <- here("../../../projects/macparland/RE/PediatricAdult")
  tar_load(meta)
  
  samples<-list.files(dataset_loc)
  samples<-samples[-grep("meta",samples)]
  print(samples)
  
  d10x.list <- sapply(1:length(samples), function(y){
    print(file.path(dataset_loc,paste(samples[y], sep=""),"filtered_feature_bc_matrix"))
    d10x <- Read10X(file.path(dataset_loc,paste(samples[y], sep=""),"filtered_feature_bc_matrix"))
    colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),samples[y],sep="-")
    # print(dim(d10x))
    #' Initialize the Seurat object with the raw (non-normalized data).
    d10x<-CreateSeuratObject(counts = d10x, project = "ped_adult_map", min.cells = 3, min.features = 0)
    d10x
    print(head(meta))
    #add meta data to each seurat object
    meta_cell<-data.frame(cell=colnames(d10x), individual=sapply(colnames(d10x), function(x) strsplit(x,"-")[[1]][2]))
    # print(head(meta_cell))
    meta_cell_add<-merge(meta_cell, meta, by.x="individual", by.y="Sample_ID")
    meta_cell_add<-meta_cell_add[match(colnames(d10x), meta_cell_add$cell),]
    # print(identical(meta_cell_add$cell, colnames(d10x)))
    rownames(meta_cell_add)<-meta_cell_add$cell
    d10x<- AddMetaData(d10x, meta_cell_add)
    d10x
  })
  
  d10x.list
}



## QC
#'The percentage of reads that map to the mitochondrial genome
#'Low-quality / dying cells often exhibit extensive mitochondrial contamination
#'We calculate mitochondrial QC metrics with the PercentageFeatureSet function, which calculates the percentage of counts originating from a set of features
#'We use the set of all genes starting with MT- as a set of mitochondrial genes
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
MT <- function(d10x.list){
  d10x.list.mt<-lapply(1:length(d10x.list), function(x){d10x.list[[x]][["percent.mt"]] <<- PercentageFeatureSet(d10x.list[[x]], pattern = "^MT-")})
  d10x.list.mt
}


QC <- function(d10x.list){
  d10x.list.QC<-lapply(1:length(d10x.list), function(x){d10x.list[[x]] <<- subset(d10x.list[[x]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 50)})
  d10x.list.QC
}



