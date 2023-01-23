meta<-read.table(here("data","data_transfer_updated_jan16_2023.csv"), header=T, sep=",")
meta<-meta[-6,]

d10x_raw_og<-readRDS(file = here("data","d10x_adult_ped_raw.rds"))
table(d10x_raw_og$individual)



count_df<-lapply(1:10, function(x){
  print(meta$Sample_ID[x])
  d10x_raw <- Read10X(file.path(paste("/media/redgar/Seagate Portable Drive/ped_liver_map_raw/",meta$file[x],"/outs",sep=""),"filtered_feature_bc_matrix"))
  colnames(d10x_raw) <- paste(sapply(strsplit(colnames(d10x_raw),split="-"),'[[',1L),"C63_caud5pr",sep="-")
  # print(dim(d10x))
  #' Initialize the Seurat object with the raw (non-normalized data).
  d10x_raw<-CreateSeuratObject(counts = d10x_raw, project = "ped_adult_map", min.cells = 0, min.features = 0)
  
  d10x_raw_og_sample<-subset(d10x_raw_og, subset = individual == meta$Sample_ID[x])
  OG<-sapply(strsplit(colnames(d10x_raw_og_sample),split="-"),'[[',1L)
  
  new<-sapply(strsplit(colnames(d10x_raw),split="-"),'[[',1L)
  data.frame(sample=meta$file[x], cells_OG=length(OG), cells_new=length(new),  overlap=length(intersect(OG,new)))
})


do.call(rbind, count_df)




meta<-read.table(here("data","data_transfer_updated_jan16_2023.csv"), header=T, sep=",")

d10x_raw <- Read10X(file.path(paste("/media/redgar/Seagate Portable Drive/ped_liver_map_raw/",meta$file[6],"/outs",sep=""),"filtered_feature_bc_matrix"))
colnames(d10x_raw) <- paste(sapply(strsplit(colnames(d10x_raw),split="-"),'[[',1L),"C85_caud3pr",sep="-")
# print(dim(d10x))
#' Initialize the Seurat object with the raw (non-normalized data).
d10x_raw<-CreateSeuratObject(counts = d10x_raw, project = "ped_adult_map", min.cells = 0, min.features = 0)
d10x_raw_og_85<-subset(d10x_raw_og, subset = individual == "C85_caud5pr")


new<-sapply(strsplit(colnames(d10x_raw),split="-"),'[[',1L)
OG85<-sapply(strsplit(colnames(d10x_raw_og_85),split="-"),'[[',1L)

data.frame(sample=meta$file[x], cells_new=length(new),  overlap=length(intersect(OG85,new)))




