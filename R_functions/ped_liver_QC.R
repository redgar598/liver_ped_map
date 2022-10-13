### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)

source("R_functions/pretty_plots.R")


dataset_loc <- here("../../../projects/macparland/RE/PediatricAdult")

samples<-list.files(dataset_loc)
samples<-samples[-grep("meta",samples)]
print(samples)

meta<-read.table(here(dataset_loc,"input_metadata.txt"), header=T)
meta$Sample_ID[which(meta$Sample_ID=="C85_caud3pr")]<-"C85_caud5pr"
#meta<-read.table(here("data/input_metadata.txt"), header=T)


d10x.list <- sapply(1:length(samples), function(y){
  print(file.path(dataset_loc,paste(samples[y], sep=""),"filtered_feature_bc_matrix"))
  d10x <- Read10X(file.path(dataset_loc,paste(samples[y], sep=""),"filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),samples[y],sep="-")
  # print(dim(d10x))
  #' Initialize the Seurat object with the raw (non-normalized data).
  d10x<-CreateSeuratObject(counts = d10x, project = "ped_adult_map", min.cells = 3, min.features = 0)
  
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



## cell counts
plt_count_raw<-lapply(1:length(d10x.list), function(x) {
  df<-data.frame(raw_cell_count=nrow(d10x.list[[x]]@meta.data),individual=unique(d10x.list[[x]]@meta.data$individual))
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
  d10x.list[[x]][["percent.mt"]] <<- PercentageFeatureSet(d10x.list[[x]], pattern = "^MT-")}))

# Show QC metrics for the first 5 cells
print(head(d10x.list[[2]]@meta.data, 5))

#'Low-quality cells or empty droplets will often have very few genes
#'Cell doublets or multiplets may exhibit an aberrantly high gene count


# Visualize QC metrics
#nFeature number of unique genes
#nCount number of total molecules
plt_QC_data<-do.call(rbind, lapply(1:4, function(x) d10x.list[[x]]@meta.data))

qc_plts<-ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th
save_plts(qc_plts, "intital_QC_plts", w=6,h=4)


MT_plt<-ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
  geom_vline(xintercept = 10)+ theme_bw()+xlab("Percent Mitochondrial")+th
save_plts(qc_plts, "percentMT_plt", w=6,h=4)





#'We filter cells that have unique feature counts over 6,000 or less than 500
#'We filter cells that have >10% mitochondrial counts
#'we will also filter doublets as called by scrublet
d10x.list.raw<-d10x.list

invisible(lapply(1:length(d10x.list), function(x){
  d10x.list[[x]] <<- subset(d10x.list[[x]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
}))

d10x.list

## cell counts after QC
plt_count_QC<-lapply(1:length(d10x.list), function(x) {
  df<-data.frame(qc_cell_count=nrow(d10x.list[[x]]@meta.data),individual=unique(d10x.list[[x]]@meta.data$individual))
  df})
plt_count_QC<-do.call(rbind, plt_count_QC)
print(plt_count_QC)

counts<-merge(plt_count_raw, plt_count_QC, by="individual")
meta<-merge(meta,counts,by.x="Sample_ID", by.y="individual")

cell_count<-grid.arrange(ggplot(meta, aes(AgeGroup, raw_cell_count,fill=AgeGroup))+
                           geom_boxplot()+geom_point()+
                           theme_bw()+geom_text(aes(label=individual), hjust=-0.25, size=3)+ylim(0, 15000)+xlab("Age Group")+
                           ylab("Total Cell Number")+th+fillscale_age+
                           theme(legend.position = "none")+ggtitle("Before Quality Control"),
                         ggplot(meta, aes(AgeGroup, qc_cell_count,fill=AgeGroup))+
                           geom_boxplot()+geom_point()+
                           theme_bw()+geom_text(aes(label=individual), hjust=-0.25, size=3)+ylim(0, 15000)+xlab("Age Group")+
                           ylab("Total Cell Number")+th+fillscale_age+
                           theme(legend.position = "none")+ggtitle("After Quality Control"), ncol=2)

save_plts(cell_count, "QC_cellcount_age", w=8,h=4)


# d10x.primary <- merge(d10x.list_primary[[1]], d10x.list_primary[[2]])
# d10x.primary <- merge(d10x.primary, d10x.list_primary[[3]])
# d10x.primary <- merge(d10x.primary, d10x.list_primary[[4]])
# 
# # An object of class Seurat 
# # 24939 features across 59320 samples within 1 assay 
# # Active assay: RNA (24939 features, 0 variable features)
# 
# d10x.primary.raw<-d10x.primary
# 
# saveRDS(d10x.primary.raw, file = here("data","d10x_primary_raw_merged.rds"))
# 




print(sessionInfo())


