library(here)
library(Seurat)


dataset_loc <- here("../../../projects/macparland/RE/PediatricAdult")

samples<-list.files(dataset_loc)
samples<-samples[-grep("meta",samples)]
print(samples)

meta<-read.table(here(dataset_loc,"input_metadata.txt"), header=T)
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
  # print(head(meta_cell_add))
  
  d10x<- AddMetaData(d10x, meta_cell_add)
  d10x
  })

d10x.list

## cell counts
plt_count_raw<-do.call(rbind,lapply(1:length(d10x.list_primary), function(x) {
  df<-as.data.frame(tapply(rownames(d10x.list_primary[[x]]@meta.data), list(d10x.list_primary[[x]]@meta.data$individual, d10x.list_primary[[x]]@meta.data$orig.ident), length))
  df$individual<-rownames(df)
  df$condition<-names(d10x.list_primary)[x]
  colnames(df)<-c("cell_count","individual","condition")
  df}))
plt_count_raw

ggplot(plt_count_raw, aes(condition, cell_count))+
  geom_boxplot(fill="lightgrey")+geom_point()+
  theme_bw()+geom_text(aes(label=individual), hjust=-0.25)+ylim(0,15000)+ylab("Total Cell Number")



#'## QC
#'The percentage of reads that map to the mitochondrial genome
#'Low-quality / dying cells often exhibit extensive mitochondrial contamination
#'We calculate mitochondrial QC metrics with the PercentageFeatureSet function, which calculates the percentage of counts originating from a set of features
#'We use the set of all genes starting with MT- as a set of mitochondrial genes

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
invisible(lapply(1:length(d10x.list_primary), function(x){
  d10x.list_primary[[x]][["percent.mt"]] <<- PercentageFeatureSet(d10x.list_primary[[x]], pattern = "^MT-")}))

# Show QC metrics for the first 5 cells
head(d10x.list_primary[[2]]@meta.data, 5)

#'Low-quality cells or empty droplets will often have very few genes
#'Cell doublets or multiplets may exhibit an aberrantly high gene count


# Visualize QC metrics
#nFeature number of unique genes
#nCount number of total molecules
plt_QC_data<-do.call(rbind, lapply(1:4, function(x) d10x.list_primary[[x]]@meta.data))

ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th

ggsave(file=here("figs","primary_ncount_nfeatures.pdf"), w=6,h=4)
ggsave(file=here("figs/jpeg","primary_ncount_nfeatures.jpeg"), w=6,h=4)


ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
  geom_vline(xintercept = 10)+ theme_bw()+xlab("Percent Mitochondrial")+th

ggsave(file=here("figs","primary_MT_features.pdf"), w=4,h=2)
ggsave(file=here("figs/jpeg","primary_MT_features.jpeg"), w=4,h=2)




#' #### scrublet results
samples<-meta_primary$Sanger.Sample.ID

doublet_scores_all<-lapply(1:length(samples), function(x){
  ## scrublet output
  doublet_scores<-read.csv(paste(here("data/"),"scrublet_output_table_",samples[x],".csv",sep=""))
  doublet_scores$index_edit<-gsub("1",samples[x],doublet_scores$index)
  doublet_scores
})
doublet_scores_all<-do.call(rbind, doublet_scores_all)

## Match doublet information
doublet_scores_list<-lapply(1:length(d10x.list_primary), function(x){
  doublet_scores_all<-doublet_scores_all[match(colnames(d10x.list_primary[[x]]), doublet_scores_all$index_edit),]
  identical(colnames(d10x.list_primary[[x]]), doublet_scores_all$index_edit)
  doublet_scores_all
})

invisible(lapply(1:length(d10x.list_primary), function(x){
  d10x.list_primary[[x]]<<- AddMetaData(d10x.list_primary[[x]], doublet_scores_list[[x]]$doublet_score, col.name = "doublet_score")
  d10x.list_primary[[x]]<<- AddMetaData(d10x.list_primary[[x]], doublet_scores_list[[x]]$predicted_doublet, col.name = "predicted_doublet")
}))

plt_QC_data<-do.call(rbind, lapply(1:4, function(x) d10x.list_primary[[x]]@meta.data))

ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=doublet_score)) + 
  geom_point(size=0.75) + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Doublet\nScore") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th

ggsave(file=here("figs","primary_doublet_features.pdf"),w=5,h=3)
ggsave(file=here("figs/jpeg","primary_doublet_features.jpeg"), w=5,h=3)

ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=doublet_score)) + 
  geom_point(size=0.75) + facet_wrap(~predicted_doublet)+
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Doublet\nScore") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th

ggsave(file=here("figs","primary_doublet_features_split.pdf"), w=12,h=4)
ggsave(file=here("figs/jpeg","primary_doublet_features_spit.jpeg"), w=12,h=4)

table(plt_QC_data$predicted_doublet)
dim(plt_QC_data)
table(plt_QC_data$predicted_doublet, plt_QC_data$individual)



## what is 3 MAD from mean
median(plt_QC_data$nFeature_RNA)+(mad(plt_QC_data$nFeature_RNA)*3)
median(plt_QC_data$nFeature_RNA)-(mad(plt_QC_data$nFeature_RNA)*3)


#'We filter cells that have unique feature counts over 6,000 or less than 500
#'We filter cells that have >10% mitochondrial counts
#'we will also filter doublets as called by scrublet
d10x.list_primary.raw<-d10x.list_primary

invisible(lapply(1:length(d10x.list_primary), function(x){
  d10x.list_primary[[x]] <<- subset(d10x.list_primary[[x]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10 & predicted_doublet=="False")
}))

d10x.list_primary

## cell counts after QC
plt_count_QC<-do.call(rbind,lapply(1:length(d10x.list_primary), function(x) {
  df<-as.data.frame(tapply(rownames(d10x.list_primary[[x]]@meta.data), list(d10x.list_primary[[x]]@meta.data$individual, d10x.list_primary[[x]]@meta.data$orig.ident), length))
  df$individual<-rownames(df)
  df$condition<-names(d10x.list_primary)[x]
  colnames(df)<-c("cell_count","individual","condition")
  df}))
plt_count_QC

#save to plot fancy
save(plt_count_raw, plt_count_QC, file=here("data","primary_cell_count.RData"))

plt_count_raw<-plt_count_raw[which(plt_count_raw$condition!="neo"),]
plt_count_raw$condition<-as.factor(plt_count_raw$condition)
levels(plt_count_raw$condition)<-c("CD","Control","UC")
plt_count_raw$condition<-factor(plt_count_raw$condition, levels=c("Control","CD","UC"))

plt_count_QC<-plt_count_QC[which(plt_count_QC$condition!="neo"),]
plt_count_QC$condition<-as.factor(plt_count_QC$condition)
levels(plt_count_QC$condition)<-c("CD","Control","UC")
plt_count_QC$condition<-factor(plt_count_QC$condition, levels=c("Control","CD","UC"))


## only samples in thesis
plt_count_raw<-plt_count_raw[-grep("036NEG|036POS",plt_count_raw$individual),]
plt_count_QC<-plt_count_QC[-grep("036NEG|036POS",plt_count_QC$individual),]

cell_count<-grid.arrange(ggplot(plt_count_raw, aes(condition, cell_count,fill=condition))+
                           geom_boxplot()+geom_point()+
                           theme_bw()+geom_text(aes(label=individual), hjust=-0.25, size=3)+ylim(0, 15000)+xlab("Diagnosis")+
                           ylab("Total Cell Number")+th+fillscale_diagnosis+
                           theme(legend.position = "none")+ggtitle("Before Quality Control"),
                         ggplot(plt_count_QC, aes(condition, cell_count,fill=condition))+
                           geom_boxplot()+geom_point()+
                           theme_bw()+geom_text(aes(label=individual), hjust=-0.25, size=3)+ylim(0, 15000)+xlab("Diagnosis")+
                           ylab("Total Cell Number")+th+fillscale_diagnosis+
                           theme(legend.position = "none")+ggtitle("After Quality Control"), ncol=2)

ggsave(cell_count,file=here("figs","primary_scCount.pdf"), w=8,h=4)
ggsave(cell_count,file=here("figs/jpeg","primary_scCount.jpeg"), w=8,h=4)


# $CD
# An object of class Seurat 
# 23944 features across 24996 samples within 1 assay 
# Active assay: RNA (23944 features, 0 variable features)
# 
# $Ctrl
# An object of class Seurat 
# 23904 features across 27622 samples within 1 assay 
# Active assay: RNA (23904 features, 0 variable features)
# 
# $UC
# An object of class Seurat 
# 19768 features across 2703 samples within 1 assay 
# Active assay: RNA (19768 features, 0 variable features)
# 
# $neo
# An object of class Seurat 
# 20196 features across 3743 samples within 1 assay 
# Active assay: RNA (20196 features, 0 variable features)


d10x.primary <- merge(d10x.list_primary[[1]], d10x.list_primary[[2]])
d10x.primary <- merge(d10x.primary, d10x.list_primary[[3]])
d10x.primary <- merge(d10x.primary, d10x.list_primary[[4]])

# An object of class Seurat 
# 24939 features across 59320 samples within 1 assay 
# Active assay: RNA (24939 features, 0 variable features)

d10x.primary.raw<-d10x.primary

saveRDS(d10x.primary.raw, file = here("data","d10x_primary_raw_merged.rds"))




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

