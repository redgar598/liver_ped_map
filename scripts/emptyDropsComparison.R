# ---
#   title: "Pipeline with Empty Drops"
# output: html_notebook
# ---
  
#This script will filter my healthy woodchuck data with Empty Drops and compare it to the dataset filtered with CellRanger and DropletQC.


library(here)
library(DropletUtils)
library(Seurat)
library(dplyr)
# library(scClustViz)
# library(viridis)
# library(presto)
# library(cluster)
# library(viridisLite)
# library(shiny)


# Set sample variables -->
  
  
# woodchuck <- "L212"
# sample <- "TLH"
lower <- 800


# Get woodchuck metadata -->
  
  #  -->
  # source("~/Dropbox/Zoe/scf_version/analysis/scripts/woodchuckMetadata.R") -->
  #  -->
  
  # Extract sample info -->
  
  #  -->
  # sampleInfo <- woodchuck_info[woodchuck_info$woodchuck == woodchuck &  -->
                                      #                                woodchuck_info$sample_type == sample,] -->
  #  -->
  
  #Read in raw data


dataset_loc <- here("/media/redgar/Seagate Portable Drive/ped_liver_map_raw")

samples<-list.files(dataset_loc)
print(samples)

meta<-read.table(here("data/data_transfer_updated_may15_2023_IFALD_PBMC.csv"), header=T, sep=",")
samples<-samples[grep("C97",samples)]


sobj.data <- Read10X(file.path(dataset_loc,paste(samples,"/outs", sep=""),"raw_feature_bc_matrix"))




#From these counts, make single cell experiment object


sobj <- CreateSeuratObject(counts = sobj.data)
sce <- as.SingleCellExperiment(sobj)
sce


#Plot counts


br.out <- barcodeRanks(counts(sce))
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
       legend=c("knee", "inflection"))


#Calculate empty drops with `lower = 800` since default `lower` returns tens of thousands of cells


set.seed(100)
e.out <- emptyDrops(counts(sce), lower=lower)
e.out


#Detect significant deviations from the ambient profile by setting FDR threshold


# Keep cells that have an FDR of less than 1%
is.cell <- e.out$FDR <= 0.01
# Count the cells that are remaining
sum(is.cell, na.rm=TRUE)


#Look for entries that have `FDR` above the desired threshold and `Limited==TRUE`, it indicates that `npts` should be increased in the `emptyDrops` function


table(Limited=e.out$Limited, Significant=is.cell)


#Visualise the results. Droplets detected as cells should have large negative log-probabilities OR very large total counts


plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
     xlab="Total UMI count", ylab="-Log Probability")


#Subset object to keep only cells


namedCounts <- counts(sce)
colnames(namedCounts) <- rownames(sce@colData)
cell.counts <- namedCounts[,which(is.cell),drop=FALSE]
dim(cell.counts)


#Convert to Seurat object, don't use `min.features` since `emptyDrops` was run.


sobj <- CreateSeuratObject(counts = cell.counts, min.cells = 3)
sobj


#Add metadata
meta_cell<-data.frame(cell=colnames(sobj), individual="C97_caud3pr")
meta_cell_add<-merge(meta_cell, meta, by.x="individual", by.y="Sample_ID")
meta_cell_add<-meta_cell_add[match(colnames(sobj), meta_cell_add$cell),]
print(identical(meta_cell_add$cell, colnames(sobj)))
rownames(meta_cell_add)<-meta_cell_add$cell
sobj<- AddMetaData(sobj, meta_cell_add)

#Store percentage of mitocondrial and ribosomal reads, and visualise


sobj <- PercentageFeatureSet(sobj, pattern = "^MT-", col.name = "percent.mt")
sobj <- PercentageFeatureSet(sobj, pattern = "^RPL|^Rpl|^RPS|^Rps", col.name = "percent.ribo")
# Visualize
VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE)


#Apply mitochondrial filter


mitoCells <- WhichCells(sobj, expression = percent.mt < 50) 
# Only keep these cells
sobj <- subset(sobj, cells = mitoCells)
sobj

save(sobj, file=here("/media/redgar/Seagate Portable Drive/processed_data/empty_drop_C97_conservative.RData"))

load(here("/media/redgar/Seagate Portable Drive/processed_data/empty_drop_C97_conservative.RData"))

##Run SCTransform

# 
# sobj <- SCTransform(sobj) %>%
#   RunPCA() %>%
#   FindNeighbors(dims = 1:30) %>% 
#   RunUMAP(dims = 1:30) %>%
#   RunTSNE()

sobj <- NormalizeData(sobj)
sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
sobj <- ScaleData(sobj) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))

# dimension reduction
sobj <- RunPCA(sobj, ndims.print = 1:10, nfeatures.print = 10)
sobj <- RunUMAP(sobj, dims = 1:30)
sobj <- RunTSNE(sobj, dims = 1:30)

sobj <- FindNeighbors(sobj, reduction = "pca", dims = 1:20)
sobj <- FindClusters(sobj, resolution = 0.2)



            # #Set up for scClustViz
            # 
            # 
            # scSeurat <- sobj
            # DE_bw_clust <- TRUE
            # FDRthresh <- 0.01
            # seurat_resolution <- 0
            # sCVdata_list <- list()
            # 
            # 
            # #Cluster
            # 
            # 
            # while(DE_bw_clust) {
            #   seurat_resolution <- seurat_resolution + 0.2
            #   # ^ iteratively incrementing resolution parameter
            #   
            #   scSeurat <- FindClusters(scSeurat,
            #                            resolution=seurat_resolution,
            #                            print.output=F)
            #   message(" ")
            #   message("------------------------------------------------------")
            #   message(paste0("--------  res.",seurat_resolution," with ",
            #                  length(levels(Idents(scSeurat)))," clusters --------"))
            #   message("------------------------------------------------------")
            #   curr_sCVdata <- CalcSCV(inD=scSeurat,
            #                           cl=Idents(scSeurat),
            #                           # ^ your most recent clustering results get stored in the Seurat "ident" slot
            #                           assayType = "SCT",
            #                           exponent=exp(1),
            #                           pseudocount=1,
            #                           DRthresh=0.1,
            #                           DRforClust="pca",
            #                           calcSil=T,
            #                           calcDEvsRest=T,
            #                           calcDEcombn=T)
            #   
            #   DE_bw_NN <- sapply(DEneighb(curr_sCVdata,FDRthresh),nrow)
            #   # ^ counts # of DE genes between neighbouring clusters at your selected FDR threshold
            #   
            #   if (min(DE_bw_NN) < 1 | seurat_resolution > 2) { DE_bw_clust <- FALSE }
            #   # ^ If no DE genes between nearest neighbours, don't loop again.
            # 
            # sCVdata_list[[paste0("res.",seurat_resolution)]] <- curr_sCVdata
            # }
            # 
            # 
            # Save object
            # 
            # 
            # save(sCVdata_list, scSeurat,
            #      file = paste(woodchuck, sample, "emptyDrops_lower", lower, "_scClustViz.RData", sep = "_"))
            # 
            # 
            # OR if the object is too large, run Seurat clustering loop
            # 

# for (res in seq(from = 0.2, to = 2.0, by = 0.2)) {
#   sobj <- FindClusters(sobj, resolution = res)
# }
# Idents(sobj) <- "SCT_snn_res.0.2"
# 
# 
# Save object


#saveRDS(sobj, file = paste(woodchuck, sample, "emptyDrops_lower", lower, "_sobj.RDS", sep = "_"))


### Visualisations

#Reassign object to new variable and load in `DropletQC` filtered object

#Add DropletQC result to `sobj`

load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label_C97<-cell_label[which(cell_label$individual=="C97_caud3pr"),]
cell_label_C97$cell<-gsub("_14", "",cell_label_C97$cell)

cell_label_C97<-cell_label_C97[which(cell_label_C97$cell%in%colnames(sobj)),]
meta_cell_add<-merge(sobj@meta.data[c("nCount_RNA","nFeature_RNA" , "individual","cell","percent.mt", "percent.ribo")], cell_label_C97[,c("cell","CellType_refined")], all.x=T, by="cell")

meta_cell_add<-meta_cell_add[match(colnames(sobj), meta_cell_add$cell),]
print(identical(meta_cell_add$cell, colnames(sobj)))
rownames(meta_cell_add)<-meta_cell_add$cell
sobj<- AddMetaData(sobj, meta_cell_add)


load("/home/redgar/Documents/liver_ped_map/data/QC_metrics.Rdata")
plt_QC_C97<-plt_QC_data[which(plt_QC_data$individual=="C97_caud3pr"),]

plt_QC_C97<-plt_QC_C97[which(plt_QC_C97$cell%in%colnames(sobj)),]
meta_cell_add<-merge(sobj@meta.data[c("nCount_RNA","nFeature_RNA" , "individual","cell","percent.mt", "percent.ribo")], plt_QC_C97, all.x=T, by="cell")

meta_cell_add<-meta_cell_add[match(colnames(sobj), meta_cell_add$cell),]
print(identical(meta_cell_add$cell, colnames(sobj)))
rownames(meta_cell_add)<-meta_cell_add$cell
sobj<- AddMetaData(sobj, meta_cell_add)


table(sobj$cell_status)


#With `lower = 800`, datasets look very similar. Try increasing lower?
  
  
status<-DimPlot(sobj, group.by = "cell_status")
save_plts(status, "cell_status_human_C97", w=8, h=6)

fraction<-FeaturePlot(sobj, reduction="umap", features = "nuclear_fraction")
save_plts(fraction, "nuclear_fraction_human_C97", w=8, h=6)

celltype<-DimPlot(sobj, group.by = "CellType_refined", label=T)+colscale_cellType
save_plts(celltype, "celltype_human_C97", w=10, h=6)


save(sobj, file=here("/media/redgar/Seagate Portable Drive/processed_data/empty_drop_C97_default.RData"))

load(here("/media/redgar/Seagate Portable Drive/processed_data/empty_drop_C97_default.RData"))
sobj_default<-sobj
load(here("/media/redgar/Seagate Portable Drive/processed_data/empty_drop_C97_conservative.RData"))
sobj_800<-sobj
rm(sobj)
gc()




sobj_default$cut_at_lower800<-"0"
sobj_default$cut_at_lower800[which(!(sobj_default$cell%in%sobj_800$cell))]<-"1"
table(sobj_default$cut_at_lower800)


umap_mat_myeloid<-as.data.frame(Embeddings(object = sobj_default, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-sobj_default@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

cut<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
    geom_point(aes(fill=CellType_refined, color=cut_at_lower800),size=1, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
    fillscale_cellType+theme_classic()+scale_color_manual(values=c("white","black"))

save_plts(cut, "cut_at_lower800_human_C97", w=10, h=6)





##################
## Multiple level comparison
##################

library(here)
library(DropletUtils)
library(Seurat)
library(dplyr)


#Read in raw data
dataset_loc <- here("/media/redgar/Seagate Portable Drive/ped_liver_map_raw")

samples<-list.files(dataset_loc)
print(samples)

meta<-read.table(here("data/data_transfer_updated_may15_2023_IFALD_PBMC.csv"), header=T, sep=",")
samples<-samples[grep("C97",samples)]

sobj.data <- Read10X(file.path(dataset_loc,paste(samples,"/outs", sep=""),"raw_feature_bc_matrix"))


#From these counts, make single cell experiment object
sobj <- CreateSeuratObject(counts = sobj.data)
sce <- as.SingleCellExperiment(sobj)
sce

rm(sobj.data)
gc()


cell_lower<-lapply(c(400, 800, 1000, 1200, 1400), function(lower){
  #Calculate empty drops with `lower = 800` since default `lower` returns tens of thousands of cells
  set.seed(100)
  print(lower)
  e.out <- emptyDrops(counts(sce), lower=lower)
  
  #Detect significant deviations from the ambient profile by setting FDR threshold
  
  # Keep cells that have an FDR of less than 1%
  is.cell <- e.out$FDR <= 0.01
  # Count the cells that are remaining
  print(sum(is.cell, na.rm=TRUE))
  
  #Subset object to keep only cells
  namedCounts <- counts(sce)
  colnames(namedCounts) <- rownames(sce@colData)
  cell.counts <- namedCounts[,which(is.cell),drop=FALSE]
  gc()
  
  #Convert to Seurat object, don't use `min.features` since `emptyDrops` was run.
  sobj <- CreateSeuratObject(counts = cell.counts, min.cells = 3)

  #Add metadata
  meta_cell<-data.frame(cell=colnames(sobj), individual="C97_caud3pr")
  meta_cell_add<-merge(meta_cell, meta, by.x="individual", by.y="Sample_ID")
  meta_cell_add<-meta_cell_add[match(colnames(sobj), meta_cell_add$cell),]
  print(identical(meta_cell_add$cell, colnames(sobj)))
  rownames(meta_cell_add)<-meta_cell_add$cell
  sobj<- AddMetaData(sobj, meta_cell_add)
  
  #Store percentage of mitocondrial and ribosomal reads, and visualise
  sobj <- PercentageFeatureSet(sobj, pattern = "^MT-", col.name = "percent.mt")
  sobj <- PercentageFeatureSet(sobj, pattern = "^RPL|^Rpl|^RPS|^Rps", col.name = "percent.ribo")
  
  #Apply mitochondrial filter
  mitoCells <- WhichCells(sobj, expression = percent.mt < 50) 
  # Only keep these cells
  sobj <- subset(sobj, cells = mitoCells)
  sobj
  gc()
  as.character(sobj$cell)})


load("/home/redgar/Documents/liver_ped_map/data/QC_metrics.Rdata")
plt_QC_C97<-plt_QC_data[which(plt_QC_data$individual=="C97_caud3pr"),]

lowers<-c(400, 800, 1000, 1200, 1400)
lapply(1:length(lowers), function(x){
  lowers[x]
  plt_QC_C97$lower<<-"Cut by emptyDrops"
  plt_QC_C97$lower[which(plt_QC_C97$cell%in%cell_lower[[x]])]<<-"Kept by emptyDrops"
  colnames(plt_QC_C97)[which(colnames(plt_QC_C97)=="lower")]<<-as.character(paste("lower",lowers[x], sep="_"))
})





library(ggsankey)

df <- plt_QC_C97 %>%
  make_long(cell_status, lower_400)

dagg <- df%>%
  dplyr::group_by(node)%>%
  tally()

df <- merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)


cut<-ggplot(df, aes(x = x,
                    next_x = next_x,
                    node = node,
                    next_node = next_node,
                    fill = factor(node),
                    label = n)) +
  geom_sankey() +xlab("")+ggtitle("Lower 400")+
  theme_sankey(base_size = 16)+
  geom_sankey(flow.alpha = 0.25,  color = "gray40", show.legend = TRUE)+
  geom_sankey_label(size = 3, color = "white", fill= "gray40", hjust = -0.2)+
  scale_fill_manual(values = c("cornflowerblue","grey40","red","grey","dodgerblue4"))
save_plts(cut, "sankey_lower400_human_C97", w=10, h=6)

df <- plt_QC_C97 %>%
  make_long(cell_status, lower_800)
dagg <- df%>%
  dplyr::group_by(node)%>%
  tally()

df <- merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)

cut<-ggplot(df, aes(x = x,
               next_x = next_x,
               node = node,
               next_node = next_node,
               fill = factor(node),
               label = n)) +
  geom_sankey() +xlab("")+ggtitle("Lower 800")+
  theme_sankey(base_size = 16)+
  geom_sankey(flow.alpha = 0.25,  color = "gray40", show.legend = TRUE)+
  geom_sankey_label(size = 3, color = "white", fill= "gray40", hjust = -0.2)+
  scale_fill_manual(values = c("cornflowerblue","grey40","red","grey","dodgerblue4"))
save_plts(cut, "sankey_lower800_human_C97", w=10, h=6)


df <- plt_QC_C97 %>%
  make_long(cell_status, lower_1000)
dagg <- df%>%
  dplyr::group_by(node)%>%
  tally()

df <- merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)

cut<-ggplot(df, aes(x = x,
                    next_x = next_x,
                    node = node,
                    next_node = next_node,
                    fill = factor(node),
                    label = n)) +
  geom_sankey() +xlab("")+ggtitle("Lower 1000")+
  theme_sankey(base_size = 16)+
  geom_sankey(flow.alpha = 0.25,  color = "gray40", show.legend = TRUE)+
  geom_sankey_label(size = 3, color = "white", fill= "gray40", hjust = -0.2)+
  scale_fill_manual(values = c("cornflowerblue","grey40","red","grey","dodgerblue4"))
save_plts(cut, "sankey_lower1000_human_C97", w=10, h=6)


df <- plt_QC_C97 %>%
  make_long(cell_status, lower_1200)
dagg <- df%>%
  dplyr::group_by(node)%>%
  tally()

df <- merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)

cut<-ggplot(df, aes(x = x,
                    next_x = next_x,
                    node = node,
                    next_node = next_node,
                    fill = factor(node),
                    label = n)) +
  geom_sankey() +xlab("")+ggtitle("Lower 1200")+
  theme_sankey(base_size = 16)+
  geom_sankey(flow.alpha = 0.25,  color = "gray40", show.legend = TRUE)+
  geom_sankey_label(size = 3, color = "white", fill= "gray40", hjust = -0.2)+
  scale_fill_manual(values = c("cornflowerblue","grey40","red","grey","dodgerblue4"))
save_plts(cut, "sankey_lower1200_human_C97", w=10, h=6)

df <- plt_QC_C97 %>%
  make_long(cell_status, lower_1400)
dagg <- df%>%
  dplyr::group_by(node)%>%
  tally()

df <- merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)

cut<-ggplot(df, aes(x = x,
                    next_x = next_x,
                    node = node,
                    next_node = next_node,
                    fill = factor(node),
                    label = n)) +
  geom_sankey() +xlab("")+ggtitle("Lower 1400")+
  theme_sankey(base_size = 16)+
  geom_sankey(flow.alpha = 0.25,  color = "gray40", show.legend = TRUE)+
  geom_sankey_label(size = 3, color = "white", fill= "gray40", hjust = -0.2)+
  scale_fill_manual(values = c("cornflowerblue","grey40","red","grey","dodgerblue4"))
save_plts(cut, "sankey_lower1400_human_C97", w=10, h=6)



########################
## C68
########################
# ---
#   title: "Pipeline with Empty Drops"
# output: html_notebook
# ---

#This script will filter my healthy woodchuck data with Empty Drops and compare it to the dataset filtered with CellRanger and DropletQC.


library(here)
library(DropletUtils)
library(Seurat)
library(dplyr)
# library(scClustViz)
# library(viridis)
# library(presto)
# library(cluster)
# library(viridisLite)
# library(shiny)


# Set sample variables -->


# woodchuck <- "L212"
# sample <- "TLH"
lower <- 800


# Get woodchuck metadata -->

#  -->
# source("~/Dropbox/Zoe/scf_version/analysis/scripts/woodchuckMetadata.R") -->
#  -->

# Extract sample info -->

#  -->
# sampleInfo <- woodchuck_info[woodchuck_info$woodchuck == woodchuck &  -->
#                                woodchuck_info$sample_type == sample,] -->
#  -->

#Read in raw data


dataset_loc <- here("/media/redgar/Seagate Portable Drive/ped_liver_map_raw")

samples<-list.files(dataset_loc)
print(samples)

meta<-read.table(here("data/data_transfer_updated_may15_2023_IFALD_PBMC.csv"), header=T, sep=",")
samples<-samples[grep("C61",samples)]


sobj.data <- Read10X(file.path(dataset_loc,paste(samples,"/outs", sep=""),"raw_feature_bc_matrix"))




#From these counts, make single cell experiment object


sobj <- CreateSeuratObject(counts = sobj.data)
sce <- as.SingleCellExperiment(sobj)
sce
# 
# 
# #Plot counts
# 
# 
# br.out <- barcodeRanks(counts(sce))
# plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
# o <- order(br.out$rank)
# lines(br.out$rank[o], br.out$fitted[o], col="red")
# abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
# abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
# legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
#        legend=c("knee", "inflection"))
# 
# 
# #Calculate empty drops with `lower = 800` since default `lower` returns tens of thousands of cells
# 
# 
# set.seed(100)
# e.out <- emptyDrops(counts(sce), lower=lower)
# e.out
# 
# 
# #Detect significant deviations from the ambient profile by setting FDR threshold
# 
# 
# # Keep cells that have an FDR of less than 1%
# is.cell <- e.out$FDR <= 0.01
# # Count the cells that are remaining
# sum(is.cell, na.rm=TRUE)
# 
# 
# #Look for entries that have `FDR` above the desired threshold and `Limited==TRUE`, it indicates that `npts` should be increased in the `emptyDrops` function
# 
# 
# table(Limited=e.out$Limited, Significant=is.cell)
# 
# 
# #Visualise the results. Droplets detected as cells should have large negative log-probabilities OR very large total counts
# 
# 
# plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
#      xlab="Total UMI count", ylab="-Log Probability")
# 
# 
# #Subset object to keep only cells
# 
# 
# namedCounts <- counts(sce)
# colnames(namedCounts) <- rownames(sce@colData)
# cell.counts <- namedCounts[,which(is.cell),drop=FALSE]
# dim(cell.counts)
# 
# 
# #Convert to Seurat object, don't use `min.features` since `emptyDrops` was run.
# 
# 
# sobj <- CreateSeuratObject(counts = cell.counts, min.cells = 3)
# sobj
# 
# 
# #Add metadata
# meta_cell<-data.frame(cell=colnames(sobj), individual="C97_caud3pr")
# meta_cell_add<-merge(meta_cell, meta, by.x="individual", by.y="Sample_ID")
# meta_cell_add<-meta_cell_add[match(colnames(sobj), meta_cell_add$cell),]
# print(identical(meta_cell_add$cell, colnames(sobj)))
# rownames(meta_cell_add)<-meta_cell_add$cell
# sobj<- AddMetaData(sobj, meta_cell_add)
# 
# #Store percentage of mitocondrial and ribosomal reads, and visualise
# 
# 
# sobj <- PercentageFeatureSet(sobj, pattern = "^MT-", col.name = "percent.mt")
# sobj <- PercentageFeatureSet(sobj, pattern = "^RPL|^Rpl|^RPS|^Rps", col.name = "percent.ribo")
# # Visualize
# VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE)
# 
# 
# #Apply mitochondrial filter
# 
# 
# mitoCells <- WhichCells(sobj, expression = percent.mt < 50)
# # Only keep these cells
# sobj <- subset(sobj, cells = mitoCells)
# sobj
# 
# save(sobj, file=here("/media/redgar/Seagate Portable Drive/processed_data/empty_drop_C97_conservative.RData"))
# 
# load(here("/media/redgar/Seagate Portable Drive/processed_data/empty_drop_C97_conservative.RData"))
# 
# ##Run SCTransform
# 
# #
# # sobj <- SCTransform(sobj) %>%
# #   RunPCA() %>%
# #   FindNeighbors(dims = 1:30) %>%
# #   RunUMAP(dims = 1:30) %>%
# #   RunTSNE()
# 
# sobj <- NormalizeData(sobj)
# sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
# sobj <- ScaleData(sobj) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))
# 
# # dimension reduction
# sobj <- RunPCA(sobj, ndims.print = 1:10, nfeatures.print = 10)
# sobj <- RunUMAP(sobj, dims = 1:30)
# sobj <- RunTSNE(sobj, dims = 1:30)
# 
# sobj <- FindNeighbors(sobj, reduction = "pca", dims = 1:20)
# sobj <- FindClusters(sobj, resolution = 0.2)
# 
# 
# 
# # #Set up for scClustViz
# #
# #
# # scSeurat <- sobj
# # DE_bw_clust <- TRUE
# # FDRthresh <- 0.01
# # seurat_resolution <- 0
# # sCVdata_list <- list()
# #
# #
# # #Cluster
# #
# #
# # while(DE_bw_clust) {
# #   seurat_resolution <- seurat_resolution + 0.2
# #   # ^ iteratively incrementing resolution parameter
# #
# #   scSeurat <- FindClusters(scSeurat,
# #                            resolution=seurat_resolution,
# #                            print.output=F)
# #   message(" ")
# #   message("------------------------------------------------------")
# #   message(paste0("--------  res.",seurat_resolution," with ",
# #                  length(levels(Idents(scSeurat)))," clusters --------"))
# #   message("------------------------------------------------------")
# #   curr_sCVdata <- CalcSCV(inD=scSeurat,
# #                           cl=Idents(scSeurat),
# #                           # ^ your most recent clustering results get stored in the Seurat "ident" slot
# #                           assayType = "SCT",
# #                           exponent=exp(1),
# #                           pseudocount=1,
# #                           DRthresh=0.1,
# #                           DRforClust="pca",
# #                           calcSil=T,
# #                           calcDEvsRest=T,
# #                           calcDEcombn=T)
# #
# #   DE_bw_NN <- sapply(DEneighb(curr_sCVdata,FDRthresh),nrow)
# #   # ^ counts # of DE genes between neighbouring clusters at your selected FDR threshold
# #
# #   if (min(DE_bw_NN) < 1 | seurat_resolution > 2) { DE_bw_clust <- FALSE }
# #   # ^ If no DE genes between nearest neighbours, don't loop again.
# #
# # sCVdata_list[[paste0("res.",seurat_resolution)]] <- curr_sCVdata
# # }
# #
# #
# # Save object
# #
# #
# # save(sCVdata_list, scSeurat,
# #      file = paste(woodchuck, sample, "emptyDrops_lower", lower, "_scClustViz.RData", sep = "_"))
# #
# #
# # OR if the object is too large, run Seurat clustering loop
# #
# 
# # for (res in seq(from = 0.2, to = 2.0, by = 0.2)) {
# #   sobj <- FindClusters(sobj, resolution = res)
# # }
# # Idents(sobj) <- "SCT_snn_res.0.2"
# #
# #
# # Save object
# 
# 
# #saveRDS(sobj, file = paste(woodchuck, sample, "emptyDrops_lower", lower, "_sobj.RDS", sep = "_"))
# 
# 
# ### Visualisations
# 
# #Reassign object to new variable and load in `DropletQC` filtered object
# 
# #Add DropletQC result to `sobj`
# 
# load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
# 
# cell_label_C97<-cell_label[which(cell_label$individual=="C97_caud3pr"),]
# cell_label_C97$cell<-gsub("_14", "",cell_label_C97$cell)
# 
# cell_label_C97<-cell_label_C97[which(cell_label_C97$cell%in%colnames(sobj)),]
# meta_cell_add<-merge(sobj@meta.data[c("nCount_RNA","nFeature_RNA" , "individual","cell","percent.mt", "percent.ribo")], cell_label_C97[,c("cell","CellType_refined")], all.x=T, by="cell")
# 
# meta_cell_add<-meta_cell_add[match(colnames(sobj), meta_cell_add$cell),]
# print(identical(meta_cell_add$cell, colnames(sobj)))
# rownames(meta_cell_add)<-meta_cell_add$cell
# sobj<- AddMetaData(sobj, meta_cell_add)
# 
# 
# load("/home/redgar/Documents/liver_ped_map/data/QC_metrics.Rdata")
# plt_QC_C97<-plt_QC_data[which(plt_QC_data$individual=="C97_caud3pr"),]
# 
# plt_QC_C97<-plt_QC_C97[which(plt_QC_C97$cell%in%colnames(sobj)),]
# meta_cell_add<-merge(sobj@meta.data[c("nCount_RNA","nFeature_RNA" , "individual","cell","percent.mt", "percent.ribo")], plt_QC_C97, all.x=T, by="cell")
# 
# meta_cell_add<-meta_cell_add[match(colnames(sobj), meta_cell_add$cell),]
# print(identical(meta_cell_add$cell, colnames(sobj)))
# rownames(meta_cell_add)<-meta_cell_add$cell
# sobj<- AddMetaData(sobj, meta_cell_add)
# 
# 
# table(sobj$cell_status)
# 
# 
# #With `lower = 800`, datasets look very similar. Try increasing lower?
# 
# 
# status<-DimPlot(sobj, group.by = "cell_status")
# save_plts(status, "cell_status_human_C97", w=8, h=6)
# 
# fraction<-FeaturePlot(sobj, reduction="umap", features = "nuclear_fraction")
# save_plts(fraction, "nuclear_fraction_human_C97", w=8, h=6)
# 
# celltype<-DimPlot(sobj, group.by = "CellType_refined", label=T)+colscale_cellType
# save_plts(celltype, "celltype_human_C97", w=10, h=6)
# 
# 
# save(sobj, file=here("/media/redgar/Seagate Portable Drive/processed_data/empty_drop_C97_default.RData"))
# 
# load(here("/media/redgar/Seagate Portable Drive/processed_data/empty_drop_C97_default.RData"))
# sobj_default<-sobj
# load(here("/media/redgar/Seagate Portable Drive/processed_data/empty_drop_C97_conservative.RData"))
# sobj_800<-sobj
# rm(sobj)
# gc()
# 
# 
# 
# 
# sobj_default$cut_at_lower800<-"0"
# sobj_default$cut_at_lower800[which(!(sobj_default$cell%in%sobj_800$cell))]<-"1"
# table(sobj_default$cut_at_lower800)
# 
# 
# umap_mat_myeloid<-as.data.frame(Embeddings(object = sobj_default, reduction = "umap"))#
# umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
# meta_myeloid<-sobj_default@meta.data
# meta_myeloid$cell<-rownames(meta_myeloid)
# plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")
# 
# cut<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
#   geom_point(aes(fill=CellType_refined, color=cut_at_lower800),size=1, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
#   fillscale_cellType+theme_classic()+scale_color_manual(values=c("white","black"))
# 
# save_plts(cut, "cut_at_lower800_human_C97", w=10, h=6)
# 
# 



##################
## Multiple level comparison
##################

library(here)
library(DropletUtils)
library(Seurat)
library(dplyr)


#Read in raw data
dataset_loc <- here("/media/redgar/Seagate Portable Drive/ped_liver_map_raw")

samples<-list.files(dataset_loc)
print(samples)

meta<-read.table(here("data/data_transfer_updated_mar20_2023_IFALD.csv"), header=T, sep=",")
samples<-samples[grep("C68",samples)]

sobj.data <- Read10X(file.path(dataset_loc,paste(samples,"/outs", sep=""),"raw_feature_bc_matrix"))


#From these counts, make single cell experiment object
sobj <- CreateSeuratObject(counts = sobj.data)
sce <- as.SingleCellExperiment(sobj)
sce

rm(sobj.data)
gc()


cell_lower<-lapply(c(400, 800, 1000, 1200, 1400), function(lower){
  #Calculate empty drops with `lower = 800` since default `lower` returns tens of thousands of cells
  set.seed(100)
  print(lower)
  e.out <- emptyDrops(counts(sce), lower=lower)
  
  #Detect significant deviations from the ambient profile by setting FDR threshold
  
  # Keep cells that have an FDR of less than 1%
  is.cell <- e.out$FDR <= 0.01
  # Count the cells that are remaining
  print(sum(is.cell, na.rm=TRUE))
  
  #Subset object to keep only cells
  namedCounts <- counts(sce)
  colnames(namedCounts) <- rownames(sce@colData)
  cell.counts <- namedCounts[,which(is.cell),drop=FALSE]
  gc()
  
  #Convert to Seurat object, don't use `min.features` since `emptyDrops` was run.
  sobj <- CreateSeuratObject(counts = cell.counts, min.cells = 3)
  
  #Add metadata
  meta_cell<-data.frame(cell=colnames(sobj), individual="C68_caud3pr")
  meta_cell_add<-merge(meta_cell, meta, by.x="individual", by.y="Sample_ID")
  meta_cell_add<-meta_cell_add[match(colnames(sobj), meta_cell_add$cell),]
  print(identical(meta_cell_add$cell, colnames(sobj)))
  rownames(meta_cell_add)<-meta_cell_add$cell
  sobj<- AddMetaData(sobj, meta_cell_add)
  
  #Store percentage of mitocondrial and ribosomal reads, and visualise
  sobj <- PercentageFeatureSet(sobj, pattern = "^MT-", col.name = "percent.mt")
  sobj <- PercentageFeatureSet(sobj, pattern = "^RPL|^Rpl|^RPS|^Rps", col.name = "percent.ribo")
  
  #Apply mitochondrial filter
  mitoCells <- WhichCells(sobj, expression = percent.mt < 50) 
  # Only keep these cells
  sobj <- subset(sobj, cells = mitoCells)
  sobj
  gc()
  as.character(sobj$cell)})


load("/home/redgar/Documents/liver_ped_map/data/QC_metrics.Rdata")
plt_QC_C68<-plt_QC_data[which(plt_QC_data$individual=="C68_caud3pr"),]

lowers<-c(400, 800, 1000, 1200, 1400)
lapply(1:length(lowers), function(x){
  lowers[x]
  plt_QC_C68$lower<<-"Cut by emptyDrops"
  plt_QC_C68$lower[which(plt_QC_C68$cell%in%cell_lower[[x]])]<<-"Kept by emptyDrops"
  colnames(plt_QC_C68)[which(colnames(plt_QC_C68)=="lower")]<<-as.character(paste("lower",lowers[x], sep="_"))
})





library(ggsankey)

df <- plt_QC_C68 %>%
  make_long(cell_status, lower_400)

dagg <- df%>%
  dplyr::group_by(node)%>%
  tally()

df <- merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)


cut<-ggplot(df, aes(x = x,
                    next_x = next_x,
                    node = node,
                    next_node = next_node,
                    fill = factor(node),
                    label = n)) +
  geom_sankey() +xlab("")+ggtitle("Lower 400")+
  theme_sankey(base_size = 16)+
  geom_sankey(flow.alpha = 0.25,  color = "gray40", show.legend = TRUE)+
  geom_sankey_label(size = 3, color = "white", fill= "gray40", hjust = -0.2)+
  scale_fill_manual(values = c("cornflowerblue","grey40","grey","dodgerblue4"))
save_plts(cut, "sankey_lower400_human_C68", w=10, h=6)

df <- plt_QC_C68 %>%
  make_long(cell_status, lower_800)
dagg <- df%>%
  dplyr::group_by(node)%>%
  tally()

df <- merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)

cut<-ggplot(df, aes(x = x,
                    next_x = next_x,
                    node = node,
                    next_node = next_node,
                    fill = factor(node),
                    label = n)) +
  geom_sankey() +xlab("")+ggtitle("Lower 800")+
  theme_sankey(base_size = 16)+
  geom_sankey(flow.alpha = 0.25,  color = "gray40", show.legend = TRUE)+
  geom_sankey_label(size = 3, color = "white", fill= "gray40", hjust = -0.2)+
  scale_fill_manual(values = c("cornflowerblue","grey40","grey","dodgerblue4"))
save_plts(cut, "sankey_lower800_human_C68", w=10, h=6)


df <- plt_QC_C68 %>%
  make_long(cell_status, lower_1000)
dagg <- df%>%
  dplyr::group_by(node)%>%
  tally()

df <- merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)

cut<-ggplot(df, aes(x = x,
                    next_x = next_x,
                    node = node,
                    next_node = next_node,
                    fill = factor(node),
                    label = n)) +
  geom_sankey() +xlab("")+ggtitle("Lower 1000")+
  theme_sankey(base_size = 16)+
  geom_sankey(flow.alpha = 0.25,  color = "gray40", show.legend = TRUE)+
  geom_sankey_label(size = 3, color = "white", fill= "gray40", hjust = -0.2)+
  scale_fill_manual(values = c("cornflowerblue","grey40","grey","dodgerblue4"))
save_plts(cut, "sankey_lower1000_human_C68", w=10, h=6)


df <- plt_QC_C68 %>%
  make_long(cell_status, lower_1200)
dagg <- df%>%
  dplyr::group_by(node)%>%
  tally()

df <- merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)

cut<-ggplot(df, aes(x = x,
                    next_x = next_x,
                    node = node,
                    next_node = next_node,
                    fill = factor(node),
                    label = n)) +
  geom_sankey() +xlab("")+ggtitle("Lower 1200")+
  theme_sankey(base_size = 16)+
  geom_sankey(flow.alpha = 0.25,  color = "gray40", show.legend = TRUE)+
  geom_sankey_label(size = 3, color = "white", fill= "gray40", hjust = -0.2)+
  scale_fill_manual(values = c("cornflowerblue","grey40","grey","dodgerblue4"))
save_plts(cut, "sankey_lower1200_human_C68", w=10, h=6)

df <- plt_QC_C68 %>%
  make_long(cell_status, lower_1400)
dagg <- df%>%
  dplyr::group_by(node)%>%
  tally()

df <- merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)

cut<-ggplot(df, aes(x = x,
                    next_x = next_x,
                    node = node,
                    next_node = next_node,
                    fill = factor(node),
                    label = n)) +
  geom_sankey() +xlab("")+ggtitle("Lower 1400")+
  theme_sankey(base_size = 16)+
  geom_sankey(flow.alpha = 0.25,  color = "gray40", show.legend = TRUE)+
  geom_sankey_label(size = 3, color = "white", fill= "gray40", hjust = -0.2)+
  scale_fill_manual(values = c("cornflowerblue","grey40","grey","dodgerblue4"))
save_plts(cut, "sankey_lower1400_human_C68", w=10, h=6)


