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



meta<-read.table(here("data/data_transfer_updated_jan16_2023.csv"), header=T, sep=",")

dataset_loc <- here("/media/redgar/Seagate Portable Drive/ped_liver_map_raw")
samples<-list.files(dataset_loc)
print(samples)
y=4

caud<-meta$Sample_ID[which(meta$file == samples[y])]
print(caud)
print(file.path(dataset_loc,paste(samples[y],"/outs", sep=""),"filtered_feature_bc_matrix"))
d10x <- Read10X(file.path(dataset_loc,paste(samples[y],"/outs", sep=""),"filtered_feature_bc_matrix"))
colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),caud,sep="-")
# print(dim(d10x))
#' Initialize the Seurat object with the raw (non-normalized data).
d10x<-CreateSeuratObject(counts = d10x, project = "ped_adult_map", min.cells = 0, min.features = 0)
d10x

## SoupX needs clusters so quickly make clusters for each sample
d10x <- NormalizeData(d10x)
d10x <- FindVariableFeatures(d10x, selection.method = "vst", nfeatures = 2000)
d10x <- ScaleData(d10x) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))
d10x    <- RunPCA(d10x, verbose = F)
d10x    <- RunUMAP(d10x, dims = 1:30, verbose = F)
d10x    <- FindNeighbors(d10x, dims = 1:30, verbose = F)
d10x    <- FindClusters(d10x, verbose = T)
meta_clusters    <- d10x@meta.data


toc=Read10X(file.path(dataset_loc,paste(samples[y],"/outs/filtered_feature_bc_matrix", sep="")))
tod = Read10X(file.path(dataset_loc,paste(samples[y],"/outs/raw_feature_bc_matrix", sep="")))
sc=SoupChannel(tod,toc)

#Adding extra metadata to the SoupChannel Object
sc = setClusters(sc, setNames(meta_clusters$seurat_clusters, rownames(meta_clusters)))

# projection <- read.csv("HCC161_Normal_5pr/analysis/umap/2_components/projection.csv", row.names=1)
# sc = setDR(sc,projection)

# sc_metadata=cbind(projection, clusters)

# #Visual Sanity checks
# dd = sc_metadata[colnames(sc$toc),]
# dd$PTPRC = sc$toc['PTPRC',]
# gg = ggplot(dd,aes(UMAP.1,UMAP.2)) +
#   geom_point(aes(colour=IGKC>0))
# plot(gg)
# 
# gg = plotMarkerMap(sc,'PTPRC')
# plot(gg)


nonExpressedGeneList = list(hemoglobin=c('HBB',"HBA1",'HBA2'),IG = c('IGKC', "IGHE", "IGHM"))
head(sc$soupProfile[order(sc$soupProfile$est,decreasing=TRUE),],n=20)

plotMarkerDistribution(sc)

useToEst = estimateNonExpressingCells(sc,nonExpressedGeneList)
sc = calculateContaminationFraction(sc,nonExpressedGeneList,useToEst=useToEst)
head(sc$metaData)

out = adjustCounts(sc)
d10x_soupx = CreateSeuratObject(out)


cntSoggy = rowSums(sc$toc>0)
cntStrained = rowSums(out>0)
mostZeroed = tail(sort((cntSoggy-cntStrained)/cntSoggy),n=10)
mostZeroed

tail(sort(rowSums(sc$toc>out)/rowSums(sc$toc>0)),n=30)

DropletUtils:::write10xCounts('./strainedCounts',out)