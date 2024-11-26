### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(reshape2)
library(gtools)
library(colorspace)
library(cowplot)
library(RColorBrewer)


source("scripts/00_pretty_plots.R")
source("scripts/00_plot_gene_exp.R")
source("scripts/00_fanciest_UMAP.R")



#################
## ASGARD
#################
#devtools::install_github("lanagarmire/Asgard")

library('Asgard')
library(Hmisc)


################
## download drug libraries
################

# cd /media/redgar/Seagate\ Portable\ Drive/L1000
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_cell_info.txt.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz

PrepareReference(cell.info="/media/redgar/Seagate Portable Drive/L1000/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt",
                 gene.info="/media/redgar/Seagate Portable Drive/L1000/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt",
                 GSE70138.sig.info = "/media/redgar/Seagate Portable Drive/L1000/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt",
                 GSE92742.sig.info = "/media/redgar/Seagate Portable Drive/L1000/GSE92742_Broad_LINCS_sig_info.txt",
                 GSE70138.gctx = "/media/redgar/Seagate Portable Drive/L1000/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx",
                 GSE92742.gctx = "/media/redgar/Seagate Portable Drive/L1000/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",
                 Output.Dir = "/media/redgar/Seagate Portable Drive/L1000/DrugReference/"
)


# on compute canada
# cd scratch
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_cell_info.txt.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz


PrepareReference(cell.info="/home/redgar25/scratch/L1000/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt",
                 gene.info="/home/redgar25/scratch/L1000/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt",
                 GSE70138.sig.info = "/home/redgar25/scratch/L1000/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt",
                 GSE92742.sig.info = "/home/redgar25/scratch/L1000/GSE92742_Broad_LINCS_sig_info.txt",
                 GSE70138.gctx = "/home/redgar25/scratch/L1000/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx",
                 GSE92742.gctx = "/home/redgar25/scratch/L1000/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",
                 Output.Dir = "/home/redgar25/scratch/L1000/DrugReference/"
)





################
## PCLS data
################
d10x<-readRDS(file = here("/media/redgar/Seagate Portable Drive/PCLS/PCLSIL4Day0to5_clean.RDS"))


d10x_day1<-subset(d10x, subset = Day %in% c("One"))
d10x_day1 <- NormalizeData(d10x_day1,scale.factor = 10000, normalization.method = "LogNormalize")

d10x_day3<-subset(d10x, subset = Day %in% c("Three"))
d10x_day3 <- NormalizeData(d10x_day3,scale.factor = 10000, normalization.method = "LogNormalize")

d10x_day5<-subset(d10x, subset = Day %in% c("Five"))
d10x_day5 <- NormalizeData(d10x_day5,scale.factor = 10000, normalization.method = "LogNormalize")

rm(d10x)
gc()

#Get differential genes from Seurat (Wilcoxon Rank Sum test)
DefaultAssay(d10x_day1) <- "RNA"
set.seed(42)
Gene.list <- list()
C_names <- NULL
for(i in unique(d10x_day1@meta.data$genlab)){
  Idents(d10x_day1) <- "genlab"
  c_cells <- subset(d10x_day1, genlab == i)
  Idents(c_cells) <- "Treatment"
  C_data <- FindMarkers(c_cells, ident.1 = "none", ident.2 = "IL_4")
  C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$avg_log2FC,adj.P.Val=C_data$p_val_adj,P.Value=C_data$p_val) ##for Seurat version > 4.0, please use avg_log2FC instead of avg_logFC
  Gene.list[[i]] <- C_data_for_drug
  C_names <- c(C_names,i)
  gc()
}
names(Gene.list.1) <- unique(d10x_day1@meta.data$genlab)
Gene.list.1<-Gene.list

#Get differential genes from Seurat (Wilcoxon Rank Sum test)
DefaultAssay(d10x_day3) <- "RNA"
set.seed(42)
Gene.list <- list()
C_names <- NULL
for(i in unique(d10x_day3@meta.data$genlab)[1:8]){
  Idents(d10x_day3) <- "genlab"
  c_cells <- subset(d10x_day3, genlab == i)
  Idents(c_cells) <- "Treatment"
  C_data <- FindMarkers(c_cells, ident.1 = "none", ident.2 = "IL_4")
  C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$avg_log2FC,adj.P.Val=C_data$p_val_adj,P.Value=C_data$p_val) ##for Seurat version > 4.0, please use avg_log2FC instead of avg_logFC
  Gene.list[[i]] <- C_data_for_drug
  C_names <- c(C_names,i)
  gc()
}
names(Gene.list.3) <- unique(d10x_day3@meta.data$genlab)[1:8]
Gene.list.3<-Gene.list

#Get differential genes from Seurat (Wilcoxon Rank Sum test)
DefaultAssay(d10x_day5) <- "RNA"
set.seed(42)
Gene.list <- list()
C_names <- NULL
for(i in unique(d10x_day5@meta.data$genlab)[1:8]){
  Idents(d10x_day5) <- "genlab"
  c_cells <- subset(d10x_day5, genlab == i)
  Idents(c_cells) <- "Treatment"
  C_data <- FindMarkers(c_cells, ident.1 = "none", ident.2 = "IL_4")
  C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$avg_log2FC,adj.P.Val=C_data$p_val_adj,P.Value=C_data$p_val) ##for Seurat version > 4.0, please use avg_log2FC instead of avg_logFC
  Gene.list[[i]] <- C_data_for_drug
  C_names <- c(C_names,i)
  gc()
}
names(Gene.list.5) <- unique(d10x_day5@meta.data$genlab)[1:8]
Gene.list.5<-Gene.list

save(Gene.list.1,Gene.list.3,Gene.list.5, file="/media/redgar/Seagate Portable Drive/PCLS/PCLS_DGE_for_ASGARD.RData")


###########
## convert DEG to human names
###########

# Basic function to convert mouse to human gene names
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(Orthology.eg.db)

mapIt <- function(mouseids, horg, morg, orth){
  mouseg <- mapIds(morg, mouseids, "ENTREZID", "SYMBOL")
  mapped <- select(orth, mouseg, "Homo_sapiens","Mus_musculus")
  names(mapped) <- c("Mus_egid","Homo_egid")
  husymb <- select(horg, as.character(mapped[,2]), "SYMBOL","ENTREZID")
  return(data.frame(Mus_symbol = mouseids,
                    mapped,
                    Homo_symbol = husymb[,2]))
}


Gene.list.1<-lapply(1:length(Gene.list.1), function(x){
  Gene.list.1[[x]]$gene<-rownames(Gene.list.1[[x]])
  
  gene_mapped<-mapIt(Gene.list.1[[x]]$gene, org.Hs.eg.db, org.Mm.eg.db, Orthology.eg.db)
  Gene.list.1[[x]]<-Gene.list.1[[x]][which(!(is.na(gene_mapped$Homo_symbol))),]
  
  rownames(Gene.list.1[[x]])<-gene_mapped$Homo_symbol[which(!(is.na(gene_mapped$Homo_symbol)))]
  Gene.list.1[[x]]$gene<-NULL
  Gene.list.1[[x]]
})

Gene.list.3<-lapply(1:length(Gene.list.3), function(x){
  Gene.list.3[[x]]$gene<-rownames(Gene.list.3[[x]])
  
  gene_mapped<-mapIt(Gene.list.3[[x]]$gene, org.Hs.eg.db, org.Mm.eg.db, Orthology.eg.db)
  Gene.list.3[[x]]<-Gene.list.3[[x]][which(!(is.na(gene_mapped$Homo_symbol))),]
  
  rownames(Gene.list.3[[x]])<-gene_mapped$Homo_symbol[which(!(is.na(gene_mapped$Homo_symbol)))]
  Gene.list.3[[x]]$gene<-NULL
  Gene.list.3[[x]]
})

Gene.list.5<-lapply(1:length(Gene.list.5), function(x){
  print(x)
  Gene.list.5[[x]]$gene<-rownames(Gene.list.5[[x]])
  
  gene_mapped<-mapIt(Gene.list.5[[x]]$gene, org.Hs.eg.db, org.Mm.eg.db, Orthology.eg.db)
  Gene.list.5[[x]]<-Gene.list.5[[x]][which(!(is.na(gene_mapped$Homo_symbol))),]
  
  rownames(Gene.list.5[[x]])<-gene_mapped$Homo_symbol[which(!(is.na(gene_mapped$Homo_symbol)))]
  Gene.list.5[[x]]$gene<-NULL
  Gene.list.5[[x]]
})



#############
## Mono-drug repurposing for every cell type
#############
#Load tissue specific drug reference produced by PrepareReference function as mentioned above. Please select proper tissue accroding to the disease.
my_gene_info<-read.table(file="/media/redgar/Seagate Portable Drive/L1000/DrugReference/liver_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info<-read.table(file="/media/redgar/Seagate Portable Drive/L1000/DrugReference/liver_drug_info.txt",sep="\t",header = T,quote = "")
drug.ref.profiles = GetDrugRef(drug.response.path = '/media/redgar/Seagate Portable Drive/L1000/DrugReference/liver_rankMatrix.txt',
                               probe.to.genes = my_gene_info, 
                               drug.info = my_drug_info)

#Repurpose mono-drugs for every cell type                               
Drug.ident.res = GetDrug(gene.data = Gene.list.1, 
                         drug.ref.profiles = drug.ref.profiles, 
                         repurposing.unit = "drug", 
                         connectivity = "negative", 
                         drug.type = "FDA")

head(Drug.ident.res[[which(names(Drug.ident.res)=="Hepatocyte")]])
head(Drug.ident.res[[which(names(Drug.ident.res)=="Cholangiocyte")]])
head(Drug.ident.res[[which(names(Drug.ident.res)=="Myeloid")]])

sig_drugs.1<-do.call(rbind, lapply(1:7, function(x){
  sig_drug<-Drug.ident.res[[x]][which(Drug.ident.res[[x]]$FDR<0.1),]
  if(nrow(sig_drug)>0){
    sig_drug$celltype<-names(Drug.ident.res)[x]
    sig_drug}
}))

head(sig_drugs.1)
sig_drugs.1$Drug.name


Drug.ident.res = GetDrug(gene.data = Gene.list.3, 
                         drug.ref.profiles = drug.ref.profiles, 
                         repurposing.unit = "drug", 
                         connectivity = "negative", 
                         drug.type = "FDA")

head(Drug.ident.res[[which(names(Drug.ident.res)=="Hepatocyte")]])
head(Drug.ident.res[[which(names(Drug.ident.res)=="Cholangiocyte")]])
head(Drug.ident.res[[which(names(Drug.ident.res)=="Myeloid")]])

sig_drugs.3<-do.call(rbind, lapply(1:7, function(x){
  sig_drug<-Drug.ident.res[[x]][which(Drug.ident.res[[x]]$FDR<0.1),]
  if(nrow(sig_drug)>0){
    sig_drug$celltype<-names(Drug.ident.res)[x]
    sig_drug}
}))

head(sig_drugs.3)
sig_drugs.3$Drug.name

Drug.ident.res = GetDrug(gene.data = Gene.list.5, 
                         drug.ref.profiles = drug.ref.profiles, 
                         repurposing.unit = "drug", 
                         connectivity = "negative", 
                         drug.type = "FDA")

head(Drug.ident.res[[which(names(Drug.ident.res)=="Hepatocyte")]])
head(Drug.ident.res[[which(names(Drug.ident.res)=="Cholangiocyte")]])
head(Drug.ident.res[[which(names(Drug.ident.res)=="Myeloid")]])

sig_drugs.5<-do.call(rbind, lapply(1:7, function(x){
  sig_drug<-Drug.ident.res[[x]][which(Drug.ident.res[[x]]$FDR<0.1),]
  if(nrow(sig_drug)>0){
    sig_drug$celltype<-names(Drug.ident.res)[x]
    sig_drug}
}))

head(sig_drugs.5)
sig_drugs.5$Drug.name








# 
# cell_metadata <- d10x_day3@meta.data
# cell_metadata$cluster <- d10x_day3@meta.data$genlab
# cell_metadata$sample<-cell_metadata$ID
# 
# Drug.score <- DrugScore(cell_metadata, cluster_degs = Gene.list, 
#                         cluster_drugs = Drug.ident.res, tissue = "liver", 
#                         case = c("D3_Ctrl"), 
#                         gse92742_gctx_path = "/media/redgar/Seagate Portable Drive/L1000/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx", 
#                         gse70138_gctx_path = "/media/redgar/Seagate Portable Drive/L1000/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx")
# 
# 
# Final.drugs<-subset(Drug.score,Drug.therapeutic.score>quantile(Drug.score$Drug.therapeutic.score, 0.99,na.rm=T) & FDR <0.05)
# 
# 
# #Select drug for individual clusters
# d10x_day3@meta.data$cluster<-d10x_day3@meta.data$genlab
# d10x_day3@meta.data$sample<-d10x_day3@meta.data$ID
# 
# Final.drugs<-TopDrug(SC.integrated=d10x_day3,
#                      Drug.data=Drug.ident.res,
#                      Drug.FDR=0.1,
#                      FDA.drug.only=TRUE,
#                      Case=c("D3_Ctrl"))
# 
