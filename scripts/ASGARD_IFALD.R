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
## IFALD data
################
d10x<-readRDS(file = here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_d10x_adult_ped_raw.rds"))
load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

d10x_raw_IFALD_PED<-subset(d10x, subset = age_condition %in% c("Ped Healthy","Ped IFALD"))
rm(d10x)
gc()

d10x_raw_IFALD_PED <- NormalizeData(d10x_raw_IFALD_PED,scale.factor = 10000, normalization.method = "LogNormalize")


#Get differential genes from Seurat (Wilcoxon Rank Sum test)
DefaultAssay(d10x_raw_IFALD_PED) <- "RNA"
set.seed(42)
Gene.list <- list()
C_names <- NULL
for(i in unique(d10x_raw_IFALD_PED@meta.data$CellType_refined)){
  Idents(d10x_raw_IFALD_PED) <- "CellType_refined"
  c_cells <- subset(d10x_raw_IFALD_PED, CellType_refined == i)
  Idents(c_cells) <- "age_condition"
  C_data <- FindMarkers(c_cells, ident.1 = "Ped IFALD", ident.2 = "Ped Healthy")
  C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$avg_log2FC,adj.P.Val=C_data$p_val_adj,P.Value=C_data$p_val) ##for Seurat version > 4.0, please use avg_log2FC instead of avg_logFC
  Gene.list[[i]] <- C_data_for_drug
  C_names <- c(C_names,i)
  gc()
}
names(Gene.list) <- C_names

save(Gene.list, file="/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/DGE_for_ASGARD.RData")




#############
## Mono-drug repurposing for every cell type
#############
#Load tissue specific drug reference produced by PrepareReference function as mentioned above. Please select proper tissue accroding to the disease.
my_gene_info<-read.table(file="DrugReference/liver_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info<-read.table(file="DrugReference/liver_drug_info.txt",sep="\t",header = T,quote = "")
drug.ref.profiles = GetDrugRef(drug.response.path = 'DrugReference/liver_rankMatrix.txt',
                               probe.to.genes = my_gene_info, 
                               drug.info = my_drug_info)

#Repurpose mono-drugs for every cell type                               
Drug.ident.res = GetDrug(gene.data = Gene.list, 
                         drug.ref.profiles = drug.ref.profiles, 
                         repurposing.unit = "drug", 
                         connectivity = "negative", 
                         drug.type = "FDA")


####################
#Select mono-drug therapies
####################
Final.drugs<-subset(Drug.score,Drug.therapeutic.score>quantile(Drug.score$Drug.therapeutic.score, 0.99,na.rm=T) & FDR <0.05)


#Select drug for individual clusters
Final.drugs<-TopDrug(SC.integrated=SC.integrated,
                     Drug.data=Drug.ident.res,
                     Drug.FDR=0.1,
                     FDA.drug.only=TRUE,
                     Case=Case.samples,
                     DrugScore=FALSE
)
