# library('variancePartition')
# library('BiocParallel')
# library('plyr')

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

library(edgeR)


d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))

load(here("data","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

# for (celltype in cts){
#   print(celltype)

d10x_raw_KC<-subset(d10x, subset = CellType_refined %in% c("KC Like"))
rm(d10x)

## filter genes with low counts for efficiency
min_percent_cells = 10  # Adjust as needed
num_cells_with_counts <- Matrix::rowSums(d10x_raw_KC@assays$RNA@counts > 0)
mini <- d10x_raw_KC@assays$RNA@counts[num_cells_with_counts >= (min_percent_cells / 100) * ncol(d10x_raw_KC@assays$RNA@counts),]


df<-as.data.frame(t(as.matrix(mini)))
df$cell<-rownames(df)
df<-merge(df, d10x_raw_KC@meta.data[,c("individual","cell")], by="cell")
df$cell<-NULL

pseudobulk_count <- df %>%
  group_by(individual) %>%
  summarize(across(!starts_with("individual"), sum))

pseudobulk_count<-as.data.frame(pseudobulk_count)
rownames(pseudobulk_count)<-pseudobulk_count$individual
pseudobulk_count$individual<-NULL

countMatrix <- as.matrix(t(pseudobulk_count))

metadata <- d10x_raw_KC@meta.data[,c("individual","Treatment","Tissue","chemistry","Sex","Age","AgeGroup", "FreshorFrozen", "Perfused", "BMI","age_condition")]
metadata<-metadata[!duplicated(metadata),]
metadata<-metadata[match(colnames(countMatrix),metadata$individual),]

# if (dim(metadata)[1] < 5){
#   print("less than 5 samples detected. skipping this celltype")
#   next
# }



metadata$Sex <- as.factor(metadata$Sex) # ensure it is a factor

metadata$age_condition <- as.factor(metadata$age_condition)
metadata$age_condition <- relevel(metadata$age_condition, ref="Ped Healthy")

# Create a DGEList object
dge = DGEList( countMatrix )

dge$design <- model.matrix(~ age_condition + Sex, data = metadata)

# Filter out low-count genes and normalize
keep <- rowSums(cpm(dge) > 1) >= 2
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

# Estimate dispersion
dge <- estimateDisp(dge)

# Fit a negative binomial model with covariates
fit <- glmFit(dge, design = dge$design)

# Perform differential expression analysis
de_results <- glmLRT(fit, coef = "age_conditionAdult Healthy") 
top_de_genes <- topTags(de_results, n = Inf)$table  


top_de_genes[which(top_de_genes$logFC>0),][1:10,]
top_de_genes[which(top_de_genes$FDR<0.05 & top_de_genes$logFC>0),]
top_de_genes[which(top_de_genes$FDR<0.05 & top_de_genes$logFC<0),]

top_de_genes[grep("CCL",rownames(top_de_genes)),]




gene<-"CCL3"
gene<-"CLEC10A"
gene<-"IL1B"

gene<-"KMT2B"
gene<-"SAA2"
gene<-"DHCR24"
gene<-"GNG10"


plt<-metadata
plt$gene<-countMatrix[which(rownames(countMatrix)==gene),]

ggplot(plt, aes(age_condition, gene))+geom_boxplot(outlier.shape = NA)+geom_point(position = position_jitter(w=0.25))
ggplot(plt, aes(age_condition, gene))+geom_boxplot()+geom_text(aes(label=individual))
