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
library(FLASHMM)


source("scRNA_seq_scripts/00_pretty_plots.R")
source("scRNA_seq_scripts/00_fanciest_UMAP.R")
source("scRNA_seq_scripts/00_plot_gene_exp.R")
source("scRNA_seq_scripts/00_entropy_d10x.R")



####################################################################################
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_d10x_adult_ped_raw.rds"))

load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
d10x$Sex[which(d10x$individual%in%c("C113","C115"))]<-"M"

d10x_raw_KC<-subset(d10x, subset = CellType_refined %in% c("KC Like"))

d10x_raw_KC_healthy<-subset(d10x_raw_KC, subset = age_condition %in% c("Adult Healthy", "Ped Healthy"))


## get counts
counts <- as.matrix(d10x_raw_KC_healthy@assays$RNA@counts)
## library size to metadata
metadata <- d10x_raw_KC_healthy@meta.data
metadata$libsize<-colSums(counts)

metadata$cell_id<-rownames(metadata)

Ped_cells <- metadata %>% filter(age_condition == "Ped Healthy") %>% pull(cell_id)
Adult_cells <- metadata %>% filter(age_condition == "Adult Healthy") %>% pull(cell_id)

counts_ped <- counts[, Ped_cells]
counts_adult <- counts[, Adult_cells]

zero_percent_ped <- rowSums(counts_ped == 0) / ncol(counts_ped)
zero_percent_adult <- rowSums(counts_adult == 0) / ncol(counts_adult)

# Keep genes expressed in >10% of cells in both sge groups
gene_to_test <- intersect(
  names(zero_percent_ped)[zero_percent_ped < 0.9],
  names(zero_percent_adult)[zero_percent_adult < 0.9]
)

counts<-counts[gene_to_test,]

metadata$age_condition<-as.character(metadata$age_condition)

##Gene expression matrix, Y = log2(1 + counts)
Y <- log2(1 + counts)

##Design matrix for fixed effects
#X <- model.matrix(~ 0 + log(libsize) + cls + cls:trt, data = meta)
X <- model.matrix(~ 0 + log(libsize) + age_condition + Sex, data = metadata)

##Design matrix for random effects
#Z <- model.matrix(~ 0 + as.factor(sam), data = dat$meta)
Z <- model.matrix(~ 0 + as.factor(file), data = metadata)

##Dimension of random effects
d <- ncol(Z)


max.iter <- 100
epsilon <- 1e-5
fit <- lmmfit(Y, X, Z, d = d, max.iter = max.iter, epsilon = epsilon)

gc <- (apply(abs(fit$dlogL), 2, max) < epsilon) 
sum(gc) 


## contrast for sex
contrast <- rep(0, ncol(X))
names(contrast) <- colnames(X)
contrast["age_conditionAdult Healthy"] <- 1
contrast["age_conditionPed Healthy"] <- -1

contrast_age <- lmmtest(fit, contrast = contrast)
contrast_age<-as.data.frame(contrast_age)


##(3) The DE genes specific to a cell-type
##Coefficients, t-values, and p-values for the genes specific to a cell-type.
index <- grep("age_conditionPed Healthy", rownames(fit$coef))
ce <- fit$coef[index, gc]
tv <- fit$t[index, gc]
pv <- fit$p[index, gc]

contrast_age<-contrast_age[which(rownames(contrast_age)%in%names(ce)),]
identical(rownames(contrast_age),names(ce))

out <- data.frame(
  gene = names(ce), 
  coef = contrast_age$`_coef` , p = contrast_age$`_p`)

##Adjust p-values by FDR.
out$FDR <- p.adjust(out$p, method = "fdr")

##The DE genes with FDR < 0.05 and abs(logFC) > 0.5
out <- out[order(out$p), ]
rownames(out) <- NULL
out[(out$FDR < 0.05) & (abs(out$coef) > 0.5) , ]

sig_de<-out[(out$FDR < 0.05) & (abs(out$coef) > 0.5) , ]
nrow(sig_de)

sig_de_flash025<-out[(out$FDR < 0.05) & (abs(out$coef) > 0.25) , ]
nrow(sig_de_flash025)

#### check key genes
sig_de[grep("CCL", sig_de$gene),]
out[grep("CCL", out$gene),]


VlnPlot(d10x_raw_KC_healthy, features = c("RPS29","NDUFB1"), group.by = "age_condition")
VlnPlot(d10x_raw_KC_healthy, features = sig_de$gene[1:10], group.by = "age_condition")
VlnPlot(d10x_raw_KC_healthy, features = c("CCL4","IL1B"), group.by = "age_condition")

ggplot(out, aes(coef, -log10(p)))+geom_point()

save(out, file="data/flashMM_KC_age.RData")

load("data/flashMM_KC_age.RData")
sig_de_flash025<-out[(out$FDR < 0.05) & (abs(out$coef) > 0.25) , ]

##############
## MAST
##############
d10x<-readRDS(file = here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_d10x_adult_ped_raw.rds"))

load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
d10x$Sex[which(d10x$individual%in%c("C113","C115"))]<-"M"

d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")
d10x_raw_KC<-subset(d10x, subset = CellType_refined %in% c("KC Like"))
rm(d10x)
gc()

Idents(d10x_raw_KC)<-d10x_raw_KC$age_condition
de<-FindMarkers(d10x_raw_KC, ident.1 = "Adult Healthy", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA","Sex"), verbose=F)
sig_de<-de[which(de$p_val_adj < 0.005 & abs(de$avg_log2FC) > 1),]


#################
## Compare
#################
out_in_MAST<-out[which(out$gene%in%rownames(de)),]
MAST_in_out<-de[which(rownames(de)%in%out$gene),]
MAST_in_out$gene<-rownames(MAST_in_out)

compare_de<-merge(MAST_in_out, out_in_MAST, by="gene")
cor_r<-cor.test(compare_de$avg_log2FC, compare_de$coef)

compare_de$Sig<-"Neither"
compare_de$Sig[which(compare_de$gene%in%rownames(sig_de))]<-"Only MAST"
compare_de$Sig[which(compare_de$gene%in%sig_de_flash025$gene)]<-"Only FLASH"
compare_de$Sig[which(compare_de$gene%in%intersect(rownames(sig_de), sig_de_flash025$gene))]<-"Both"


de_combine_plt_cor<-ggplot(compare_de, aes(avg_log2FC,coef))+
  geom_vline(xintercept = c(1,-1), color="grey")+stat_smooth(method="lm", se=F, color="black")+
  geom_point(shape=21, aes(fill=Sig))+theme_bw()+
  scale_fill_manual(values=c("#f09c16","grey","#7fc97f","#386cb0"),name="Significant in which cohort")+
  theme(legend.position = "bottom")+guides(fill=guide_legend(ncol=2))+
  xlab("Fold Change (log2)\nFull (5' and 3') Cohort")+  ylab("Fold Change (log2)\n3' Samples Only")+
  geom_hline(yintercept = c(1,-1), color="grey")+xlim(-4,4)+ylim(-4,4)+
  annotate("text", x=3, y=-3, label= paste("italic(r[s]) ==  ", round(cor_r$estimate,2)),
           parse = TRUE)
de_combine_plt_cor

####################################################################################
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_d10x_adult_ped_raw.rds"))

load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)
d10x$Sex[which(d10x$individual%in%c("C113","C115"))]<-"M"

d10x_raw_KC<-subset(d10x, subset = CellType_refined %in% c("KC Like"))

d10x_raw_KC_ped<-subset(d10x_raw_KC, subset = age_condition %in% c("Ped IFALD", "Ped Healthy"))


## get counts
counts <- as.matrix(d10x_raw_KC_ped@assays$RNA@counts)
## library size to metadata
metadata <- d10x_raw_KC_ped@meta.data
metadata$libsize<-colSums(counts)

metadata$cell_id<-rownames(metadata)

Ped_cells <- metadata %>% filter(age_condition == "Ped Healthy") %>% pull(cell_id)
IFLAD_cells <- metadata %>% filter(age_condition == "Ped IFALD") %>% pull(cell_id)

counts_ped <- counts[, Ped_cells]
counts_IFALD <- counts[, IFLAD_cells]

zero_percent_ped <- rowSums(counts_ped == 0) / ncol(counts_ped)
zero_percent_IFALD <- rowSums(counts_IFALD == 0) / ncol(counts_IFALD)

# Keep genes expressed in >10% of cells in both sge groups
gene_to_test <- intersect(
  names(zero_percent_ped)[zero_percent_ped < 0.9],
  names(zero_percent_IFALD)[zero_percent_IFALD < 0.9]
)

counts<-counts[gene_to_test,]

metadata$age_condition<-as.character(metadata$age_condition)

##Gene expression matrix, Y = log2(1 + counts)
Y <- log2(1 + counts)

##Design matrix for fixed effects
#X <- model.matrix(~ 0 + log(libsize) + cls + cls:trt, data = meta)
X <- model.matrix(~ 0 + log(libsize) + age_condition + Sex, data = metadata)

##Design matrix for random effects
#Z <- model.matrix(~ 0 + as.factor(sam), data = dat$meta)
Z <- model.matrix(~ 0 + as.factor(file), data = metadata)

##Dimension of random effects
d <- ncol(Z)


max.iter <- 100
epsilon <- 1e-5
fit <- lmmfit(Y, X, Z, d = d, max.iter = max.iter, epsilon = epsilon)

gc <- (apply(abs(fit$dlogL), 2, max) < epsilon) 
sum(gc) 


## contrast for sex
contrast <- rep(0, ncol(X))
names(contrast) <- colnames(X)
contrast["age_conditionAdult Healthy"] <- 1
contrast["age_conditionPed Healthy"] <- -1

contrast_age <- lmmtest(fit, contrast = contrast)
contrast_age<-as.data.frame(contrast_age)


##(3) The DE genes specific to a cell-type
##Coefficients, t-values, and p-values for the genes specific to a cell-type.
index <- grep("age_conditionPed Healthy", rownames(fit$coef))
ce <- fit$coef[index, gc]
tv <- fit$t[index, gc]
pv <- fit$p[index, gc]

contrast_age<-contrast_age[which(rownames(contrast_age)%in%names(ce)),]
identical(rownames(contrast_age),names(ce))

out <- data.frame(
  gene = names(ce), 
  coef = contrast_age$`_coef` , t = tv, p = pv)

##Adjust p-values by FDR.
out$FDR <- p.adjust(out$p, method = "fdr")

##The DE genes with FDR < 0.05 and abs(logFC) > 0.5
out <- out[order(out$p), ]
rownames(out) <- NULL
out[(out$FDR < 0.05) & (abs(out$coef) > 0.5) , ]

sig_de<-out[(out$FDR < 0.05) & (abs(out$coef) > 0.5) , ]
nrow(sig_de)

sig_de<-out[(out$FDR < 0.05) & (abs(out$coef) > 0.25) , ]
nrow(sig_de)

#### check key genes
sig_de[grep("CCL", sig_de$gene),]
out[grep("CCL", out$gene),]


VlnPlot(d10x_raw_KC_healthy, features = c("CAP1","DDIT4"), group.by = "age_condition")
VlnPlot(d10x_raw_KC_healthy, features = sig_de$gene[1:10], group.by = "age_condition")


MAST_de<-read.csv(here("data/differential_age_KC.csv"))

sig_de<-out[(out$FDR < 0.05) , ]

MAST_de_tested<-MAST_de[which(MAST_de$X%in%out$gene),]

length(which(MAST_de_tested$X%in%sig_de$gene))

