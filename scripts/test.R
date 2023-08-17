
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
library(SCINA)




source("scripts/00_pretty_plots.R")
source("scripts/00_entropy_d10x.R")
source("scripts/00_fanciest_UMAP.R")

#################
## SCINA
#################
#d10x<-readRDS(file = here("data","IFALD_d10x_adult_ped_raw.rds"))
d10x<-readRDS(file = here("../../../projects/macparland/RE/PediatricAdult/processed_data","IFALD_d10x_adult_ped_raw.rds"))

d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

signatures<-read.csv(here("data/Liver_Markers_with_citations - Human_for_SCINA.csv"))

d10x_exp <- GetAssayData(d10x)
results = SCINA(d10x_exp, signatures, max_iter = 100, convergence_n = 10,
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
d10x$SCINA_broad<-results$cell_labels
gc()

no_refined<-data.frame(cell=colnames(d10x), SCINA_broad=results$cell_labels, SCINA_refined=NA)


### Cell subsets
RBCsignatures<-read.csv(here("data/Liver_Markers_with_citations - Erythrocytes.csv"))
Neutrosignatures<-read.csv(here("data/Liver_Markers_with_citations - Neutrophil.csv"))
Tsignatures<-read.csv(here("data/Liver_Markers_with_citations - Tcell.csv"))
Bsignatures<-read.csv(here("data/Liver_Markers_with_citations - Bcell.csv"))
myeloidsignatures<-read.csv(here("data/Liver_Markers_with_citations - Myeloid.csv"))

# d10x.combined_RBC<-subset(d10x, subset = SCINA_broad == "Erythrocytes")
# d10x_exp <- GetAssayData(d10x.combined_RBC)
# results = SCINA(d10x_exp, RBCsignatures, max_iter = 100, convergence_n = 10,
#                 convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
# d10x.combined_RBC$SCINA_refined<-results$cell_labels

d10x.combined_myeloid<-subset(d10x, subset = SCINA_broad == "Myeloid")
d10x_exp <- GetAssayData(d10x.combined_myeloid)
results = SCINA(d10x_exp, myeloidsignatures, max_iter = 100, convergence_n = 10,
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
d10x.combined_myeloid$SCINA_refined<-results$cell_labels

d10x.combined_bcell<-subset(d10x, subset = SCINA_broad == "B_cell")
d10x_exp <- GetAssayData(d10x.combined_bcell)
results = SCINA(d10x_exp, Bsignatures, max_iter = 100, convergence_n = 10,
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
d10x.combined_bcell$SCINA_refined<-results$cell_labels

d10x.combined_tcell<-subset(d10x, subset = SCINA_broad == "T_cell")
d10x_exp <- GetAssayData(d10x.combined_tcell)
results = SCINA(d10x_exp, Tsignatures, max_iter = 100, convergence_n = 10,
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
d10x.combined_tcell$SCINA_refined<-results$cell_labels

# d10x.combined_neutro<-subset(d10x, subset = SCINA_broad == "Neutrophil")
# d10x_exp <- GetAssayData(d10x.combined_neutro)
# results = SCINA(d10x_exp, Neutrosignatures, max_iter = 100, convergence_n = 10,
#                 convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
# d10x.combined_neutro$SCINA_refined<-results$cell_labels



SCINA_cell_labels<-rbind(#d10x.combined_RBC@meta.data[,c("cell","SCINA_broad","SCINA_refined")],d10x.combined_neutro@meta.data[,c("cell","SCINA_broad","SCINA_refined")]
                         d10x.combined_myeloid@meta.data[,c("cell","SCINA_broad","SCINA_refined")],
                         d10x.combined_bcell@meta.data[,c("cell","SCINA_broad","SCINA_refined")],
                         d10x.combined_tcell@meta.data[,c("cell","SCINA_broad","SCINA_refined")])
SCINA_cell_labels$cell<-rownames(SCINA_cell_labels)

SCINA_cell_labels<-rbind(SCINA_cell_labels, no_refined[which(!(no_refined$cell%in%SCINA_cell_labels$cell)),])
save(SCINA_cell_labels, file=here("data","IFALD_adult_ped_SCINA_markers_withcitations_cell_labels.RData"))

