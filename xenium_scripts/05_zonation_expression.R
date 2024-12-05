## Load Libraries
library(here)
library(Seurat)

library(SCINA)
library(reshape2)
library(dplyr)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

library(sp)
library(rgeos)
library(scales)
library(viridis)
library(ggridges)

library(lme4)
library(tiff)


source("xenium_scripts/00_pretty_plots.R")
source("xenium_scripts/00_long_functions.R")



## the panel comes with cell type labels
load(file=here("data/cell_type_labels_BIDCell.RData"))





count_files<-c("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_2/cell_gene_matrices/2024_03_21_16_22_10/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_3/cell_gene_matrices/2024_03_21_16_38_48/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_4/cell_gene_matrices/2024_03_21_17_03_05/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C105/cell_gene_matrices/2024_05_17_13_32_03/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C85/cell_gene_matrices/2024_05_17_13_29_40/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C95/cell_gene_matrices/2024_07_04_12_35_29/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C101/cell_gene_matrices/2024_06_28_09_23_40/expr_mat.csv")


samples<-c("C94_2", "C94_3","C94_4","C105","C85","C95","C101")
models<-c("2024_03_21_16_22_10","2024_03_21_16_38_48","2024_03_21_17_03_05","2024_05_17_13_32_03","2024_05_17_13_29_40","2024_07_04_12_35_29","2024_06_28_09_23_40")


load(here("data","zonation_allcells.RData"))

######################
## Significance each sample indivudally
######################
de.list <- sapply(1:7, function(x){
  counts<-read.csv(count_files[x])
  print(samples[x])
  counts$X<-NULL
  rownames(counts) <- counts$cell_id
  counts$cell_id <- NULL

  # Transpose the data and convert to sparse matrix.
  mat <- as(t(as.matrix(counts)), "sparseMatrix")

  seu <- CreateSeuratObject(counts=mat)
  seu$sample<-strsplit(strsplit(count_files[x],"/")[[1]][6],"__")[[1]][3]
  seu

  ## zonation distance

  plt_zonation_xenium<-zonation_scores[[x]]

  plt_zonation_xenium<-plt_zonation_xenium[match(colnames(seu), plt_zonation_xenium$cell),]
  identical(colnames(seu), plt_zonation_xenium$cell)
  rownames(plt_zonation_xenium)<-plt_zonation_xenium$cell

  seu<-AddMetaData(seu, plt_zonation_xenium)
  seu <- NormalizeData(seu)

  ##################
  ## Differential Expression with portal distance
  ##################

  celltypes<-names(table(seu$CellType))[which(table(seu$CellType)>50)]

  lapply(1:length(celltypes), function(y){
    celltype<-celltypes[y]
    print(celltype)

    ## data looks pretty zero inflated so will fit a negative binomial
    ## normalize
    seu_celltype<-subset(seu, subset = CellType == celltype)
    gene_exp<-as.matrix(GetAssayData(object = seu_celltype, assay = "RNA", slot = "data"))

    meta<-seu_celltype@meta.data

    identical(rownames(meta), colnames(gene_exp))

    ## lm negative binomial on all genes
    pval_distance<-do.call(rbind, lapply(1:nrow(gene_exp), function(x) {
      gene<-rownames(gene_exp)[x]
      #print(gene)
      meta$exp<-gene_exp[which(rownames(gene_exp)==gene),]

      if((sum(meta$exp==0)/nrow(meta))<0.85){ #more that 15% of cell express
        # Negative binomial GLMM using the function glmer.nb()

        testFunction_portal <- function (meta) {return(tryCatch(glm.nb(exp ~ distance_to_periportal, data = meta), error=function(e) NULL)) }
        testFunction_central <- function (meta) {return(tryCatch(glm.nb(exp ~ distance_to_pericentral, data = meta), error=function(e) NULL)) }

        if(is.null(testFunction_portal(meta))&is.null(testFunction_central(meta))){
          data.frame(gene = gene,pval_periportal = "Errored",estimate_periportal = "Errored",pval_pericentral = "Errored",estimate_pericentral = "Errored")}else{
            if(is.null(testFunction_portal(meta))){
              mnb_pericentral <- glm.nb(exp ~ distance_to_pericentral, data = meta)
              data.frame(gene = gene,pval_periportal = "Errored",estimate_periportal = "Errored",
                         pval_pericentral = summary(mnb_pericentral)$coefficients["distance_to_pericentral","Pr(>|z|)"],
                         estimate_pericentral = summary(mnb_pericentral)$coefficients["distance_to_pericentral","Estimate"])}else{
                           if(is.null(testFunction_central(meta))){
                             mnb_periportal <- glm.nb(exp ~ distance_to_periportal, data = meta)
                             data.frame(gene = gene,
                                        pval_periportal = summary(mnb_periportal)$coefficients["distance_to_periportal","Pr(>|z|)"],
                                        estimate_periportal = summary(mnb_periportal)$coefficients["distance_to_periportal","Estimate"],
                                        pval_pericentral = "Errored",estimate_pericentral = "Errored")}else{
                                          mnb_periportal <- glm.nb(exp ~ distance_to_periportal, data = meta)
                                          mnb_pericentral <- glm.nb(exp ~ distance_to_pericentral, data = meta)
                                          data.frame(gene = gene,
                                                     pval_periportal = summary(mnb_periportal)$coefficients["distance_to_periportal","Pr(>|z|)"],
                                                     estimate_periportal = summary(mnb_periportal)$coefficients["distance_to_periportal","Estimate"],
                                                     pval_pericentral = summary(mnb_pericentral)$coefficients["distance_to_pericentral","Pr(>|z|)"],
                                                     estimate_pericentral = summary(mnb_pericentral)$coefficients["distance_to_pericentral","Estimate"])}}}


      }
    }))

    pval_distance$p_adjusted_periportal<-p.adjust(pval_distance$pval_periportal, method = "bonferroni")
    pval_distance$p_adjusted_pericentral<-p.adjust(pval_distance$pval_pericentral, method = "bonferroni")
    pval_distance$celltype<-celltype

    write.csv(pval_distance, file=paste(here("data/"),samples[x],"_", gsub("/","_",celltype),"_DEG_distance.csv",sep=""), quote=F, row.names = F)
  })


 })


####################
## interpretation
####################
cell_markers_general<-read.csv(here("data/Xenium_CombinedPanel.csv"))
custom_markers<-read.csv(here("data/Xenium_MacParlandGeneListUpdated.csv"))
custom_markers$Gene[grep("HAL", custom_markers$Gene)]<-"HAL"

xenium<-read.csv(file=here("/home/redgar/Documents/xenium_liver/data/Xenium_CombinedPanel.csv"))


sig_de_age_RR<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_RR.csv"))
sig_de_age_KC<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_KC.csv"))
sig_de_age_MHCII<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_MHCII.csv"))

sig_de_age_KC[which(sig_de_age_KC$X%in%xenium$Gene),]
sig_de_age_RR[which(sig_de_age_RR$X%in%xenium$Gene),]
sig_de_age_MHCII[which(sig_de_age_MHCII$X%in%xenium$Gene),]

intersect(xenium$Gene, sig_de_age_KC$X)
intersect(xenium$Gene, sig_de_age_RR$X)
intersect(xenium$Gene, sig_de_age_MHCII$X)

ped_map_KC<-intersect(xenium$Gene, sig_de_age_KC$X)
ped_map_RR<-intersect(xenium$Gene, sig_de_age_RR$X)
ped_map_MHCII<-intersect(xenium$Gene, sig_de_age_MHCII$X)



unique(zonation_scores[[1]]$CellType)

DEG_celltype<-function(celltype){
  do.call(rbind,
          lapply(samples[2:7], function(sample){
            DEG<-read.csv(file=here("data", paste(sample,"_",celltype,"_DEG_distance.csv", sep="")))
            DEG$sample<-sample
            DEG
          }))}

inped_notadult<-function(deg){
  inboth_ped<-names(which(table(deg[which(deg$sample%in%c("C105","C85")),]$gene)>1))
  adult<-deg[which(!(deg$sample%in%c("C105","C85"))),]$gene
  inboth_ped[which(!(inboth_ped%in%adult))]
}
inadult_notped<-function(deg){
  in_ped<-deg[which(deg$sample%in%c("C105","C85")),]$gene
  adult<-deg[which(!(deg$sample%in%c("C105","C85"))),]$gene
  adult[which(!(adult%in%in_ped))]
}

MHCII_DEG<-DEG_celltype("Macrophage MHCII High")
MHCII_DEG<-MHCII_DEG[order(MHCII_DEG$gene),]
peri_sig<-MHCII_DEG[which(MHCII_DEG$p_adjusted_periportal<0.005),]
peri_sig[order(peri_sig$gene),]
central_sig<-MHCII_DEG[which(MHCII_DEG$p_adjusted_pericentral<0.005),]
central_sig[order(central_sig$gene),]
MHCII_DEG[grep("C1QC|CXCL|APOE",MHCII_DEG$gene),]
inped_notadult(peri_sig)
inped_notadult(central_sig)
inadult_notped(peri_sig)
inadult_notped(central_sig)


diff_pedmap<-MHCII_DEG[which(MHCII_DEG$gene%in%ped_map_MHCII),]
central_sig_ped_map<-diff_pedmap[which(diff_pedmap$p_adjusted_pericentral<0.005),]
central_sig_ped_map[order(central_sig_ped_map$gene),]
diff_pedmap<-MHCII_DEG[which(MHCII_DEG$gene%in%ped_map_MHCII),]
peri_sig_ped_map<-diff_pedmap[which(diff_pedmap$p_adjusted_periportal<0.005),]
peri_sig_ped_map[order(peri_sig_ped_map$gene),]



Mono_mac_DEG<-DEG_celltype("Mono-Mac")
Mono_mac_DEG<-Mono_mac_DEG[order(Mono_mac_DEG$gene),]
Mono_mac_DEG[which(Mono_mac_DEG$p_adjusted_periportal<0.005),]
peri_sig<-Mono_mac_DEG[which(Mono_mac_DEG$p_adjusted_periportal<0.005),]
peri_sig[order(peri_sig$gene),]
central_sig<-Mono_mac_DEG[which(Mono_mac_DEG$p_adjusted_pericentral<0.005),]
central_sig[order(central_sig$gene),]
Mono_mac_DEG[grep("C1QC|CXCL|APOE",Mono_mac_DEG$gene),]
inped_notadult(peri_sig)
inped_notadult(central_sig)
inadult_notped(peri_sig)
inadult_notped(central_sig)

diff_pedmap<-Mono_mac_DEG[which(Mono_mac_DEG$gene%in%ped_map_RR),]
central_sig_ped_map<-diff_pedmap[which(diff_pedmap$p_adjusted_pericentral<0.005),]
central_sig_ped_map[order(central_sig_ped_map$gene),]
diff_pedmap<-Mono_mac_DEG[which(Mono_mac_DEG$gene%in%ped_map_RR),]
peri_sig_ped_map<-diff_pedmap[which(diff_pedmap$p_adjusted_periportal<0.005),]
peri_sig_ped_map[order(peri_sig_ped_map$gene),]


KC_DEG<-DEG_celltype("KC Like")
KC_DEG<-KC_DEG[order(KC_DEG$gene),]
peri_sig<-KC_DEG[which(KC_DEG$p_adjusted_periportal<0.005),]
peri_sig[order(peri_sig$gene),]
central_sig<-KC_DEG[which(KC_DEG$p_adjusted_pericentral<0.005),]
central_sig[order(central_sig$gene),]
KC_DEG[grep("C1QC|CXCL|APOE",KC_DEG$gene),]

diff_pedmap<-KC_DEG[which(KC_DEG$gene%in%ped_map_KC),]
central_sig_ped_map<-diff_pedmap[which(diff_pedmap$p_adjusted_pericentral<0.005),]
central_sig_ped_map[order(central_sig_ped_map$gene),]
diff_pedmap<-KC_DEG[which(KC_DEG$gene%in%ped_map_KC),]
peri_sig_ped_map<-diff_pedmap[which(diff_pedmap$p_adjusted_periportal<0.005),]
peri_sig_ped_map[order(peri_sig_ped_map$gene),]

inped_notadult(peri_sig)
inped_notadult(central_sig)
inadult_notped(peri_sig)
inadult_notped(central_sig)

KC_DEG[grep("CD14",KC_DEG$gene),]
peri_sig[grep("IGF",peri_sig$gene),]
central_sig[grep("IGF",central_sig$gene),]



HSCact_DEG<-DEG_celltype("HSC (Activated)")
HSCact_DEG<-HSCact_DEG[order(HSCact_DEG$gene),]
peri_sig<-HSCact_DEG[which(HSCact_DEG$p_adjusted_periportal<0.005),]
peri_sig[order(peri_sig$gene),]
central_sig<-HSCact_DEG[which(HSCact_DEG$p_adjusted_pericentral<0.005),]
central_sig[order(central_sig$gene),]
inped_notadult(peri_sig)
inped_notadult(central_sig)
inadult_notped(peri_sig)
inadult_notped(central_sig)

HSCque_DEG<-DEG_celltype("HSC (Quiescent)")
HSCque_DEG<-HSCque_DEG[order(HSCque_DEG$gene),]
peri_sig<-HSCque_DEG[which(HSCque_DEG$p_adjusted_periportal<0.005),]
peri_sig[order(peri_sig$gene),]
central_sig<-HSCque_DEG[which(HSCque_DEG$p_adjusted_pericentral<0.005),]
central_sig[order(central_sig$gene),]
inped_notadult(peri_sig)
inped_notadult(central_sig)
inadult_notped(peri_sig)
inadult_notped(central_sig)



custom_markers[which(custom_markers$Gene%in%peri_sig$gene),]
cell_markers_general[which(cell_markers_general$Gene%in%peri_sig$gene),]
cell_markers_general[which(cell_markers_general$Gene%in%Mono_mac_DEG[which(Mono_mac_DEG$p_adjusted_pericentral<0.005),]$gene),]


## shared periportal
intersect(C94_4_DEG[which(C94_4_DEG$p_adjusted_periportal<0.005),]$gene, C94_3_DEG[which(C94_3_DEG$p_adjusted_periportal<0.005),]$gene)
intersect(C85_DEG[which(C85_DEG$p_adjusted_periportal<0.005),]$gene, C105_DEG[which(C105_DEG$p_adjusted_periportal<0.005),]$gene)

## shared pericentral
intersect(C94_4_DEG[which(C94_4_DEG$p_adjusted_pericentral<0.005),]$gene, C94_3_DEG[which(C94_3_DEG$p_adjusted_pericentral<0.005),]$gene)
intersect(C85_DEG[which(C85_DEG$p_adjusted_pericentral<0.005),]$gene, C105_DEG[which(C105_DEG$p_adjusted_pericentral<0.005),]$gene)

ped_pericentral<-intersect(C85_DEG[which(C85_DEG$p_adjusted_pericentral<0.005),]$gene, C105_DEG[which(C105_DEG$p_adjusted_pericentral<0.005),]$gene)

cell_markers_general[which(cell_markers_general$Gene%in%ped_pericentral),]
custom_markers[which(custom_markers$Gene%in%ped_pericentral),]

cell_markers_general[which(cell_markers_general$Gene%in%C85_DEG[which(C85_DEG$p_adjusted_periportal<0.005),]$gene),]
custom_markers[which(custom_markers$Gene%in%C85_DEG[which(C85_DEG$p_adjusted_periportal<0.005),]$gene),]




########################
## plot gene
########################
plot_gene_expression("C1QC","KC Like")
plot_gene_expression("APOE","KC Like")
plot_gene_expression("CD14","KC Like")

plot_gene_expression("FCGR3A","Mono-Mac")


########################
## plot gene expression over distance
########################
distance_gene_expression("CD14","Portal","KC Like")
distance_gene_expression("CD14","Central","KC Like")

distance_gene_expression("APOE","Portal","KC Like")
distance_gene_expression("APOE","Central","KC Like")
distance_gene_expression("C1QC","Portal","KC Like")
distance_gene_expression("C1QC","Central","KC Like")


distance_gene_expression("APOE","Portal","KC Like")
distance_gene_expression("CD5L","Portal","KC Like")
distance_gene_expression("MARCO","Portal","KC Like")






######################
## Significance all sample combined
######################

## Excluding C94_2

d10x.list <- sapply(2:7, function(x){
  print(x)
  counts<-read.csv(count_files[x])
  print(samples[x])
  counts$X<-NULL
  rownames(counts) <- counts$cell_id
  counts$cell_id <- NULL
  
  # Transpose the data and convert to sparse matrix.
  mat <- as(t(as.matrix(counts)), "sparseMatrix")
  
  seu <- CreateSeuratObject(counts=mat)
  seu$sample<-strsplit(strsplit(count_files[x],"/")[[1]][6],"__")[[1]][3]
  seu
  
  ## zonation distance
  
  plt_zonation_xenium<-zonation_scores[[x]]
  
  plt_zonation_xenium<-plt_zonation_xenium[match(colnames(seu), plt_zonation_xenium$cell),]
  identical(colnames(seu), plt_zonation_xenium$cell)
  rownames(plt_zonation_xenium)<-plt_zonation_xenium$cell
  
  seu<-AddMetaData(seu, plt_zonation_xenium)
  seu})

d10x.list

xenium.obj <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "xenium_liver")
rm(d10x.list)
gc()


xenium.obj<-JoinLayers(xenium.obj)
xenium.obj <- NormalizeData(xenium.obj)


#######
## Differential Expression with portal distance
#######

celltypes<-names(table(xenium.obj$CellType))[which(table(xenium.obj$CellType)>50)]

pval_distance_allsamples<-lapply(1:length(celltypes), function(y){
  celltype<-celltypes[y]
  print(celltype)
  
  ## data looks pretty zero inflated so will fit a negative binomial
  ## normalize
  seu_celltype<-subset(xenium.obj, subset = CellType == celltype)
  gene_exp<-as.matrix(GetAssayData(object = seu_celltype, assay = "RNA", slot = "data"))
  
  meta<-seu_celltype@meta.data
  
  identical(rownames(meta), colnames(gene_exp))
  
  ## lm negative binomial on all genes
  pval_distance<-do.call(rbind, lapply(1:nrow(gene_exp), function(x) {
    gene<-rownames(gene_exp)[x]
    #print(gene)
    meta$exp<-gene_exp[which(rownames(gene_exp)==gene),]
    
    if((sum(meta$exp==0)/nrow(meta))<0.85){ #more that 15% of cell express
      # Negative binomial GLMM using the function glmer.nb()
      
      testFunction_portal <- function (meta) {return(tryCatch(glm.nb(exp ~ distance_to_periportal, data = meta), error=function(e) NULL)) }  
      testFunction_central <- function (meta) {return(tryCatch(glm.nb(exp ~ distance_to_pericentral, data = meta), error=function(e) NULL)) }  
      
      if(is.null(testFunction_portal(meta))&is.null(testFunction_central(meta))){
        data.frame(gene = gene,pval_periportal = "Errored",estimate_periportal = "Errored",pval_pericentral = "Errored",estimate_pericentral = "Errored")}else{
          if(is.null(testFunction_portal(meta))){
            mnb_pericentral <- glm.nb(exp ~ distance_to_pericentral, data = meta)
            data.frame(gene = gene,pval_periportal = "Errored",estimate_periportal = "Errored",
                       pval_pericentral = summary(mnb_pericentral)$coefficients["distance_to_pericentral","Pr(>|z|)"],
                       estimate_pericentral = summary(mnb_pericentral)$coefficients["distance_to_pericentral","Estimate"])}else{
                         if(is.null(testFunction_central(meta))){
                           mnb_periportal <- glm.nb(exp ~ distance_to_periportal, data = meta)
                           data.frame(gene = gene,
                                      pval_periportal = summary(mnb_periportal)$coefficients["distance_to_periportal","Pr(>|z|)"],
                                      estimate_periportal = summary(mnb_periportal)$coefficients["distance_to_periportal","Estimate"],
                                      pval_pericentral = "Errored",estimate_pericentral = "Errored")}else{
                                        mnb_periportal <- glm.nb(exp ~ distance_to_periportal, data = meta)
                                        mnb_pericentral <- glm.nb(exp ~ distance_to_pericentral, data = meta)
                                        data.frame(gene = gene,
                                                   pval_periportal = summary(mnb_periportal)$coefficients["distance_to_periportal","Pr(>|z|)"],
                                                   estimate_periportal = summary(mnb_periportal)$coefficients["distance_to_periportal","Estimate"],
                                                   pval_pericentral = summary(mnb_pericentral)$coefficients["distance_to_pericentral","Pr(>|z|)"],
                                                   estimate_pericentral = summary(mnb_pericentral)$coefficients["distance_to_pericentral","Estimate"])}}}
      
      
    }
  }))
  
  pval_distance$p_adjusted_periportal<-p.adjust(pval_distance$pval_periportal, method = "bonferroni")
  pval_distance$p_adjusted_pericentral<-p.adjust(pval_distance$pval_pericentral, method = "bonferroni")
  pval_distance$celltype<-celltype
  pval_distance
})


pval_distance_allsamples<-do.call(rbind, pval_distance_allsamples)

save(pval_distance_allsamples, file=here("data/pval_zonation_allsamples_combined.RData"))

load(here("data/pval_zonation_allsamples_combined.RData"))

pval_distance_allsamples_myeloid<-pval_distance_allsamples[which(pval_distance_allsamples$celltype%in%c("Macrophage MHCII High","KC Like","Mono-Mac")),]
pval_distance_allsamples_myeloid<-pval_distance_allsamples_myeloid[which(pval_distance_allsamples_myeloid$p_adjusted_periportal<0.05 | pval_distance_allsamples_myeloid$p_adjusted_pericentral<0.05),]
write.csv(pval_distance_allsamples_myeloid, file=here("data/sig_myeloid_zonated.csv"), quote=F)

####################
## ped differential genes
####################
cell_markers_general<-read.csv(here("data/Xenium_CombinedPanel.csv"))
custom_markers<-read.csv(here("data/Xenium_MacParlandGeneListUpdated.csv"))
custom_markers$Gene[grep("HAL", custom_markers$Gene)]<-"HAL"

xenium<-read.csv(file=here("/home/redgar/Documents/xenium_liver/data/Xenium_CombinedPanel.csv"))


sig_de_age_RR<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_RR.csv"))
sig_de_age_KC<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_KC.csv"))
sig_de_age_MHCII<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_MHCII.csv"))

sig_de_age_KC[which(sig_de_age_KC$X%in%xenium$Gene),]
sig_de_age_RR[which(sig_de_age_RR$X%in%xenium$Gene),]
sig_de_age_MHCII[which(sig_de_age_MHCII$X%in%xenium$Gene),]

intersect(xenium$Gene, sig_de_age_KC$X)
intersect(xenium$Gene, sig_de_age_RR$X)
intersect(xenium$Gene, sig_de_age_MHCII$X)

ped_map_KC<-intersect(xenium$Gene, sig_de_age_KC$X)
ped_map_RR<-intersect(xenium$Gene, sig_de_age_RR$X)
ped_map_MHCII<-intersect(xenium$Gene, sig_de_age_MHCII$X)


KC_zonated<-pval_distance_allsamples[which(pval_distance_allsamples$celltype=="KC Like"),]
KC_periportal<-KC_zonated[which(KC_zonated$p_adjusted_periportal<0.05),]
KC_pericentral<-KC_zonated[which(KC_zonated$p_adjusted_pericentral<0.05),]

KC_periportal[which(KC_periportal$gene%in%ped_map_KC),]
KC_pericentral[which(KC_pericentral$gene%in%ped_map_KC),]

length(unique(KC_periportal[which(KC_periportal$gene%in%ped_map_KC),]$gene, KC_pericentral[which(KC_pericentral$gene%in%ped_map_KC),]$gene))#8
length(unique(KC_periportal$gene, KC_pericentral$gene)) #76



distance_gene_expression_combined("APOE","Portal","KC Like")
distance_gene_expression_combined("IL10","Portal","KC Like")


distance_gene_expression_combined("APOE","Central","KC Like")
distance_gene_expression_combined("IL10","Central","KC Like")

distance_gene_expression_combined("MARCO","Portal","KC Like")
distance_gene_expression_combined("CD5L","Portal","KC Like")



distance_gene_expression_combined("APOE","Central",NA)
distance_gene_expression_combined("IL10","Central",NA)

distance_gene_expression_combined("APOE","Portal",NA)
distance_gene_expression_combined("IL10","Portal",NA)






gene_myeloid_portal<-function(gene){
  celltype<-c("Macrophage MHCII High","KC Like","Mono-Mac")
  count_files<-c("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_2/cell_gene_matrices/2024_03_21_16_22_10/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_3/cell_gene_matrices/2024_03_21_16_38_48/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_4/cell_gene_matrices/2024_03_21_17_03_05/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C105/cell_gene_matrices/2024_05_17_13_32_03/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C85/cell_gene_matrices/2024_05_17_13_29_40/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C95/cell_gene_matrices/2024_07_04_12_35_29/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C101/cell_gene_matrices/2024_06_28_09_23_40/expr_mat.csv")
  samples<-c("C94_2", "C94_3","C94_4","C105","C85","C95","C101")
  models<-c("2024_03_21_16_22_10","2024_03_21_16_38_48","2024_03_21_17_03_05","2024_05_17_13_32_03","2024_05_17_13_29_40","2024_07_04_12_35_29","2024_06_28_09_23_40")
  
  load(here("data","zonation_allcells.RData"))
  zonation_scores<-do.call(rbind, zonation_scores)
  
  meta_exp<-lapply(2:length(samples), function(x){
    data_dir<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x],"/cell_gene_matrices/", models[x], "/expr_mat.csv", sep="")
    tiff_path<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x], "/model_outputs/",models[x],"/test_output/epoch_1_step_4000_connected.tif", sep="")
    smpl<-samples[x]
    print(smpl)
    
    ## gene expression
    counts<-read.csv(count_files[x])
    print(samples[x])
    counts$X<-NULL
    rownames(counts) <- counts$cell_id
    counts$cell_id <- NULL
    
    # Transpose the data and convert to sparse matrix.
    mat <- as(t(as.matrix(counts)), "sparseMatrix")
    
    seu <- CreateSeuratObject(counts=mat)
    seu$sample<-strsplit(strsplit(count_files[x],"/")[[1]][6],"__")[[1]][3]
    seu
    
    zonation_scores_sample<-zonation_scores[which(zonation_scores$sample==smpl),]
    zonation_scores_sample<-zonation_scores_sample[match(colnames(seu), zonation_scores_sample$cell),]
    identical(colnames(seu), zonation_scores_sample$cell)
    rownames(zonation_scores_sample)<-zonation_scores_sample$cell
    
    seu<-AddMetaData(seu, zonation_scores_sample)
    seu <- NormalizeData(seu)
    
    gene_exp<-as.matrix(GetAssayData(object = seu, assay = "RNA", slot = "data"))
    
    meta<-seu@meta.data
    meta$exp<-gene_exp[which(rownames(gene_exp)==gene),]
    
    meta})
  
  meta_exp<-do.call(rbind, meta_exp)
  
  
  # distance_plt<- ggplot()+
  #   geom_line(aes(distance_to_periportal, exp, color=CellType, group=interaction(CellType, sample)),
  #             meta_exp[which((meta_exp$CellType%in%celltype)),],
  #             stat="smooth",method = "glm.nb", 
  #             size = 0.5,
  #             alpha = 0.5)+
  #   stat_smooth(aes(distance_to_periportal, exp, color=CellType),meta_exp[which((meta_exp$CellType%in%celltype)),] ,method="glm.nb", se=T,fill = "grey90", linewidth=1)+
  #   theme_bw()+
  #   theme(strip.background = element_rect(fill="white"))+
  #   xlab("Periportal Distance")+ylab(paste(gene, "Expression"))+
  #   colscale_cellType+theme(strip.text.y = element_text(angle = 0))
  # distance_plt
  # 
  
  # distance_plt<- ggplot(meta_exp[which((meta_exp$CellType%in%celltype)),], aes(distance_to_periportal, exp, color=CellType))+
  #   #stat_smooth(aes(distance_to_periportal, exp, color=CellType, group=interaction(CellType, sample)),meta_exp[which((meta_exp$CellType%in%celltype)),] ,method="glm.nb", se=F,linewidth=0.5, alpha=0.1)+
  #   geom_line(aes(distance_to_periportal, exp, color=CellType, group=interaction(CellType, sample)),
  #             meta_exp[which((meta_exp$CellType%in%celltype)),],
  #             stat="smooth",method = "glm.nb", 
  #             size = 0.5,
  #             alpha = 0.4)+ylim(0, max(meta_exp$exp))+
  #   stat_smooth(method="glm.nb", se=T,fill = "grey90", linewidth=1.5)+
  #   theme_bw()+
  #   theme(strip.background = element_rect(fill="white"))+
  #   xlab("Periportal Distance")+ylab(paste(gene, "Expression"))+
  #   colscale_cellType+theme(strip.text.y = element_text(angle = 0))
  # distance_plt
  # 
  meta_exp$AgeGroup<-as.factor(meta_exp$sample)
  levels(meta_exp$AgeGroup)<-c( "Adult" ,"Adult", "Adult", "Pediatric" , "Pediatric" ,  "Adult",   "Adult" )
  
  distance_plt<- ggplot(meta_exp[which((meta_exp$CellType%in%celltype)),], aes(distance_to_periportal, exp, color=CellType))+
    #stat_smooth(aes(distance_to_periportal, exp, color=CellType, group=interaction(CellType, sample)),meta_exp[which((meta_exp$CellType%in%celltype)),] ,method="glm.nb", se=F,linewidth=0.5, alpha=0.1)+
    geom_line(aes(distance_to_periportal, exp, color=CellType, group =interaction(CellType, AgeGroup), linetype=AgeGroup),
              meta_exp[which((meta_exp$CellType%in%celltype)),],
              stat="smooth",method = "glm.nb", 
              size = 0.75,
              alpha = 0.5)+
    stat_smooth(method="glm.nb", se=T,fill = "grey90", linewidth=1.5)+
    theme_bw()+
    theme(strip.background = element_rect(fill="white"))+
    xlab("Periportal Distance")+ylab(paste(gene, "Expression"))+
    colscale_cellType+theme(strip.text.y = element_text(angle = 0))+
    scale_linetype_manual(values = c("11", "dotted"))
  distance_plt
  
}



gene_myeloid_portal("MARCO")
save_plts(gene_myeloid_portal("MARCO"),"MARCO_periportal_myeloid" ,  w=6, h=4)

gene_myeloid_portal("APOA1")
save_plts(gene_myeloid_portal("APOA1"),"APOA1_periportal_myeloid" ,  w=6, h=4)
    
gene_myeloid_portal("APOE")
save_plts(gene_myeloid_portal("APOE"),"APOE_periportal_myeloid" ,  w=6, h=4)





### cd5l limit weird
#gene_myeloid_portal("CD5L")

load(here("data","zonation_allcells.RData"))
zonation_scores<-do.call(rbind, zonation_scores)
gene<-"CD5L"
celltype<-c("Macrophage MHCII High","KC Like","Mono-Mac")

meta_exp<-lapply(2:length(samples), function(x){
  data_dir<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x],"/cell_gene_matrices/", models[x], "/expr_mat.csv", sep="")
  tiff_path<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x], "/model_outputs/",models[x],"/test_output/epoch_1_step_4000_connected.tif", sep="")
  smpl<-samples[x]
  print(smpl)
  
  ## gene expression
  counts<-read.csv(count_files[x])
  print(samples[x])
  counts$X<-NULL
  rownames(counts) <- counts$cell_id
  counts$cell_id <- NULL
  
  # Transpose the data and convert to sparse matrix.
  mat <- as(t(as.matrix(counts)), "sparseMatrix")
  
  seu <- CreateSeuratObject(counts=mat)
  seu$sample<-strsplit(strsplit(count_files[x],"/")[[1]][6],"__")[[1]][3]
  seu
  
  zonation_scores_sample<-zonation_scores[which(zonation_scores$sample==smpl),]
  zonation_scores_sample<-zonation_scores_sample[match(colnames(seu), zonation_scores_sample$cell),]
  identical(colnames(seu), zonation_scores_sample$cell)
  rownames(zonation_scores_sample)<-zonation_scores_sample$cell
  
  seu<-AddMetaData(seu, zonation_scores_sample)
  seu <- NormalizeData(seu)
  
  gene_exp<-as.matrix(GetAssayData(object = seu, assay = "RNA", slot = "data"))
  
  meta<-seu@meta.data
  meta$exp<-gene_exp[which(rownames(gene_exp)==gene),]
  
  meta})

meta_exp<-do.call(rbind, meta_exp)

meta_exp$AgeGroup<-as.factor(meta_exp$sample)
levels(meta_exp$AgeGroup)<-c( "Adult" ,"Adult", "Adult", "Pediatric" , "Pediatric" ,  "Adult",   "Adult" )

distance_plt<- ggplot(meta_exp[which((meta_exp$CellType%in%celltype)),], aes(distance_to_periportal, exp, color=CellType))+
  #stat_smooth(aes(distance_to_periportal, exp, color=CellType, group=interaction(CellType, sample)),meta_exp[which((meta_exp$CellType%in%celltype)),] ,method="glm.nb", se=F,linewidth=0.5, alpha=0.1)+
  geom_line(aes(distance_to_periportal, exp, color=CellType, group =interaction(CellType, AgeGroup), linetype=AgeGroup),
            meta_exp[which((meta_exp$CellType%in%celltype)),],
            stat="smooth",method = "glm.nb", 
            size = 0.75,
            alpha = 0.5)+ylim(0, 6)+
  stat_smooth(method="glm.nb", se=T,fill = "grey90", linewidth=1.5)+
  theme_bw()+
  theme(strip.background = element_rect(fill="white"))+
  xlab("Periportal Distance")+ylab(paste(gene, "Expression"))+
  colscale_cellType+theme(strip.text.y = element_text(angle = 0))+
  scale_linetype_manual(values = c("11", "dotted"))
distance_plt

save_plts(distance_plt,"CD5L_periportal_myeloid" ,  w=6, h=4)


######################
## Significance ped and adult seperate
######################

## Excluding C94_2

d10x.list.adult <- sapply(c(2,3,6,7), function(x){
  print(x)
  counts<-read.csv(count_files[x])
  print(samples[x])
  counts$X<-NULL
  rownames(counts) <- counts$cell_id
  counts$cell_id <- NULL
  
  # Transpose the data and convert to sparse matrix.
  mat <- as(t(as.matrix(counts)), "sparseMatrix")
  
  seu <- CreateSeuratObject(counts=mat)
  seu$sample<-strsplit(strsplit(count_files[x],"/")[[1]][6],"__")[[1]][3]
  seu
  
  ## zonation distance
  
  plt_zonation_xenium<-zonation_scores[[x]]
  
  plt_zonation_xenium<-plt_zonation_xenium[match(colnames(seu), plt_zonation_xenium$cell),]
  identical(colnames(seu), plt_zonation_xenium$cell)
  rownames(plt_zonation_xenium)<-plt_zonation_xenium$cell
  
  seu<-AddMetaData(seu, plt_zonation_xenium)
  seu})

d10x.list.adult

xenium.obj <- merge(d10x.list.adult[[1]], y= d10x.list.adult[2:length(d10x.list.adult)], merge.data=TRUE, project = "xenium_liver")
rm(d10x.list.adult)
gc()


xenium.obj<-JoinLayers(xenium.obj)
xenium.obj <- NormalizeData(xenium.obj)


#######
## Differential Expression with portal distance
#######

celltypes<-names(table(xenium.obj$CellType))[which(table(xenium.obj$CellType)>50)]



##### adults
pval_distance_allsamples<-lapply(1:length(celltypes), function(y){
  celltype<-celltypes[y]
  print(celltype)
  
  ## data looks pretty zero inflated so will fit a negative binomial
  ## normalize
  seu_celltype<-subset(xenium.obj, subset = CellType == celltype)
  gene_exp<-as.matrix(GetAssayData(object = seu_celltype, assay = "RNA", slot = "data"))
  
  meta<-seu_celltype@meta.data
  
  identical(rownames(meta), colnames(gene_exp))
  
  ## lm negative binomial on all genes
  pval_distance<-do.call(rbind, lapply(1:nrow(gene_exp), function(x) {
    gene<-rownames(gene_exp)[x]
    #print(gene)
    meta$exp<-gene_exp[which(rownames(gene_exp)==gene),]
    
    if((sum(meta$exp==0)/nrow(meta))<0.85){ #more that 15% of cell express
      # Negative binomial GLMM using the function glmer.nb()
      
      testFunction_portal <- function (meta) {return(tryCatch(glm.nb(exp ~ distance_to_periportal, data = meta), error=function(e) NULL)) }  
      testFunction_central <- function (meta) {return(tryCatch(glm.nb(exp ~ distance_to_pericentral, data = meta), error=function(e) NULL)) }  
      
      if(is.null(testFunction_portal(meta))&is.null(testFunction_central(meta))){
        data.frame(gene = gene,pval_periportal = "Errored",estimate_periportal = "Errored",pval_pericentral = "Errored",estimate_pericentral = "Errored")}else{
          if(is.null(testFunction_portal(meta))){
            mnb_pericentral <- glm.nb(exp ~ distance_to_pericentral, data = meta)
            data.frame(gene = gene,pval_periportal = "Errored",estimate_periportal = "Errored",
                       pval_pericentral = summary(mnb_pericentral)$coefficients["distance_to_pericentral","Pr(>|z|)"],
                       estimate_pericentral = summary(mnb_pericentral)$coefficients["distance_to_pericentral","Estimate"])}else{
                         if(is.null(testFunction_central(meta))){
                           mnb_periportal <- glm.nb(exp ~ distance_to_periportal, data = meta)
                           data.frame(gene = gene,
                                      pval_periportal = summary(mnb_periportal)$coefficients["distance_to_periportal","Pr(>|z|)"],
                                      estimate_periportal = summary(mnb_periportal)$coefficients["distance_to_periportal","Estimate"],
                                      pval_pericentral = "Errored",estimate_pericentral = "Errored")}else{
                                        mnb_periportal <- glm.nb(exp ~ distance_to_periportal, data = meta)
                                        mnb_pericentral <- glm.nb(exp ~ distance_to_pericentral, data = meta)
                                        data.frame(gene = gene,
                                                   pval_periportal = summary(mnb_periportal)$coefficients["distance_to_periportal","Pr(>|z|)"],
                                                   estimate_periportal = summary(mnb_periportal)$coefficients["distance_to_periportal","Estimate"],
                                                   pval_pericentral = summary(mnb_pericentral)$coefficients["distance_to_pericentral","Pr(>|z|)"],
                                                   estimate_pericentral = summary(mnb_pericentral)$coefficients["distance_to_pericentral","Estimate"])}}}
      
      
    }
  }))
  
  pval_distance$p_adjusted_periportal<-p.adjust(pval_distance$pval_periportal, method = "bonferroni")
  pval_distance$p_adjusted_pericentral<-p.adjust(pval_distance$pval_pericentral, method = "bonferroni")
  pval_distance$celltype<-celltype
  pval_distance
})


pval_distance_adults<-do.call(rbind, pval_distance_allsamples)

save(pval_distance_adults, file=here("data/pval_zonation_adults_combined.RData"))


####################
## ped differential genes
####################
cell_markers_general<-read.csv(here("data/Xenium_CombinedPanel.csv"))
custom_markers<-read.csv(here("data/Xenium_MacParlandGeneListUpdated.csv"))
custom_markers$Gene[grep("HAL", custom_markers$Gene)]<-"HAL"

xenium<-read.csv(file=here("/home/redgar/Documents/xenium_liver/data/Xenium_CombinedPanel.csv"))


sig_de_age_RR<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_RR.csv"))
sig_de_age_KC<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_KC.csv"))
sig_de_age_MHCII<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_MHCII.csv"))

sig_de_age_KC[which(sig_de_age_KC$X%in%xenium$Gene),]
sig_de_age_RR[which(sig_de_age_RR$X%in%xenium$Gene),]
sig_de_age_MHCII[which(sig_de_age_MHCII$X%in%xenium$Gene),]

intersect(xenium$Gene, sig_de_age_KC$X)
intersect(xenium$Gene, sig_de_age_RR$X)
intersect(xenium$Gene, sig_de_age_MHCII$X)

ped_map_KC<-intersect(xenium$Gene, sig_de_age_KC$X)
ped_map_RR<-intersect(xenium$Gene, sig_de_age_RR$X)
ped_map_MHCII<-intersect(xenium$Gene, sig_de_age_MHCII$X)


KC_zonated<-pval_distance_adults[which(pval_distance_adults$celltype=="KC Like"),]
KC_periportal<-KC_zonated[which(KC_zonated$p_adjusted_periportal<0.05),]
KC_pericentral<-KC_zonated[which(KC_zonated$p_adjusted_pericentral<0.05),]

KC_periportal[which(KC_periportal$gene%in%ped_map_KC),]
KC_pericentral[which(KC_pericentral$gene%in%ped_map_KC),]



distance_gene_expression_combined("APOE","Portal","KC Like")
distance_gene_expression_combined("MARCO","Portal","KC Like")
distance_gene_expression_combined("CD5L","Portal","KC Like")




######### peds
d10x.list.peds <- sapply(c(4,5), function(x){
  print(x)
  counts<-read.csv(count_files[x])
  print(samples[x])
  counts$X<-NULL
  rownames(counts) <- counts$cell_id
  counts$cell_id <- NULL
  
  # Transpose the data and convert to sparse matrix.
  mat <- as(t(as.matrix(counts)), "sparseMatrix")
  
  seu <- CreateSeuratObject(counts=mat)
  seu$sample<-strsplit(strsplit(count_files[x],"/")[[1]][6],"__")[[1]][3]
  seu
  
  ## zonation distance
  
  plt_zonation_xenium<-zonation_scores[[x]]
  
  plt_zonation_xenium<-plt_zonation_xenium[match(colnames(seu), plt_zonation_xenium$cell),]
  identical(colnames(seu), plt_zonation_xenium$cell)
  rownames(plt_zonation_xenium)<-plt_zonation_xenium$cell
  
  seu<-AddMetaData(seu, plt_zonation_xenium)
  seu})

d10x.list.peds

xenium.obj <- merge(d10x.list.peds[[1]], y= d10x.list.peds[2:length(d10x.list.peds)], merge.data=TRUE, project = "xenium_liver")
rm(d10x.list.peds)
gc()


xenium.obj<-JoinLayers(xenium.obj)
xenium.obj <- NormalizeData(xenium.obj)



pval_distance_peds<-lapply(1:length(celltypes), function(y){
  celltype<-celltypes[y]
  print(celltype)
  
  ## data looks pretty zero inflated so will fit a negative binomial
  ## normalize
  seu_celltype<-subset(xenium.obj, subset = CellType == celltype)
  gene_exp<-as.matrix(GetAssayData(object = seu_celltype, assay = "RNA", slot = "data"))
  
  meta<-seu_celltype@meta.data
  
  identical(rownames(meta), colnames(gene_exp))
  
  ## lm negative binomial on all genes
  pval_distance<-do.call(rbind, lapply(1:nrow(gene_exp), function(x) {
    gene<-rownames(gene_exp)[x]
    #print(gene)
    meta$exp<-gene_exp[which(rownames(gene_exp)==gene),]
    
    if((sum(meta$exp==0)/nrow(meta))<0.85){ #more that 15% of cell express
      # Negative binomial GLMM using the function glmer.nb()
      
      testFunction_portal <- function (meta) {return(tryCatch(glm.nb(exp ~ distance_to_periportal, data = meta), error=function(e) NULL)) }  
      testFunction_central <- function (meta) {return(tryCatch(glm.nb(exp ~ distance_to_pericentral, data = meta), error=function(e) NULL)) }  
      
      if(is.null(testFunction_portal(meta))&is.null(testFunction_central(meta))){
        data.frame(gene = gene,pval_periportal = "Errored",estimate_periportal = "Errored",pval_pericentral = "Errored",estimate_pericentral = "Errored")}else{
          if(is.null(testFunction_portal(meta))){
            mnb_pericentral <- glm.nb(exp ~ distance_to_pericentral, data = meta)
            data.frame(gene = gene,pval_periportal = "Errored",estimate_periportal = "Errored",
                       pval_pericentral = summary(mnb_pericentral)$coefficients["distance_to_pericentral","Pr(>|z|)"],
                       estimate_pericentral = summary(mnb_pericentral)$coefficients["distance_to_pericentral","Estimate"])}else{
                         if(is.null(testFunction_central(meta))){
                           mnb_periportal <- glm.nb(exp ~ distance_to_periportal, data = meta)
                           data.frame(gene = gene,
                                      pval_periportal = summary(mnb_periportal)$coefficients["distance_to_periportal","Pr(>|z|)"],
                                      estimate_periportal = summary(mnb_periportal)$coefficients["distance_to_periportal","Estimate"],
                                      pval_pericentral = "Errored",estimate_pericentral = "Errored")}else{
                                        mnb_periportal <- glm.nb(exp ~ distance_to_periportal, data = meta)
                                        mnb_pericentral <- glm.nb(exp ~ distance_to_pericentral, data = meta)
                                        data.frame(gene = gene,
                                                   pval_periportal = summary(mnb_periportal)$coefficients["distance_to_periportal","Pr(>|z|)"],
                                                   estimate_periportal = summary(mnb_periportal)$coefficients["distance_to_periportal","Estimate"],
                                                   pval_pericentral = summary(mnb_pericentral)$coefficients["distance_to_pericentral","Pr(>|z|)"],
                                                   estimate_pericentral = summary(mnb_pericentral)$coefficients["distance_to_pericentral","Estimate"])}}}
      
      
    }
  }))
  
  pval_distance$p_adjusted_periportal<-p.adjust(pval_distance$pval_periportal, method = "bonferroni")
  pval_distance$p_adjusted_pericentral<-p.adjust(pval_distance$pval_pericentral, method = "bonferroni")
  pval_distance$celltype<-celltype
  pval_distance
})


pval_distance_peds<-do.call(rbind, pval_distance_peds)

save(pval_distance_peds, file=here("data/pval_zonation_peds_combined.RData"))


## ped differential genes
cell_markers_general<-read.csv(here("data/Xenium_CombinedPanel.csv"))
custom_markers<-read.csv(here("data/Xenium_MacParlandGeneListUpdated.csv"))
custom_markers$Gene[grep("HAL", custom_markers$Gene)]<-"HAL"

xenium<-read.csv(file=here("/home/redgar/Documents/xenium_liver/data/Xenium_CombinedPanel.csv"))


sig_de_age_RR<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_RR.csv"))
sig_de_age_KC<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_KC.csv"))
sig_de_age_MHCII<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_MHCII.csv"))

sig_de_age_KC[which(sig_de_age_KC$X%in%xenium$Gene),]
sig_de_age_RR[which(sig_de_age_RR$X%in%xenium$Gene),]
sig_de_age_MHCII[which(sig_de_age_MHCII$X%in%xenium$Gene),]

intersect(xenium$Gene, sig_de_age_KC$X)
intersect(xenium$Gene, sig_de_age_RR$X)
intersect(xenium$Gene, sig_de_age_MHCII$X)

ped_map_KC<-intersect(xenium$Gene, sig_de_age_KC$X)
ped_map_RR<-intersect(xenium$Gene, sig_de_age_RR$X)
ped_map_MHCII<-intersect(xenium$Gene, sig_de_age_MHCII$X)


KC_zonated<-pval_distance_peds[which(pval_distance_peds$celltype=="KC Like"),]
KC_periportal<-KC_zonated[which(KC_zonated$p_adjusted_periportal<0.05),]
KC_pericentral<-KC_zonated[which(KC_zonated$p_adjusted_pericentral<0.05),]

KC_periportal[which(KC_periportal$gene%in%ped_map_KC),]
KC_pericentral[which(KC_pericentral$gene%in%ped_map_KC),]




#######
## zonation of gene in hepatocytes
#######

celltypes<-c("Hepatocyte (Periportal)", "Hepatocyte (Pericentral)", "Hepatocyte (Cycling)")

pval_distance_allsamples_hepatocytes<-lapply(1:length(celltypes), function(y){
  celltype<-celltypes[y]
  print(celltype)
  
  ## data looks pretty zero inflated so will fit a negative binomial
  ## normalize
  seu_celltype<-subset(xenium.obj, subset = CellType == celltype)
  gene_exp<-as.matrix(GetAssayData(object = seu_celltype, assay = "RNA", slot = "data"))
  
  meta<-seu_celltype@meta.data
  
  identical(rownames(meta), colnames(gene_exp))
  
  ## lm negative binomial on all genes
  pval_distance<-do.call(rbind, lapply(1:nrow(gene_exp), function(x) {
    gene<-rownames(gene_exp)[x]
    #print(gene)
    meta$exp<-gene_exp[which(rownames(gene_exp)==gene),]
    
    if((sum(meta$exp==0)/nrow(meta))<0.85){ #more that 15% of cell express
      # Negative binomial GLMM using the function glmer.nb()
      
      testFunction_zonation_zscore <- function (meta) {return(tryCatch(glm.nb(exp ~ zonation_zscore, data = meta), error=function(e) NULL)) }  

      if(is.null(testFunction_zonation_zscore(meta))){
        data.frame(gene = gene,pval_periportal = "Errored",estimate_periportal = "Errored",pval_pericentral = "Errored",estimate_pericentral = "Errored")}else{
          if(is.null(testFunction_portal(meta))){
            mnb_zonation_zscore <- glm.nb(exp ~ zonation_zscore, data = meta)
            data.frame(gene = gene,
                       pval_zonation_zscore = summary(mnb_zonation_zscore)$coefficients["zonation_zscore","Pr(>|z|)"],
                       estimate_zonation_zscore = summary(mnb_zonation_zscore)$coefficients["zonation_zscore","Estimate"])}}
          }}))
  
  pval_distance$p_adjusted_zonation_zscore<-p.adjust(pval_distance$pval_zonation_zscore, method = "bonferroni")
  pval_distance$celltype<-celltype
  pval_distance
})

pval_distance_allsamples_hepatocytes<-do.call(rbind, pval_distance_allsamples_hepatocytes)


save(pval_distance_allsamples_hepatocytes, file=here("data/pval_zonation_hepatocytes_allsamples_combined.RData"))

