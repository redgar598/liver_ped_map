#'---
#'title: scRNAseq Differential Expression
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---



#'### Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(here)
library(ggplot2)
library(reshape2)
library(gridExtra)
#library(limma)
library(cowplot)
library(gtools)
#library(ggsignif)


options(stringsAsFactors = FALSE)

source("R_functions/pretty_plots.R")



## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","d10x_adult_ped_raw.rds"))


######
## add cell type labels
######
load(here("data","adult_ped_cellRefined.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)




##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")


## testing factor
d10x$cell_sex<-paste(as.character(d10x$CellType_refined), d10x$Sex, sep = "_")
Idents(d10x) <- "cell_sex"

table(d10x$CellType_refined, d10x$Sex)


#MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene,
#consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
#and linear regression for the continuous process (i.e., the expression level). 

cell_types<-unique(as.character(d10x$CellType_refined))
cell_types<-cell_types[-grep("Hepatocyte Like",cell_types)]
cell_types[grep("CD3",cell_types)]<-"CD3"
cell_types[grep("DEFA",cell_types)]<-"DEFA"

contrasts_celltype_age<-do.call(rbind,lapply(1:length(cell_types), function(x){
  combinations(n = 2, r = 2, v = d10x$cell_sex[grep(cell_types[x],d10x$cell_sex)], repeats.allowed = FALSE)
}))

contrasts_celltype_age

nrow(contrasts_celltype_age)


diff_exp_all_sex<-lapply(1:nrow(contrasts_celltype_age), function(x){
  de<-FindMarkers(d10x, ident.1 = contrasts_celltype_age[x,1], ident.2 = contrasts_celltype_age[x,2], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
  print(paste(contrasts_celltype_age[x,1],"vs", contrasts_celltype_age[x,2],":", nrow(de), sep=" "))
  de$gene<-rownames(de)
  rownames(de)<-NULL
  de<-de[,c(6,1:5)]
  de$cell.1<-contrasts_celltype_age[x,1]
  de$cell.2<-contrasts_celltype_age[x,2]
  de})


diff_exp_all_sex<-do.call(rbind, diff_exp_all_sex)

save(diff_exp_all_sex, file=here("data","adult_ped_sex_diff_genes.RData"))
#load(file=here("data","adult_ped_sex_diff_genes.RData"))



##################
## IFN signalling
##################
#Lists from GSEA "HALLMARK_INTERFERON_GAMMA_RESPONSE" and "HALLMARK_INTERFERON_ALPHA_RESPONSE"
IFNa<-c("ADAR","B2M","BATF2","BST2","C1S","CASP1","CASP8","CCRL2","CD47","CD74","CMPK2","CNP","CSF1","CXCL10","CXCL11","DDX60","DHX58","EIF2AK2","ELF1","EPSTI1","MVB12A","TENT5A","CMTR1","GBP2","GBP4","GMPR","HERC6","HLA-C","IFI27","IFI30","IFI35","IFI44","IFI44L","IFIH1","IFIT2","IFIT3","IFITM1","IFITM2","IFITM3","IL15","IL4R","IL7","IRF1","IRF2","IRF7","IRF9","ISG15","ISG20","LAMP3","LAP3","LGALS3BP","LPAR6","LY6E","MOV10","MX1","NCOA7","NMI","NUB1","OAS1","OASL","OGFR","PARP12","PARP14","PARP9","PLSCR1","PNPT1","HELZ2","PROCR","PSMA3","PSMB8","PSMB9","PSME1","PSME2","RIPK2","RNF31","RSAD2","RTP4","SAMD9","SAMD9L","SELL","SLC25A28","SP110","STAT2","TAP1","TDRD7","TMEM140","TRAFD1","TRIM14","TRIM21","TRIM25","TRIM26","TRIM5","TXNIP","UBA7","UBE2L6","USP18","WARS1")
IFNg<-c("ADAR","APOL6","ARID5B","ARL4A","AUTS2","B2M","BANK1","BATF2","BPGM","BST2","BTG1","C1R","C1S","CASP1","CASP3","CASP4","CASP7","CASP8","CCL2","CCL5","CCL7","CD274","CD38","CD40","CD69","CD74","CD86","CDKN1A","CFB","CFH","CIITA","CMKLR1","CMPK2","CSF2RB","CXCL10","CXCL11","CXCL9","DDX58","DDX60","DHX58","EIF2AK2","EIF4E3","EPSTI1","FAS","FCGR1A","FGL2","FPR1","CMTR1","GBP4","GBP6","GCH1","GPR18","GZMA","HERC6","HIF1A","HLA-A","HLA-B","HLA-DMA","HLA-DQA1","HLA-DRB1","HLA-G","ICAM1","IDO1","IFI27","IFI30","IFI35","IFI44","IFI44L","IFIH1","IFIT1","IFIT2","IFIT3","IFITM2","IFITM3","IFNAR2","IL10RA","IL15","IL15RA","IL18BP","IL2RB","IL4R","IL6","IL7","IRF1","IRF2","IRF4","IRF5","IRF7","IRF8","IRF9","ISG15","ISG20","ISOC1","ITGB7","JAK2","KLRK1","LAP3","LATS2","LCP2","LGALS3BP","LY6E","LYSMD2","MARCHF1","METTL7B","MT2A","MTHFD2","MVP","MX1","MX2","MYD88","NAMPT","NCOA3","NFKB1","NFKBIA","NLRC5","NMI","NOD1","NUP93","OAS2","OAS3","OASL","OGFR","P2RY14","PARP12","PARP14","PDE4B","PELI1","PFKP","PIM1","PLA2G4A","PLSCR1","PML","PNP","PNPT1","HELZ2","PSMA2","PSMA3","PSMB10","PSMB2","PSMB8","PSMB9","PSME1","PSME2","PTGS2","PTPN1","PTPN2","PTPN6","RAPGEF6","RBCK1","RIPK1","RIPK2","RNF213","RNF31","RSAD2","RTP4","SAMD9L","SAMHD1","SECTM1","SELP","SERPING1","SLAMF7","SLC25A28","SOCS1","SOCS3","SOD2","SP110","SPPL2A","SRI","SSPN","ST3GAL5","ST8SIA4","STAT1","STAT2","STAT3","STAT4","TAP1","TAPBP","TDRD7","TNFAIP2","TNFAIP3","TNFAIP6","TNFSF10","TOR1B","TRAFD1","TRIM14","TRIM21","TRIM25","TRIM26","TXNIP","UBE2L6","UPP1","USP18","VAMP5","VAMP8","VCAM1","WARS1","XAF1","XCL1","ZBP1","ZNFX1")

#list from http://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=GOBP_RESPONSE_TO_TYPE_I_INTERFERON
# GOBP_RESPONSE_TO_TYPE_I_INTERFERON GO:0034340 GO:0060337,GO:0071357
type1_IFN<-c("TANK","ADAR","IFITM3","IFITM2","CDC37","USP18","TREX1","TRIM6","DCST1","TTLL12","YTHDF3","SAMHD1","LSM14A","SETD2","TBK1","CNOT7","UBE2K","IFNE","STING1","IFI27","IFIT1","IFNA1","IFNA2","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA21","IFNAR1","IFNAR2","IFNB1","IFNW1","IRAK1","IRF3","IRF7","JAK1","USP27X","MIR21","MMP12","IFIT1B","MX1","MYD88","OAS1","OAS2","OAS3","YTHDF2","SHFL","METTL3","IFNK","MAVS","USP29","PTPN1","PTPN2","PTPN6","PTPN11","CACTIN","AZI2","SHMT2","SMPD1","SP100","STAT1","STAT2","TRAF3","TYK2","WNT5A","MUL1","ZBP1","TRIM56","NLRC5","IFITM1","FADD","CH25H","TRIM41","RNF185","ISG15","IKBKE","TBKBP1")

diff_exp_sig_sex<-diff_exp_all_sex[which(diff_exp_all_sex$p_val_adj < 0.005 & abs(diff_exp_all_sex$avg_log2FC) > 1),]

diff_exp_sig_sex[which(diff_exp_sig_sex$gene%in%IFNa),]
diff_exp_sig_sex[which(diff_exp_sig_sex$gene%in%IFNg),]
diff_exp_sig_sex[which(diff_exp_sig_sex$gene%in%type1_IFN),]

#################
## Look at some interesting markers
#################
diff_exp_all[which(diff_exp_all$gene%in%c("KLRG1", "B3GAT1", "CD69","ITGAE")),]
diff_exp_all[which(diff_exp_all$gene%in%c("LYZ", "MARCO", "MRC1","PTPRC")),]

keygenes<-unique(diff_exp_sig_sex[which(diff_exp_sig_sex$gene%in%type1_IFN),]$gene)
diff_exp_all[which(diff_exp_all$gene%in%keygenes),]

Idents(d10x) <- "Sex"

cell_types<-unique(as.character(d10x$CellType_refined))

all_plots<-lapply(1:length(keygenes), function(y){
  plots<-lapply(1:length(cell_types),function(x){
    p<-VlnPlot(subset(d10x, subset = CellType_refined == cell_types[x]) , features = keygenes[y], pt.size = 0, log=T)
    p<-if(length(grep(cell_types[x], diff_exp_sig_sex[which(diff_exp_sig_sex$gene==keygenes[y]),]$cell.1))!=0){
      p+theme(plot.background = element_rect(color = "black",size = 2)) +xlab("") + ylab("")+ theme(legend.position="none")}else{  
        p +xlab("") + ylab("")+ theme(legend.position="none")
      }
    p})
  plot_grid(plotlist = plots, ncol=1)})


label_blank<-lapply(1:length(cell_types), function(x){
  ggplot()+geom_blank()+theme_void()+ggtitle(cell_types[x])+ theme(plot.title = element_text(hjust = 0.5,vjust = -30))  })
label_blank<-plot_grid(plotlist = label_blank, ncol=1)

plot_grid(label_blank, plot_grid(plotlist=all_plots, ncol=length(keygenes)), rel_widths=c(0.1,1))
ggsave2(here("figures", "sex_IFN1_genes.pdf"), w=20,h=45)
ggsave2(here("figures/jpeg", "sex_IFN1_genes.jpeg"), w=20,h=20,bg="white")


# 
# # ### Top DE genes
# diff_exp_all %>%
#   group_by(cell.1) %>%
#   top_n(n = 10, wt = abs(avg_log2FC)) -> top10
# 
# top_DE<-as.data.frame(top10)
# 
# Idents(d10x) <- d10x$AgeGroup
# 
# label_blank<-lapply(1:length(cell_types), function(x){
#   ggplot()+geom_blank()+theme_void()+ggtitle(cell_types[x])+ theme(plot.title = element_text(hjust = 0.5,vjust = -30))  })
# label_blank<-plot_grid(plotlist = label_blank, ncol=1)
# 
# plot_list_top<-lapply(1:length(cell_types), function(x){
#   plots <- VlnPlot(subset(d10x, subset = CellType_rough == cell_types[x]) , features = top_DE[grep(cell_types[x],top_DE$cell.1),"gene"], pt.size = 0, log=T)
#   plots <- lapply(X = plots, FUN = function(p) p + fillscale_age +xlab("")+ theme(plot.title = element_text(size = 15)))
#   plot_grid(plotlist = plots, nrow=1)})
# top_DE_plot<-plot_grid(plotlist = plot_list_top, nrow=length(cell_types))
# 
# plot_grid(label_blank, top_DE_plot, rel_widths=c(0.1,1))
# 
# ggsave2(here("figures", "TopDE_adult_ped.pdf"), w=20,h=20)
# ggsave2(here("figures/jpeg", "TopDE_adult_ped.jpeg"), w=20,h=20,bg="white")
