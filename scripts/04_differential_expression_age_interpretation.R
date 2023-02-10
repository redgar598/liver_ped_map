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
library(scales)


options(stringsAsFactors = FALSE)

source("scripts/00_pretty_plots.R")



## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(file = here("data","d10x_adult_ped_raw.rds"))

######
## add cell type labels
######
load(here("data","adult_ped_cellRefined_withDropletQC.rds"))

cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
identical(colnames(d10x), cell_label$index)

d10x <- AddMetaData(d10x, metadata = cell_label)

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")

## testing factor
d10x$cell_age<-paste(d10x$CellType_refined, d10x$AgeGroup, sep = "_")
Idents(d10x) <- "cell_age"

table(d10x$CellType_refined, d10x$AgeGroup)

## testing factor
levels(d10x$CellType_refined)[which(levels(d10x$CellType_refined)=="LSEC\n(Hepatocyte Like)")]<-"LSEC_hep"
levels(d10x$CellType_refined)[which(levels(d10x$CellType_refined)=="LSEC")]<-"LSEC_nothep"
levels(d10x$CellType_refined)[which(levels(d10x$CellType_refined)=="Neutrophil\n(DEFA+)")]<-"Neutrophil_DEFA"
levels(d10x$CellType_refined)[which(levels(d10x$CellType_refined)=="Neutrophil")]<-"Neutrophil_notDEFA"


##########
## Load DE Results
##########
cell_types<-unique(as.character(d10x$CellType_refined))
cell_types<-cell_types[-grep("Hepatocyte Like",cell_types)]
#no neutrophils in peds
cell_types<-cell_types[-grep("Neutrophil",cell_types)]
cell_types<-cell_types[-grep("Low Quality",cell_types)]
cell_types

DE_monte_carlo<-do.call(rbind, lapply(cell_types, function(celltype){
  print(celltype)
  load(here("data",paste(celltype,"adult_ped_diff_motecarlo_1000.RData",sep="_")))
  DE_monte_carlo}))

DE_monte_carlo_sig<-DE_monte_carlo[which(DE_monte_carlo$monte_carlo_sig<0.001),]

table(DE_monte_carlo_sig$cell)
table(DE_monte_carlo_sig$gene)[which(table(DE_monte_carlo_sig$gene)>0)]
length(unique(DE_monte_carlo_sig$gene))
sig_genes<-table(DE_monte_carlo_sig$gene)[which(table(DE_monte_carlo_sig$gene)>0)]
ord_name<-names(sig_genes[rev(order(sig_genes))])

summary_tbl<-as.data.frame(DE_monte_carlo_sig %>%
  select(gene, cell) %>% 
  group_by(gene) %>%
  mutate(all_cells = paste(cell, collapse = " | "))%>%
  select(gene, all_cells))
summary_tbl<-summary_tbl[!duplicated(summary_tbl),]
summary_tbl$all_cells<-gsub("\n"," ", summary_tbl$all_cells)

write.csv(file=here("data","Significant_genes_adult_ped.csv"),summary_tbl[match(ord_name,summary_tbl$gene),])

############
## Look at hits
############
## Signatures
myeloid_immune_supressive<-c("CTSB","CD163","MS4A7","FOLR2","GPNMB","VSIG4","HMOX1","MSR1")
inflammatory_macs<-c("CD74","HLA-DRA","TYROBP","C1QC","HLA-DPA1","HLA-DPB1","LYZ","S100A6")
exhausted_tcells<-c("TOX","PDCD1","LAG3","TNFRSF9","CXCL13","ENTPD1","HAVCR2","CD38")


DE_monte_carlo_sig[which(DE_monte_carlo_sig$gene%in%c("KLRG1", "B3GAT1", "CD69","ITGAE")),]
DE_monte_carlo_sig[which(DE_monte_carlo_sig$gene%in%c("LYZ", "MARCO", "MRC1","PTPRC")),]

keygenes<-c("KLRG1", "B3GAT1", "CD69","ITGAE","LYZ", "MARCO", "MRC1","PTPRC")
DE_monte_carlo_sig[which(DE_monte_carlo_sig$gene%in%keygenes),]


##################
## IFN signalling
##################
#Lists from GSEA "HALLMARK_INTERFERON_GAMMA_RESPONSE" and "HALLMARK_INTERFERON_ALPHA_RESPONSE"
IFNa<-c("ADAR","B2M","BATF2","BST2","C1S","CASP1","CASP8","CCRL2","CD47","CD74","CMPK2","CNP","CSF1","CXCL10","CXCL11","DDX60","DHX58","EIF2AK2","ELF1","EPSTI1","MVB12A","TENT5A","CMTR1","GBP2","GBP4","GMPR","HERC6","HLA-C","IFI27","IFI30","IFI35","IFI44","IFI44L","IFIH1","IFIT2","IFIT3","IFITM1","IFITM2","IFITM3","IL15","IL4R","IL7","IRF1","IRF2","IRF7","IRF9","ISG15","ISG20","LAMP3","LAP3","LGALS3BP","LPAR6","LY6E","MOV10","MX1","NCOA7","NMI","NUB1","OAS1","OASL","OGFR","PARP12","PARP14","PARP9","PLSCR1","PNPT1","HELZ2","PROCR","PSMA3","PSMB8","PSMB9","PSME1","PSME2","RIPK2","RNF31","RSAD2","RTP4","SAMD9","SAMD9L","SELL","SLC25A28","SP110","STAT2","TAP1","TDRD7","TMEM140","TRAFD1","TRIM14","TRIM21","TRIM25","TRIM26","TRIM5","TXNIP","UBA7","UBE2L6","USP18","WARS1")
IFNg<-c("ADAR","APOL6","ARID5B","ARL4A","AUTS2","B2M","BANK1","BATF2","BPGM","BST2","BTG1","C1R","C1S","CASP1","CASP3","CASP4","CASP7","CASP8","CCL2","CCL5","CCL7","CD274","CD38","CD40","CD69","CD74","CD86","CDKN1A","CFB","CFH","CIITA","CMKLR1","CMPK2","CSF2RB","CXCL10","CXCL11","CXCL9","DDX58","DDX60","DHX58","EIF2AK2","EIF4E3","EPSTI1","FAS","FCGR1A","FGL2","FPR1","CMTR1","GBP4","GBP6","GCH1","GPR18","GZMA","HERC6","HIF1A","HLA-A","HLA-B","HLA-DMA","HLA-DQA1","HLA-DRB1","HLA-G","ICAM1","IDO1","IFI27","IFI30","IFI35","IFI44","IFI44L","IFIH1","IFIT1","IFIT2","IFIT3","IFITM2","IFITM3","IFNAR2","IL10RA","IL15","IL15RA","IL18BP","IL2RB","IL4R","IL6","IL7","IRF1","IRF2","IRF4","IRF5","IRF7","IRF8","IRF9","ISG15","ISG20","ISOC1","ITGB7","JAK2","KLRK1","LAP3","LATS2","LCP2","LGALS3BP","LY6E","LYSMD2","MARCHF1","METTL7B","MT2A","MTHFD2","MVP","MX1","MX2","MYD88","NAMPT","NCOA3","NFKB1","NFKBIA","NLRC5","NMI","NOD1","NUP93","OAS2","OAS3","OASL","OGFR","P2RY14","PARP12","PARP14","PDE4B","PELI1","PFKP","PIM1","PLA2G4A","PLSCR1","PML","PNP","PNPT1","HELZ2","PSMA2","PSMA3","PSMB10","PSMB2","PSMB8","PSMB9","PSME1","PSME2","PTGS2","PTPN1","PTPN2","PTPN6","RAPGEF6","RBCK1","RIPK1","RIPK2","RNF213","RNF31","RSAD2","RTP4","SAMD9L","SAMHD1","SECTM1","SELP","SERPING1","SLAMF7","SLC25A28","SOCS1","SOCS3","SOD2","SP110","SPPL2A","SRI","SSPN","ST3GAL5","ST8SIA4","STAT1","STAT2","STAT3","STAT4","TAP1","TAPBP","TDRD7","TNFAIP2","TNFAIP3","TNFAIP6","TNFSF10","TOR1B","TRAFD1","TRIM14","TRIM21","TRIM25","TRIM26","TXNIP","UBE2L6","UPP1","USP18","VAMP5","VAMP8","VCAM1","WARS1","XAF1","XCL1","ZBP1","ZNFX1")

#list from http://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=GOBP_RESPONSE_TO_TYPE_I_INTERFERON
# GOBP_RESPONSE_TO_TYPE_I_INTERFERON GO:0034340 GO:0060337,GO:0071357
type1_IFN<-c("TANK","ADAR","IFITM3","IFITM2","CDC37","USP18","TREX1","TRIM6","DCST1","TTLL12","YTHDF3","SAMHD1","LSM14A","SETD2","TBK1","CNOT7","UBE2K","IFNE","STING1","IFI27","IFIT1","IFNA1","IFNA2","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA21","IFNAR1","IFNAR2","IFNB1","IFNW1","IRAK1","IRF3","IRF7","JAK1","USP27X","MIR21","MMP12","IFIT1B","MX1","MYD88","OAS1","OAS2","OAS3","YTHDF2","SHFL","METTL3","IFNK","MAVS","USP29","PTPN1","PTPN2","PTPN6","PTPN11","CACTIN","AZI2","SHMT2","SMPD1","SP100","STAT1","STAT2","TRAF3","TYK2","WNT5A","MUL1","ZBP1","TRIM56","NLRC5","IFITM1","FADD","CH25H","TRIM41","RNF185","ISG15","IKBKE","TBKBP1")

DE_monte_carlo_sig[which(DE_monte_carlo_sig$gene%in%IFNa),]
DE_monte_carlo_sig[which(DE_monte_carlo_sig$gene%in%IFNg),]
DE_monte_carlo_sig[which(DE_monte_carlo_sig$gene%in%type1_IFN),]

summary_tbl[which(summary_tbl$gene%in%c(type1_IFN,IFNa,IFNg)),]

###########
## Plot violins
###########
Idents(d10x) <- "AgeGroup"

DE_monte_carlo_sig$cell<-as.factor(DE_monte_carlo_sig$cell)
levels(DE_monte_carlo_sig$cell)<-c("CD3+ T-cells","Cholangiocytes","gd T-cells","Hepatocytes","HSC","KC Like",
                                  "Mature B-cells", "NK-like cells",
                                   "Plasma cells","RR Myeloid") 

cell_types<-as.factor(cell_types)
levels(cell_types)<-c("CD3+ T-cells","Cholangiocytes","gd T-cells","Hepatocytes","HSC","KC Like",
                      "LSEC","Mature B-cells", "NK-like cells",
                      "Plasma cells","RR Myeloid") 

d10x$CellType_refined<-as.factor(d10x$CellType_refined)
levels(d10x$CellType_refined)<-c("Mature B-cells","Plasma cells", "CD3+ T-cells","gd T-cells","NK-like cells","Cholangiocytes",
                      "LSEC","LSEC\n(Hepatocyte Like)", "Myeloid cells", "Myeloid cells\n(Hepatocyte Like)", "RR Myeloid", "KC Like",
                      "Neutrophil","Neutrophil\n(DEFA+)","HSC","Hepatocytes", "Low Quality") 



plot_key_genes<-function(keygenes, label){
  all_plots<-lapply(1:length(keygenes), function(y){
    plots<-lapply(1:length(cell_types),function(x){
      p<-VlnPlot(subset(d10x, subset = CellType_refined == as.character(cell_types[x])) , features = keygenes[y], pt.size = 0, log=T)
      p<-if( (as.character(cell_types[x]) %in% as.character(DE_monte_carlo_sig[which(DE_monte_carlo_sig$gene==keygenes[y]),]$cell))  ){
        p+theme(plot.background = element_rect(color = "black",size = 2)) +fillscale_age +xlab("") + ylab("")+ theme(legend.position="none")}else{
          p+fillscale_age +xlab("") + ylab("")+ theme(legend.position="none")
        }
      p})
    plot_grid(plotlist = plots, ncol=1)})
  
  
  label_blank<-lapply(1:length(cell_types), function(x){
    ggplot()+geom_blank()+theme_void()+ggtitle(cell_types[x])+ theme(plot.title = element_text(hjust = 0.5,vjust = -30))  })
  label_blank<-plot_grid(plotlist = label_blank, ncol=1)
  
  plot_grid(label_blank, plot_grid(plotlist=all_plots, ncol=length(keygenes)), rel_widths=c(0.1,1))
  wid<-length(keygenes)*3
  ggsave2(paste0(here("figures/"), label, "_adult_ped.pdf"), w=wid,h=30)
  ggsave2(paste0(here("figures/jpeg/"),label, "_adult_ped.jpeg"), w=wid,h=30,bg="white")}


#### IFNg
IFNg_sig<-unique(as.character(DE_monte_carlo_sig[which(DE_monte_carlo_sig$gene%in%IFNg),]$gene))
plot_key_genes(IFNg_sig, "IFNg_sig")

#### IFNa
IFNa_sig<-unique(as.character(DE_monte_carlo_sig[which(DE_monte_carlo_sig$gene%in%IFNa),]$gene))
plot_key_genes(IFNa_sig, "IFNa_sig")


table(DE_monte_carlo_sig$gene)[order(table(DE_monte_carlo_sig$gene))]
plot_key_genes(c("S100A9","S100A8","CCL4","CCL3"), "DE_unique_to_KClike")

plot_key_genes(myeloid_immune_supressive, "myeloid_immune_supressive_montecarlo")
plot_key_genes(inflammatory_macs, "inflammatory_macs_montecarlo")
plot_key_genes(exhausted_tcells, "exhausted_tcells_montecarlo")





########
## interpret DGE
########
## testing factor
levels(d10x$CellType_refined)[which(levels(d10x$CellType_refined)=="LSEC\n(Hepatocyte Like)")]<-"LSEC_hep"
levels(d10x$CellType_refined)[which(levels(d10x$CellType_refined)=="LSEC")]<-"LSEC_nothep"
levels(d10x$CellType_refined)[which(levels(d10x$CellType_refined)=="Neutrophil\n(DEFA+)")]<-"Neutrophil_DEFA"
levels(d10x$CellType_refined)[which(levels(d10x$CellType_refined)=="Neutrophil")]<-"Neutrophil_notDEFA"

d10x$cell_age<-paste(d10x$CellType_refined, d10x$AgeGroup, sep = "_")
Idents(d10x) <- "cell_age"

cell_types<-unique(as.character(d10x$CellType_refined))
cell_types<-cell_types[-grep("Hepatocyte Like",cell_types)]
#no neutrophils in peds
cell_types<-cell_types[-grep("Neutrophil",cell_types)]
cell_types<-cell_types[-grep("Low Quality",cell_types)]
cell_types[grep("CD3",cell_types)]<-"CD3"

contrasts_celltype_age<-do.call(rbind,lapply(1:length(cell_types), function(x){
  combinations(n = 2, r = 2, v = d10x$cell_age[grep(cell_types[x],d10x$cell_age)], repeats.allowed = FALSE)}))
contrasts_celltype_age
nrow(contrasts_celltype_age)


### get fold change in whole cohort for pathway analysis
## run DE 

diff_exp_all<-lapply(1:nrow(contrasts_celltype_age), function(x){
  de<-FindMarkers(d10x, ident.1 = contrasts_celltype_age[x,1], ident.2 = contrasts_celltype_age[x,2], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
  print(paste(contrasts_celltype_age[x,1],"vs", contrasts_celltype_age[x,2],":", nrow(de), sep=" "))
  de$gene<-rownames(de)
  rownames(de)<-NULL
  de<-de[,c(6,1:5)]
  de$cell.1<-contrasts_celltype_age[x,1]
  de$cell.2<-contrasts_celltype_age[x,2]
  de})


diff_exp_all<-do.call(rbind, diff_exp_all)

save(file=here("data","diff_exp_all_fullcohort.RData"),diff_exp_all)

#########
## pathway adult versus ped in KC and RR
#########
source("scripts/00_GSEA_function.R")
GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

pathway_plt<-function(de){
  gene_list = de$avg_log2FC
  names(gene_list) = de$gene
  gene_list = sort(gene_list, decreasing = TRUE)
  gene_list = gene_list[!duplicated(names(gene_list))]
  
  res = GSEA(gene_list, GO_file, pval = 0.05)
  
  plt_path<-res$Results
  plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
  plt_path$Enrichment_Cell<-"Up-regulated in \nAdult"
  plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \nPediatric"
  
  plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))
  
  plt_path$direction_label<-as.factor(plt_path$Enrichment)
  levels(plt_path$direction_label)<-c(0.1,-0.1)
  plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))
  
  # top and bottom 15
  plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])
  
  ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
    theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
    geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
    geom_vline(xintercept = 0, color="grey40")+scale_fill_manual(values=c("#fd8d3c","#6baed6"))+ 
    guides(fill = guide_legend(override.aes = list(size=5)))}

RR_GSEA<-pathway_plt(diff_exp_all[which(diff_exp_all$cell.1=="RR Myeloid_Adult"),])
save_plts(RR_GSEA, "GSEA_adult_ped_recently_recruited", w=20,h=10)

KC_GSEA<-pathway_plt(diff_exp_all[which(diff_exp_all$cell.1=="KC Like_Adult"),])
save_plts(KC_GSEA, "GSEA_adult_ped_KClike", w=20,h=10)

NK_GSEA<-pathway_plt(diff_exp_all[which(diff_exp_all$cell.1=="NK-like cells_Adult"),])
save_plts(NK_GSEA, "GSEA_adult_ped_NKlike", w=20,h=10)

HSC_GSEA<-pathway_plt(diff_exp_all[which(diff_exp_all$cell.1=="HSC_Adult"),])
save_plts(HSC_GSEA, "GSEA_adult_ped_HSClike", w=20,h=10)


##############
## Volcano
##############
vol_celltype<-function(cellType){
  diff_exp_all_celltype<-diff_exp_all[grep(cellType, diff_exp_all$cell.1),]
  volcano<-data.frame(gene=diff_exp_all_celltype$gene,Pvalue=diff_exp_all_celltype$p_val, Delta_Beta=diff_exp_all_celltype$avg_log2FC)
  
  #Thresholds 
  dB<-1 #delta beta cutoff
  Pv<-1.8e-07 #Pvalue cutoff
  
  volcano<-volcano[complete.cases(volcano),]
  
  ## positive delta beta is hypomethylated (code for volcano should be right now, should colors change?)
  color3<-sapply(1:nrow(volcano), function(x) if(volcano$Pvalue[x]<=Pv){
    if(abs(volcano$Delta_Beta[x])>dB){
      if(volcano$Delta_Beta[x]>dB){"Higher Expression in Adults\n(with Potential Biological Impact)"}else{"Higher Expression in Pediatric\n(with Potential Biological Impact)"}
    }else{if(volcano$Delta_Beta[x]>0){"Higher Expression in Adults"}else{"Higher Expression in Pediatric"}}}else{"Not Significantly Different"})
  
  volcano$Interesting_CpG3<-color3
  
  
  # COLORS! define here so they are consistent between plots
  # so even if you don't have CpGs in a color catagory the pattern will be maintained
  myColors <- c(muted("red", l=80, c=30),"red",muted("blue", l=70, c=40),"blue", "grey")
  
  color_possibilities<-c("Higher Expression in Adults",
                         "Higher Expression in Adults\n(with Potential Biological Impact)",
                         "Higher Expression in Pediatric",
                         "Higher Expression in Pediatric\n(with Potential Biological Impact)",
                         "Not Significantly Different")
  
  names(myColors) <- color_possibilities
  colscale <- scale_color_manual(name = "Direction of Change",
                                 values = myColors, drop = FALSE)
  fillscale <- scale_fill_manual(name = "Direction of Change",
                                 values = myColors, drop = FALSE)
  
  
  sig_MC_celltype<-DE_monte_carlo_sig[grep(cellType,DE_monte_carlo_sig$cell),]
  volcano_label<-volcano[which(volcano$gene%in%sig_MC_celltype$gene),]
  volcano$sig<-"not_sig"
  volcano$sig[which(volcano$gene%in%sig_MC_celltype$gene)]<-"MC_sig"
  
  ggplot(volcano, aes(Delta_Beta, -log10(Pvalue), fill=Interesting_CpG3, color=sig))+
    geom_point(shape=21, size=1)+theme_bw()+
    fillscale+scale_color_manual(values=c("black","white"))+guides(color = "none")+
    geom_vline(xintercept=c(-dB,dB), color="grey60")+
    geom_hline(yintercept=-log10(Pv), color="grey60")+
    ylab("P Value (-log10)")+xlab("Differential Expression\n(Fold change)")+
    theme(plot.margin=unit(c(1,1,1,2),"cm"))+ 
    guides(fill = guide_legend(override.aes = list(size = 4)))+
    geom_text(aes(label=gene),volcano_label,color="black",vjust=4, hjust=1,size=3)
  
  ggsave(file=paste(here("figures/"),cellType,"_differential_Ped_adult.pdf", sep=""), h=8, w=10)
  ggsave(file=paste(here("figures/jpeg/"),cellType,"_differential_Ped_adult.jpeg", sep=""), h=8, w=10)}


vol_celltype("RR")
vol_celltype("KC")
vol_celltype("NK")
vol_celltype("HSC")

######
## FC correlation plot RR-KC
######
diff_exp_all_celltype<-merge(diff_exp_all[grep("RR", diff_exp_all$cell.1),], diff_exp_all[grep("KC", diff_exp_all$cell.1),], by="gene",  suffixes = c("_RR","_KC"))

sig_MC_RR<-DE_monte_carlo_sig[grep("RR",DE_monte_carlo_sig$cell),]
sig_MC_KC<-DE_monte_carlo_sig[grep("KC",DE_monte_carlo_sig$cell),]

diff_exp_all_celltype$sig<-sapply(1:nrow(diff_exp_all_celltype), function(x){
  if(diff_exp_all_celltype$gene[x]%in%sig_MC_RR$gene & diff_exp_all_celltype$gene[x]%in%sig_MC_KC$gene){"both"}else{
    if(diff_exp_all_celltype$gene[x]%in%sig_MC_RR$gene){"RR only"}else{
      if(diff_exp_all_celltype$gene[x]%in%sig_MC_KC$gene){"KC only"}else{"NS"}
    }
  }
})

diff_exp_all_celltype_label<-diff_exp_all_celltype[which(diff_exp_all_celltype$gene%in%c(sig_MC_KC$gene,sig_MC_RR$gene)),]


ggplot(diff_exp_all_celltype, aes(avg_log2FC_RR, avg_log2FC_KC, color=sig))+geom_point()+th_present+theme_bw()+
  ylab("KC Differential Expression\n(Fold change)")+xlab("RR Differential Expression\n(log2 Fold change)")+
  scale_color_manual(values = c("red","#f7057d","grey","#b80783"))+
  geom_text(aes(label=gene),diff_exp_all_celltype_label,color="black",vjust=-0.75, hjust=1,size=3)+
  geom_vline(xintercept = c(-1,1), color="grey")+  geom_hline(yintercept = c(-1,1), color="grey")+
  ylim(c(min(diff_exp_all$avg_log2FC),max(diff_exp_all$avg_log2FC)))+
  xlim(c(min(diff_exp_all$avg_log2FC),max(diff_exp_all$avg_log2FC)))+
  annotate("text", x=2, y=3.6, label="Higher Expression in Adults")+
  annotate("text", x=-3.5, y=-3.8, label="Higher Expression in Peds")
ggsave(file=here("figures/FC_correlation_KC_RR_differential_Ped_adult.pdf"), h=7, w=8)
ggsave(file=here("figures/jpeg/FC_correlation_KC_RR_differential_Ped_adult.jpeg"), h=7, w=8)




######
## FC correlation plot NK-KC
######
diff_exp_all_celltype<-merge(diff_exp_all[grep("NK", diff_exp_all$cell.1),], diff_exp_all[grep("KC", diff_exp_all$cell.1),], by="gene",  suffixes = c("_NK","_KC"))

sig_MC_NK<-DE_monte_carlo_sig[grep("NK",DE_monte_carlo_sig$cell),]
sig_MC_KC<-DE_monte_carlo_sig[grep("KC",DE_monte_carlo_sig$cell),]

diff_exp_all_celltype$sig<-sapply(1:nrow(diff_exp_all_celltype), function(x){
  if(diff_exp_all_celltype$gene[x]%in%sig_MC_NK$gene & diff_exp_all_celltype$gene[x]%in%sig_MC_KC$gene){"both"}else{
    if(diff_exp_all_celltype$gene[x]%in%sig_MC_NK$gene){"NK only"}else{
      if(diff_exp_all_celltype$gene[x]%in%sig_MC_KC$gene){"KC only"}else{"NS"}
    }
  }
})

diff_exp_all_celltype_label<-diff_exp_all_celltype[which(diff_exp_all_celltype$gene%in%c(sig_MC_KC$gene,sig_MC_NK$gene)),]

ggplot(diff_exp_all_celltype, aes(avg_log2FC_NK, avg_log2FC_KC, color=sig))+geom_point()+th_present+theme_bw()+
  ylab("KC Differential Expression\n(Fold change)")+xlab("NK Differential Expression\n(log2 Fold change)")+
  scale_color_manual(values = c("red","#f7057d","#b80783","grey"))+
  geom_text(aes(label=gene),diff_exp_all_celltype_label,color="black",vjust=-0.75, hjust=1,size=3)+
  geom_vline(xintercept = c(-1,1), color="grey")+  geom_hline(yintercept = c(-1,1), color="grey")+
  ylim(c(min(diff_exp_all$avg_log2FC),max(diff_exp_all$avg_log2FC)))+
  xlim(c(min(diff_exp_all$avg_log2FC),max(diff_exp_all$avg_log2FC)))+
  annotate("text", x=2, y=3.6, label="Higher Expression in Adults")+
  annotate("text", x=-3.5, y=-3.8, label="Higher Expression in Peds")
ggsave(file=here("figures/FC_correlation_KC_NK_differential_Ped_adult.pdf"), h=7, w=8)
ggsave(file=here("figures/jpeg/FC_correlation_KC_NK_differential_Ped_adult.jpeg"), h=7, w=8)

######
## FC correlation plot RR-NK
######
diff_exp_all_celltype<-merge(diff_exp_all[grep("NK", diff_exp_all$cell.1),], diff_exp_all[grep("RR", diff_exp_all$cell.1),], by="gene",  suffixes = c("_NK","_RR"))

sig_MC_NK<-DE_monte_carlo_sig[grep("NK",DE_monte_carlo_sig$cell),]
sig_MC_RR<-DE_monte_carlo_sig[grep("RR",DE_monte_carlo_sig$cell),]

diff_exp_all_celltype$sig<-sapply(1:nrow(diff_exp_all_celltype), function(x){
  if(diff_exp_all_celltype$gene[x]%in%sig_MC_NK$gene & diff_exp_all_celltype$gene[x]%in%sig_MC_RR$gene){"both"}else{
    if(diff_exp_all_celltype$gene[x]%in%sig_MC_NK$gene){"NK only"}else{
      if(diff_exp_all_celltype$gene[x]%in%sig_MC_RR$gene){"RR only"}else{"NS"}
    }
  }
})

diff_exp_all_celltype_label<-diff_exp_all_celltype[which(diff_exp_all_celltype$gene%in%c(sig_MC_RR$gene,sig_MC_NK$gene)),]

ggplot(diff_exp_all_celltype, aes(avg_log2FC_NK, avg_log2FC_RR, color=sig))+geom_point()+th_present+theme_bw()+
  ylab("RR Differential Expression\n(Fold change)")+xlab("NK Differential Expression\n(log2 Fold change)")+
  scale_color_manual(values = c("red","#f7057d","grey"))+
  geom_text(aes(label=gene),diff_exp_all_celltype_label,color="black",vjust=-0.75, hjust=1,size=3)+
  geom_vline(xintercept = c(-1,1), color="grey")+  geom_hline(yintercept = c(-1,1), color="grey")+
  ylim(c(min(diff_exp_all$avg_log2FC),max(diff_exp_all$avg_log2FC)))+
  xlim(c(min(diff_exp_all$avg_log2FC),max(diff_exp_all$avg_log2FC)))+
  annotate("text", x=2, y=3.6, label="Higher Expression in Adults")+
  annotate("text", x=-3.5, y=-3.8, label="Higher Expression in Peds")
ggsave(file=here("figures/FC_correlation_RR_NK_differential_Ped_adult.pdf"), h=7, w=8)
ggsave(file=here("figures/jpeg/FC_correlation_RR_NK_differential_Ped_adult.jpeg"), h=7, w=8)


